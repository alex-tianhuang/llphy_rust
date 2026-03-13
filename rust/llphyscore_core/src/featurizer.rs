//! Module defining [`Featurizer::featurize`], a function which
//! computes the biophysical features of many sequences.
//!
//! The featurization process is organized into constructing
//! a `GridScore`, and decoding that grid score using a [`GridDecoder`].
//!
//! These two steps are based on the `seqs2grids` + `get_g2w_scores`
//! functions in the standalone LLPhyScore executable
//! found in Julie Forman-Kay's lab's [python package],
//! written by Hao Cai (@haocai1992).
//!
//! [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore

use crate::{
    datatypes::{
        Aminoacid, FEATURE_NAMES, FeatureMatrix, GridDecoderPair, ModelTrainingBase, aa_canonical_str,
    },
    load_pkg_data::{load_grid_decoders, load_grid_scorer},
    utils::{lookup_by_pair, lookup_f2p, lookup_findex},
};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use std::fmt::Write;
use std::path::Path;
/// A struct helping compute the biophysical features of many sequences.
///
/// This struct is just my builder pattern struct so that I can lump all the
/// configurable elements of a computation (does the user want a progress
/// bar, does it need to check for keyboard interrupts, etc.) into an
/// `&self` and the signature of [`Featurizer::featurize`] can remain stable.
pub struct Featurizer<'a> {
    pkg_data_root: &'a Path,
    grid_decoders: &'a [(&'a str, GridDecoderPair)],
    /// Closure that checks for interrupts in the middle of [`Featurizer::featurize`],
    /// for example for python code where keyboard interrupts do not immediately
    /// abort the thread.
    interrupter: Option<&'a dyn Fn() -> Result<(), Error>>,
    /// Flag to turn on/off progress bar during [`Featurizer::featurize`].
    with_pbar: bool,
}
impl<'a> Featurizer<'a> {
    /// Load a featurizer with no interrupter and pbar disabled (default)
    /// using only the package data needed for [`load_grid_decoders`].
    pub fn load_new(
        pkg_data_root: &'a Path,
        model_train_base: ModelTrainingBase,
        arena: &'a Bump,
    ) -> Result<Self, Error> {
        let grid_decoders = load_grid_decoders(pkg_data_root, model_train_base, arena)?;
        Ok(Self::from_parts(pkg_data_root, grid_decoders))
    }
    /// Make a featurizer with no interrupter and pbar disabled (default)
    /// using a `pkg_data_root` and `grid_decoders`.
    pub fn from_parts(
        pkg_data_root: &'a Path,
        grid_decoders: &'a [(&'a str, GridDecoderPair)],
    ) -> Self {
        Self {
            pkg_data_root,
            grid_decoders,
            interrupter: None,
            with_pbar: false,
        }
    }
    /// Add a `pbar` to this featurizer if requested.
    pub fn with_pbar(mut self, with_pbar: bool) -> Self {
        self.with_pbar = with_pbar;
        self
    }
    /// Check for interruption signals while running [`Self::featurize`].
    /// 
    /// Currently only used to check for KeyboardInterrupts from python.
    pub fn with_interrupter(mut self, interrupter: &'a dyn Fn() -> Result<(), Error>) -> Self {
        self.interrupter = Some(interrupter);
        self
    }
    /// Given sequences and PDB statistics,
    /// compute all biophysical sequence features + a `feature_sum`.
    pub fn featurize<'b, 'c>(
        &self,
        sequences: impl ExactSizeIterator<Item = &'c aa_canonical_str> + Clone,
        arena: &'b Bump,
    ) -> Result<FeatureMatrix<'b, i64>, Error> {
        let num_features = FEATURE_NAMES.len();
        let num_sequences = sequences.len();
        let Some(max_seq_len) = sequences.clone().map(|seq| seq.len()).max() else {
            return Ok(FeatureMatrix {
                feature_names: &FEATURE_NAMES,
                data: &mut [],
            });
        };
        let row_size = num_features + 1;
        const UNINIT_SENTINEL: i64 = i64::MIN;
        let data = arena.alloc_slice_fill_copy(num_sequences * row_size, UNINIT_SENTINEL);
        let completed = arena.alloc_slice_fill_copy(FEATURE_NAMES.len(), false);
        let mut temp_arena = Bump::new();
        for i in 0..num_features {
            if completed[i] {
                continue;
            }
            temp_arena.reset();
            let feature_name = FEATURE_NAMES[i];
            let (pair_name, feature_names_for_pair) = lookup_f2p(feature_name).unwrap();
            let decoder_pair = lookup_by_pair(self.grid_decoders, pair_name).ok_or_else(|| {
                Error::msg(format!("unknown/unsupported feature name {}", feature_name))
            })?;
            let [feature_idx_a, feature_idx_b] = feature_names_for_pair.map(lookup_findex);
            let mut sequence_buffer = Vec::with_capacity_in(max_seq_len, &temp_arena);
            let grid_scorer = load_grid_scorer(self.pkg_data_root, pair_name, &temp_arena)?;
            let pbar = self.with_pbar.then(|| pbar(sequences.len()));
            pbar.as_ref().map(|pbar| {
                pbar.println(
                    bumpalo::format!(in &temp_arena, "COMPUTING FEATURE PAIR {}", pair_name),
                )
            });
            for (j, seq) in sequences.clone().enumerate() {
                self.interrupter.map(|f| f()).unwrap_or(Ok(()))?;
                sequence_buffer.clear();
                sequence_buffer.extend(seq.into_iter().map(Aminoacid::to_aaindex));
                let scored = grid_scorer.score_sequence(&sequence_buffer, &temp_arena);
                if let Some(slot) = feature_idx_a {
                    let slot = j * row_size + slot;
                    let value = decoder_pair.decoder_a.decode(&scored.feature_a_scores);
                    debug_assert_ne!(value, UNINIT_SENTINEL);
                    data[slot] = value;
                }
                if let Some(slot) = feature_idx_b {
                    let slot = j * row_size + slot;
                    let value = decoder_pair.decoder_b.decode(&scored.feature_b_scores);
                    debug_assert_ne!(value, UNINIT_SENTINEL);
                    data[slot] = value;
                }
                pbar.as_ref().map(|pbar| pbar.inc(1));
            }
            pbar.as_ref().map(ProgressBar::abandon);
            if let Some(i) = feature_idx_a {
                completed[i] = true;
            }
            if let Some(i) = feature_idx_b {
                completed[i] = true;
            }
        }
        debug_assert!(completed.iter().all(|b| *b));
        for feature_vec in data.chunks_exact_mut(row_size) {
            let [subfeatures @ .., sumfeatures] = feature_vec else {
                unreachable!()
            };
            *sumfeatures = subfeatures.iter().sum::<i64>();
        }
        debug_assert!(!data.iter().any(|x| *x == UNINIT_SENTINEL));
        Ok(FeatureMatrix {
            feature_names: &FEATURE_NAMES,
            data,
        })
    }
}

/// Pbar in the tqdm style.
pub fn pbar(n: usize) -> ProgressBar {
    ProgressBar::new(n as u64).with_style(
        ProgressStyle::with_template(
            "{percent}%[{wide_bar}] {pos}/{len} [{elapsed}<{eta}, {rate}]",
        )
        .unwrap()
        .with_key("rate", |state: &ProgressState, w: &mut dyn Write| {
            write!(w, "{:.2}it/s", state.per_sec()).unwrap();
        }),
    )
}
