//! Module defining [`featurize`], a function which
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
    datatypes::{FastaEntry, FeatureMatrix, find_pair_and_features_from_one_feature_name},
    load_pkg_data::load_grid_scorer,
};
use anyhow::Error;
use bumpalo::Bump;
pub use grid_decoder::GridDecoder;
pub use grid_scorer::{GridScore, GridScorer};
use indicatif::ProgressBar;
use pyo3::Python;
pub use thresholds::Thresholds;
mod grid_decoder;
pub(crate) mod grid_scorer;
mod thresholds;

/// Given sequences and PDB statistics,
/// compute all biophysical sequence features + a `feature_sum`.
///
/// IO
/// --
/// If `PBAR` is true, this function
/// prints the grid that is currently being worked on along with a
/// progress bar to whatever the default terminal for [`indicatif`] is.
pub fn featurize<'a, const PBAR: bool>(
    sequences: &[FastaEntry<'_>],
    feature_names: &'a [&'a str],
    grid_decoders: &[(&str, [GridDecoder; 2])],
    arena: &'a Bump,
    py: Python,
) -> Result<FeatureMatrix<'a, i64>, Error> {
    let num_features = feature_names.len();
    let num_sequences = sequences.len();
    let row_size = num_features + 1;
    const UNINIT_SENTINEL: i64 = i64::MIN;
    let data = arena.alloc_slice_fill_copy(num_sequences * row_size, UNINIT_SENTINEL);
    let completed = arena.alloc_slice_fill_copy(feature_names.len(), false);
    let mut per_feature_pair_arena = Bump::new();
    for i in 0..num_features {
        if completed[i] {
            continue;
        }
        per_feature_pair_arena.reset();
        let feature_name = feature_names[i];
        let (pair_name, feature_names) = find_pair_and_features_from_one_feature_name(feature_name)
            .ok_or_else(|| {
                Error::msg(format!("unknown/unsupported feature name {}", feature_name))
            })?;
        let (_, [grid_decoder_a, grid_decoder_b]) = grid_decoders
            .iter()
            .find(|(n, _)| *n == pair_name)
            .ok_or_else(|| {
                Error::msg(format!("unknown/unsupported feature name {}", feature_name))
            })?;
        let [feature_idx_a, feature_idx_b] = feature_names.map(|feature_name| {
            let location = feature_names.iter().find(|n| **n == feature_name)?;
            let idx = (location as *const &str as usize - feature_names.as_ptr() as usize)
                / size_of::<&str>();
            Some(idx)
        });
        // `GridScorer` is a 1MiB struct.
        // This is me trying my best not to put it on the stack.
        let grid_scorer = per_feature_pair_arena
            .alloc_try_with(|| load_grid_scorer(pair_name, &per_feature_pair_arena))?;
        if PBAR {
            let pbar = ProgressBar::new(sequences.len() as _);
            pbar.println(bumpalo::format!(in &per_feature_pair_arena, "COMPUTING FEATURE PAIR {}", pair_name));
            for (j, entry) in sequences.iter().enumerate() {
                let scored = grid_scorer.score_sequence(
                    entry.sequence,
                    feature_idx_a.is_some(),
                    feature_idx_b.is_some(),
                    &per_feature_pair_arena,
                );
                py.check_signals()?;
                if let Some((slot, grid_score)) =
                    feature_idx_a.zip(scored.feature_a_scores.as_ref())
                {
                    let slot = j * row_size + slot;
                    data[slot] = grid_decoder_a.compute(grid_score);
                }
                if let Some((slot, grid_score)) =
                    feature_idx_b.zip(scored.feature_b_scores.as_ref())
                {
                    let slot = j * row_size + slot;
                    data[slot] = grid_decoder_b.compute(grid_score);
                }
                pbar.inc(1);
            }
            pbar.abandon();
        } else {
            for (j, entry) in sequences.iter().enumerate() {
                let scored = grid_scorer.score_sequence(
                    entry.sequence,
                    feature_idx_a.is_some(),
                    feature_idx_b.is_some(),
                    &per_feature_pair_arena,
                );
                py.check_signals()?;
                if let Some((slot, grid_score)) =
                    feature_idx_a.zip(scored.feature_a_scores.as_ref())
                {
                    let slot = j * row_size + slot;
                    let value = grid_decoder_a.compute(grid_score);
                    debug_assert_ne!(value, UNINIT_SENTINEL);
                    data[slot] = value;
                }
                if let Some((slot, grid_score)) =
                    feature_idx_b.zip(scored.feature_b_scores.as_ref())
                {
                    let slot = j * row_size + slot;
                    let value = grid_decoder_b.compute(grid_score);
                    debug_assert_ne!(value, UNINIT_SENTINEL);
                    data[slot] = value;
                }
            }
        }
        if let Some(i) = feature_idx_a {
            completed[i] = true;
        }
        if let Some(i) = feature_idx_b {
            completed[i] = true;
        }
    }
    for feature_vec in data.chunks_exact_mut(row_size) {
        let [subfeatures @ .., sumfeatures] = feature_vec else {
            unreachable!()
        };
        *sumfeatures = subfeatures.iter().sum::<i64>();
    }
    debug_assert!(!data.iter().any(|x| *x == UNINIT_SENTINEL));
    Ok(FeatureMatrix {
        feature_names,
        data,
    })
}
