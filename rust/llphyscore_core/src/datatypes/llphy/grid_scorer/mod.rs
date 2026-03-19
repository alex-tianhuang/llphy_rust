//! Module defining [`GridScorer`] and [`GridScore`].
//!
//! Structs that turn a sequence into a residue-level feature grids.
use crate::datatypes::{AAIndex, AAMap, Aminoacid, MAX_XMER, aa_canonical_str};
use anyhow::Error;
pub use avg_sdev_db::AvgSdevDB;
use borsh::BorshSerialize;
use bumpalo::Bump;
use bumpalo::collections::Vec;
pub use pair_freq_db::PairFreqDB;
use std::{mem::MaybeUninit, ptr::addr_of_mut};
pub use xmer::{XmerIndexableArray, XmerSize};
pub use z_grid_db::ZGridDB;
mod avg_sdev_db;
mod pair_freq_db;
#[cfg(feature = "simd")]
mod simd;
mod xmer;
mod z_grid_db;
#[cfg(feature = "simd")]
pub use simd::{
    AvgSdevDBEntry, PairFreqDBEntry, PairFreqEntrySum, ZGridDBEntry, ZGridEntrySum, ZGridSubtable,
};
#[cfg(not(feature = "simd"))]
mod no_simd;
#[cfg(not(feature = "simd"))]
pub use no_simd::{
    AvgSdevDBEntry, PairFreqDBEntry, PairFreqEntrySum, ZGridDBEntry, ZGridEntrySum, ZGridSubtable,
};
/// A struct that contains all the necessary data to
/// make biophysical feature grids ([`GridScore`]s)
/// from sequences.
///
/// Dev note
/// --------
/// The computation logic is in the [`crate::featurizer`] module.
/// This struct is basically just a bunch of lookup tables.
#[derive(PartialEq, BorshSerialize)]
pub struct GridScorer<'a> {
    pub pair_freqs: PairFreqDB,
    pub avg_sdevs: AvgSdevDB,
    pub z_grid: ZGridDB<'a>,
}
/// A helper struct for [`GridScorer::score_sequence`]
/// to reuse memory.
pub struct GridScoringBuffer<'a> {
    sequence_buffer: Vec<'a, AAIndex>,
    sequence_transposed_buffer: AAMap<Vec<'a, usize>>,
    grid_score_buffer: AAMap<[Vec<'a, f64>; 2]>,
}
/// A collection of biophysical feature scores for
/// each residue in a sequence.
pub struct GridScore<'a> {
    pub feature_a_scores: AAMap<&'a [f64]>,
    pub feature_b_scores: AAMap<&'a [f64]>,
}
impl<'a> GridScorer<'a> {
    /// Moral equivalent of implementing deserialization on [`GridScorer`],
    /// but does it in a memory arena and returns a reference to it.
    ///
    /// Dev note
    /// --------
    /// Memory arena is for:
    /// 1. `GridScorer` is 1MB (huge) and causes segfaults if I put it on the stack.
    /// 2. `ZGridDB` has a dynamically allocated size but don't want to use global allocator for that
    pub fn deserialize(buf: &mut &[u8], arena: &'a Bump) -> Result<&'a Self, Error> {
        let grid_scorer = arena.alloc_with(<MaybeUninit<GridScorer>>::uninit);
        let target = grid_scorer.as_mut_ptr();
        // SAFETY: target is a valid memory addr for a `GridScorer`,
        //         therefore so is the calculated pointer to `pair_freqs`.
        let pair_freqs =
            unsafe { &mut *addr_of_mut!((*target).pair_freqs).cast::<MaybeUninit<PairFreqDB>>() };
        PairFreqDB::deserialize_into(pair_freqs, buf)?;
        // SAFETY: target is a valid memory addr for a `GridScorer`,
        //         therefore so is the calculated pointer to `avg_sdevs`.
        let avg_sdevs =
            unsafe { &mut *addr_of_mut!((*target).avg_sdevs).cast::<MaybeUninit<AvgSdevDB>>() };
        AvgSdevDB::deserialize_into(avg_sdevs, buf)?;
        // SAFETY: target is a valid memory addr for a `GridScorer`,
        //         therefore so is the calculated pointer to `z_grid`.
        let z_grid = unsafe { &mut *addr_of_mut!((*target).z_grid).cast::<MaybeUninit<ZGridDB>>() };
        ZGridDB::deserialize_into(z_grid, buf, arena)?;
        // SAFETY: just initialized all fields above
        Ok(unsafe { grid_scorer.assume_init_ref() })
    }
    /// Turn a sequence into a biophysical feature grid (currently [`GridScore`])
    /// by looking at the average value of some given biophysical feature on windows
    /// centered on each residue type.
    pub fn score_sequence<'b>(&self, sequence: &aa_canonical_str, buffer: &'b mut GridScoringBuffer<'_>) -> GridScore<'b> {
        let Some(n_sites) = sequence.len().checked_sub(2) else {
            return GridScore {
                feature_a_scores: AAMap([&[]; 20]),
                feature_b_scores: AAMap([&[]; 20]),
            };
        };
        buffer.sequence_buffer.clear();
        buffer.sequence_buffer.extend(sequence.into_iter().map(Aminoacid::to_aaindex));
        let sequence = &*buffer.sequence_buffer;
        let mut site_counts = <AAMap<usize>>::default();
        for &aa in &sequence[1..=n_sites] {
            site_counts[aa] += 1;
        }
        let site_indexes = &mut buffer.sequence_transposed_buffer;
        for (buffer, cap) in site_indexes.0.iter_mut().zip(site_counts.0) {
            buffer.clear();
            buffer.reserve(cap);
        }
        for (i, &aa_i) in sequence[1..=n_sites].iter().enumerate() {
            site_indexes[aa_i].push(i + 1);
        }
        let mut feature_a_scores = AAMap::<&[f64]>([&[]; 20]);
        let mut feature_b_scores = AAMap::<&[f64]>([&[]; 20]);
        for (aa_i, ([feature_a_scores_i, feature_b_scores_i], sites_containing_aa_i)) in buffer.grid_score_buffer.0.iter_mut().zip(site_indexes.0.iter_mut()).enumerate() {
            let aa_i = unsafe { AAIndex::from_byte_unchecked(aa_i as u8) };
            let n_sites_i = sites_containing_aa_i.len();
            feature_a_scores_i.clear();
            feature_a_scores_i.reserve(n_sites_i);
            feature_b_scores_i.clear();
            feature_b_scores_i.reserve(n_sites_i);
            for &site in sites_containing_aa_i.iter() {
                let wingspan = wingspan_of(site, sequence.len());
                let mut outer_accumulator = ZGridEntrySum::new_zeroed();
                let mut inner_accumulator = PairFreqEntrySum::new_zeroed();
                let num_windows = std::cmp::min(wingspan, MAX_XMER);
                for j in 0..num_windows {
                    let xmer = unsafe { XmerSize::new_unchecked(j + 1) };
                    let n_term_position = site - xmer.get();
                    inner_accumulator +=
                        &self.pair_freqs[aa_i][j].n_terminal_mapping[sequence[n_term_position]];
                    let c_term_position = site + xmer.get();
                    inner_accumulator +=
                        &self.pair_freqs[aa_i][j].c_terminal_mapping[sequence[c_term_position]];
                    let freqs = inner_accumulator.as_frequencies();
                    let zscores = self.avg_sdevs[aa_i][xmer].freqs_to_zscores(freqs);
                    outer_accumulator += self.z_grid[aa_i][xmer].lookup(zscores);
                }
                let [freq_a, freq_b] = outer_accumulator.as_frequencies();
                feature_a_scores_i.push(freq_a);
                feature_b_scores_i.push(freq_b);
            }
            feature_a_scores[aa_i] = feature_a_scores_i;
            feature_b_scores[aa_i] = feature_b_scores_i;
        }
        GridScore {
            feature_a_scores,
            feature_b_scores,
        }
    }
}
impl<'a> GridScoringBuffer<'a> {
    /// Make a new [`GridScoringBuffer`].
    pub fn new(arena: &'a Bump) -> Self {
        Self { sequence_buffer: Vec::new_in(arena), sequence_transposed_buffer: AAMap(std::array::from_fn(|_| Vec::new_in(arena))), grid_score_buffer: AAMap(std::array::from_fn(|_| [Vec::new_in(arena), Vec::new_in(arena)])) }
    }
}

/// The maximum number `n` such that a residue at
/// `center` has `n` residues on either side of it
/// in a sequence of given length.
fn wingspan_of(center: usize, sequence_len: usize) -> usize {
    let space_on_left = center;
    let space_on_right = sequence_len - 1 - center;
    std::cmp::min(space_on_left, space_on_right)
}
