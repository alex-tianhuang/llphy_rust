//! Module defining [`GridScorer`] and [`GridScore`].
use crate::{
    datatypes::{AAMap, MAX_XMER, aa_canonical_str},
    leak_vec,
};
use bumpalo::{Bump, collections::Vec};
use std::cmp;
pub use xmer::{XmerIndexableArray, XmerSize, xmer_sizes};
mod avg_sdev_db;
pub use avg_sdev_db::AvgSdevDB;
mod xmer;
mod z_grid_db;
mod pair_freq_db;
pub use pair_freq_db::PairFreqDB;
pub use z_grid_db::ZGridDB;
#[cfg(feature = "simd")]
mod simd;
#[cfg(feature = "simd")]
pub use simd::{
    PairFreqEntrySum, ZGridDBEntry, ZGridEntrySum, ZGridSubtable,AvgSdevDBEntry, PairFreqDBEntry
};
#[cfg(not(feature = "simd"))]
mod no_simd;
#[cfg(not(feature = "simd"))]
pub use no_simd::{
    PairFreqEntrySum, ZGridDBEntry, ZGridEntrySum, ZGridSubtable, AvgSdevDBEntry, PairFreqDBEntry
};
/// A struct that contains all the necessary data to
/// make biophysical feature grids ([`GridScore`]s)
/// from sequences.
pub struct GridScorer<'a> {
    pub pair_freqs: PairFreqDB,
    pub avg_sdevs: AvgSdevDB,
    pub z_grid: ZGridDB<'a>,
}
/// A collection of biophysical feature scores for
/// each residue in a sequence.
pub struct GridScore<'a> {
    pub feature_a_scores: AAMap<&'a [f64]>,
    pub feature_b_scores: AAMap<&'a [f64]>,
}
impl GridScorer<'_> {
    /// Turn a sequence into a biophysical feature grid (currently [`GridScore`])
    /// by looking at the average value of some given biophysical feature on windows
    /// centered on each residue type.
    pub fn score_sequence<'a>(
        &self,
        sequence: &aa_canonical_str,
        arena: &'a Bump,
    ) -> GridScore<'a> {
        let Some(n_sites) = sequence.len().checked_sub(2) else {
            return GridScore {
                feature_a_scores: AAMap([&[]; 20]),
                feature_b_scores: AAMap([&[]; 20]),
            };
        };
        let mut trimmed_residue_counts = AAMap::default();
        for aa in &sequence[1..=n_sites] {
            trimmed_residue_counts[aa] += 1;
        }
        let mut feature_a_scores = AAMap(std::array::from_fn(|aaindex| {
            let cap = trimmed_residue_counts.0[aaindex];
            Vec::with_capacity_in(cap, arena)
        }));
        let mut feature_b_scores = AAMap(std::array::from_fn(|aaindex| {
            let cap = trimmed_residue_counts.0[aaindex];
            Vec::with_capacity_in(cap, arena)
        }));
        for i in 1..=n_sites {
            let aa = sequence[i];
            let subseq = get_subseq_centered_at(sequence, i);
            debug_assert!(subseq.len() >= 3);
            debug_assert!(subseq.len() % 2 == 1);
            let mut outer_accumulator = ZGridEntrySum::new_zeroed();
            let mut inner_accumulator = PairFreqEntrySum::new_zeroed();
            let relative_midpoint = subseq.len() / 2;
            let num_windows = cmp::min(relative_midpoint, MAX_XMER);
            for i in 0..num_windows {
                let xmer_int = i + 1;
                let xmer = unsafe { XmerSize::new_unchecked(xmer_int) };
                let n_term_position = relative_midpoint - xmer.get();
                inner_accumulator +=
                    &self.pair_freqs[aa][i].n_terminal_mapping[subseq[n_term_position]];
                let c_term_position = relative_midpoint + xmer.get();
                inner_accumulator +=
                    &self.pair_freqs[aa][i].c_terminal_mapping[subseq[c_term_position]];
                let freqs = inner_accumulator.as_frequencies();
                let zscores = self.avg_sdevs[aa][xmer].freqs_to_zscores(freqs);
                outer_accumulator += self.z_grid[aa][xmer].lookup(zscores);
            }
            let [freq_a, freq_b] = outer_accumulator.as_frequencies();
            feature_a_scores[aa].push(freq_a);
            feature_b_scores[aa].push(freq_b);
        }
        GridScore {
            feature_a_scores: AAMap(feature_a_scores.0.map(|v| leak_vec(v) as &[_])),
            feature_b_scores: AAMap(feature_b_scores.0.map(|v| leak_vec(v) as &[_])),
        }
    }
}

/// Try and get a subsequence centered at the given `center` index,
/// starting with spans of `MAX_XMER` at shrinking until it fits
/// inside the sequence.
fn get_subseq_centered_at(sequence: &aa_canonical_str, center: usize) -> &aa_canonical_str {
    let space_on_left = center;
    let space_on_right = sequence.len() - 1 - center;
    let min_space = cmp::min(space_on_left, space_on_right);
    let xmer = cmp::min(min_space, MAX_XMER);
    &sequence[center - xmer..center + xmer + 1]
}
