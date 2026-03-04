//! Module defining [`GridScorer`] and [`GridScore`].
use std::cmp;

pub use avg_sdev_db::AvgSdevDB;
use bumpalo::{Bump, collections::Vec};
pub use pair_freq_db::PairFreqDB;
pub use xmer::{XmerIndexableArray, XmerSize, xmer_sizes};
pub use z_grid_db::ZGridDB;

use crate::{
    datatypes::{
        AAMap, MAX_XMER, aa_canonical_str,
        llphy::features::grid_scorer::{
            avg_sdev_db::AvgSdevDBEntry, pair_freq_db::PairFreqDBEntry, z_grid_db::ZGridDBEntry,
        },
    },
    leak_vec,
};
mod avg_sdev_db;
mod pair_freq_db;
mod xmer;
mod z_grid_db;
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
///
/// Depending on which feature is to be calculated
/// there will/won't be a `Some` variant option.
pub struct GridScore<'a> {
    pub feature_a_scores: Option<AAMap<&'a [f64]>>,
    pub feature_b_scores: Option<AAMap<&'a [f64]>>,
}
impl GridScorer<'_> {
    /// Turn a sequence into a biophysical feature grid (currently [`GridScore`])
    /// by looking at the average value of some given biophysical feature on windows
    /// centered on each residue type.
    pub fn score_sequence<'a>(
        &self,
        sequence: &aa_canonical_str,
        include_a: bool,
        include_b: bool,
        arena: &'a Bump,
    ) -> GridScore<'a> {
        let Some(n_sites) = sequence.len().checked_sub(2) else {
            return GridScore {
                feature_a_scores: include_a.then(|| AAMap([&[] as _; 20])),
                feature_b_scores: include_b.then(|| AAMap([&[] as _; 20])),
            };
        };
        let mut trimmed_residue_counts = AAMap::default();
        for aa in &sequence[1..=n_sites] {
            trimmed_residue_counts[aa] += 1;
        }
        let mut feature_a_scores = include_a.then(|| {
            AAMap(std::array::from_fn(|aaindex| {
                let cap = trimmed_residue_counts.0[aaindex];
                Vec::with_capacity_in(cap, arena)
            }))
        });
        let mut feature_b_scores = include_b.then(|| {
            AAMap(std::array::from_fn(|aaindex| {
                let cap = trimmed_residue_counts.0[aaindex];
                Vec::with_capacity_in(cap, arena)
            }))
        });
        for i in 1..=n_sites {
            let aa = sequence[i];
            let subseq = get_subseq_centered_at(sequence, i);
            debug_assert!(subseq.len() >= 3);
            let all_xmer_scores = self.score_subseq_smart(subseq);
            let avg_sdevs = &self.avg_sdevs[aa];
            let mut outer_total = 0;
            let mut outer_weight_a = 0;
            let mut outer_weight_b = 0;
            for (xmer, (freq_a, freq_b)) in xmer_sizes().zip(all_xmer_scores) {
                let AvgSdevDBEntry {
                    avg_a,
                    std_a,
                    avg_b,
                    std_b,
                } = avg_sdevs[xmer];
                let zscore_a = (freq_a - avg_a) / std_a;
                let zscore_b = (freq_b - avg_b) / std_b;
                let ZGridDBEntry {
                    weight_total,
                    weight_a,
                    weight_b,
                } = *self.z_grid[aa][xmer].lookup(zscore_a, zscore_b);
                outer_total += weight_total;
                outer_weight_a += weight_a;
                outer_weight_b += weight_b;
            }
            if let Some(g) = feature_a_scores.as_mut() {
                let outer_freq_a = if outer_total == 0 {
                    0.0
                } else {
                    outer_weight_a as f64 / outer_total as f64
                };
                g[aa].push(outer_freq_a)
            }
            if let Some(g) = feature_b_scores.as_mut() {
                let outer_freq_b = if outer_total == 0 {
                    0.0
                } else {
                    outer_weight_b as f64 / outer_total as f64
                };
                g[aa].push(outer_freq_b)
            }
        }
        GridScore {
            feature_a_scores: feature_a_scores.map(|m| AAMap(m.0.map(|v| leak_vec(v) as &[_]))),
            feature_b_scores: feature_b_scores.map(|m| AAMap(m.0.map(|v| leak_vec(v) as &[_]))),
        }
    }
    
    /// Get feature `a` and `b` statistics for the given subsequence.
    fn score_subseq_smart<'a>(
        &self,
        subseq: &aa_canonical_str,
    ) -> impl ExactSizeIterator<Item = (f64, f64)> {
        debug_assert!(subseq.len() % 2 == 1);
        let midpoint = subseq.len() / 2;
        debug_assert!(midpoint <= MAX_XMER);
        let mut weight_sum_a = 0.0;
        let mut weight_sum_b = 0.0;
        let mut denom_total_a = 0.0;
        let mut denom_total_b = 0.0;
        let aa = subseq[midpoint];
        let cap = cmp::min(midpoint, MAX_XMER);
        let pair_freqs = &self.pair_freqs[aa];
        (0..cap).map(move |p| {
            let xmer = p + 1;
            let pair_freqs_subtable = &pair_freqs[p];
            let n_term_position = midpoint - xmer;
            let PairFreqDBEntry {
                weight_a,
                weight_b,
                total_a,
                total_b,
            } = pair_freqs_subtable.n_terminal_mapping[subseq[n_term_position]];
            weight_sum_a += weight_a;
            weight_sum_b += weight_b;
            denom_total_a += total_a;
            denom_total_b += total_b;
            let c_term_position = midpoint + xmer;
            let PairFreqDBEntry {
                weight_a,
                weight_b,
                total_a,
                total_b,
            } = pair_freqs_subtable.c_terminal_mapping[subseq[c_term_position]];
            weight_sum_a += weight_a;
            weight_sum_b += weight_b;
            denom_total_a += total_a;
            denom_total_b += total_b;
            let final_frequency_a = weight_sum_a / denom_total_a;
            let final_frequency_b = weight_sum_b / denom_total_b;
            (final_frequency_a, final_frequency_b)
        })
    }
}

/// Try and get a subsequence centered at the given `center` index,
/// starting with spans of `MAX_XMER` at shrinking until it fits
/// inside the sequence.
fn get_subseq_centered_at(
    sequence: &aa_canonical_str,
    center: usize,
) -> &aa_canonical_str {
    let space_on_left = center;
    let space_on_right = sequence.len() - 1 - center;
    let min_space = cmp::min(space_on_left, space_on_right);
    let xmer = cmp::min(min_space, MAX_XMER);
    &sequence[center - xmer..center + xmer + 1]
}