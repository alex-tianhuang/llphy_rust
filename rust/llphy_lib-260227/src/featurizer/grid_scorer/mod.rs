//! Module defining [`GridScorer`] and [`GridScore`].
use crate::{
    datatypes::{AAIndex, AAMap, MAX_XMER}, featurizer::grid_scorer::z_grid_db::deserialize_z_grid_db, leak_vec
};
use anyhow::Error;
use borsh::BorshSerialize;
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
use borsh::BorshDeserialize;
#[cfg(feature = "simd")]
mod simd;
#[cfg(feature = "simd")]
pub use simd::{
    PairFreqEntrySum, ZGridDBEntry, ZGridEntrySum, ZGridSubtable, AvgSdevDBEntry, PairFreqDBEntry
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
        sequence: &[AAIndex],
        arena: &'a Bump,
    ) -> GridScore<'a> {
        let Some(n_sites) = sequence.len().checked_sub(2) else {
            return GridScore {
                feature_a_scores: AAMap([&[]; 20]),
                feature_b_scores: AAMap([&[]; 20]),
            };
        };
        let mut trimmed_residue_counts = AAMap::default();
        for &aa in &sequence[1..=n_sites] {
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
            for j in 0..num_windows {
                let xmer = unsafe { XmerSize::new_unchecked(j + 1) };
                let n_term_position = relative_midpoint - xmer.get();
                inner_accumulator +=
                    &self.pair_freqs[aa][j].n_terminal_mapping[subseq[n_term_position]];
                let c_term_position = relative_midpoint + xmer.get();
                inner_accumulator +=
                    &self.pair_freqs[aa][j].c_terminal_mapping[subseq[c_term_position]];
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
/// Moral equivalent of implementing deserialization on [`GridScorer`],
/// but with a memory arena to put dynamically allocated subtables into.
pub fn deserialize_grid_scorer<'a>(buf: &mut &[u8], arena: &'a Bump) -> Result<GridScorer<'a>, Error> {
    let pair_freqs = <PairFreqDB>::deserialize(buf)?;
    let avg_sdevs = <AvgSdevDB>::deserialize(buf)?;
    let z_grid = deserialize_z_grid_db(buf, arena)?;
    Ok(GridScorer { pair_freqs, avg_sdevs, z_grid })
}
impl BorshSerialize for GridScorer<'_> {
    fn serialize<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        self.pair_freqs.serialize(writer)?;
        self.avg_sdevs.serialize(writer)?;
        self.z_grid.serialize(writer)?;
        Ok(())
    }
}
/// Try and get a subsequence centered at the given `center` index,
/// starting with spans of `MAX_XMER` at shrinking until it fits
/// inside the sequence.
fn get_subseq_centered_at(sequence: &[AAIndex], center: usize) -> &[AAIndex] {
    let space_on_left = center;
    let space_on_right = sequence.len() - 1 - center;
    let min_space = cmp::min(space_on_left, space_on_right);
    let xmer = cmp::min(min_space, MAX_XMER);
    &sequence[center - xmer..center + xmer + 1]
}
