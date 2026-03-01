//! Module defining a PDB-statistics data struct
//! (currently [`GridScorer`]) that converts sequences
//! to feature grids.
use crate::{
    datatypes::{AAMap, Aminoacid, aa_canonical_str},
    leak_vec,
};
use bumpalo::{Bump, collections::Vec};
use std::{array, cmp};
const MIN_XMER: usize = 1;
/// Placeholder struct that uses PDB-statistics to
/// convert sequences to feature grids.
/// 
/// This struct is meant to replicate the `GridScore`
/// class defined in the standalone LLPhyScore executable
/// found in Julie Forman-Kay's lab's [python package],
/// written by Hao Cai (@haocai1992).
/// 
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
/// 
/// It is currently missing a lot of details obviously.
pub struct GridScorer {}
/// A mapping from residue-type to a list of
/// short-range and long-range statistics calculated at
/// each subsequence centered at an aminoacid of that type.
pub type FeatureGrid<'a> = AAMap<&'a [FeatureGridEntry]>;
/// Short-range and long-range grid statistics for one
/// subsequence.
#[derive(Clone)]
pub struct FeatureGridEntry {
    /// Short range statistic.
    pub sr: f64,
    /// Long range statistic.
    pub lr: f64,
}
/// A helper struct for [`GridScorer`].
/// 
/// This type exists as it is right now because of the `GridScore`
/// class's `ZGridDB`, which is a dictionary with keys being floats
/// rounded to `0.5`. It will probably change in the future, because
/// why would you index using floats?!
/// 
/// You can find the `GridScore` class in Julie Forman-Kay's lab's
/// [python package], written by Hao Cai (@haocai1992).
/// 
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
#[derive(Clone, Copy)]
struct LineKey(f64);
impl GridScorer {
    /// Turn a sequence into a biophysical feature grid (currently [`FeatureGrid`])
    /// by looking at the average value of some given biophysical feature on windows
    /// centered on each residue type.
    pub fn score_sequence<'a>(
        &self,
        sequence: &aa_canonical_str,
        arena: &'a Bump,
    ) -> FeatureGrid<'a> {
        let Some(n_sites) = sequence.len().checked_sub(2 * MIN_XMER) else {
            return AAMap([&[]; 20]);
        };
        let mut residue_counts = AAMap([0_usize; 20]);
        for aa in sequence {
            residue_counts[aa] += 1;
        }
        let mut grid = AAMap(array::from_fn(|aaindex| {
            let cap_bound = residue_counts.0[aaindex];
            Vec::with_capacity_in(cap_bound, arena)
        }));
        for i in 1..=n_sites {
            let aa = sequence[i];
            let mut grid_counts = [0.0, 0.0, 0.0];
            let subseq = self.get_seq_centered_at(sequence, i);
            debug_assert!(subseq.len() >= 2 * MIN_XMER + 1);
            let all_xmer_scores = self.score_subseq_smart(subseq);
            for (xmer, (freq_sr, freq_lr)) in (MIN_XMER..).zip(all_xmer_scores) {
                let (mean_sr, std_sr) = self.query_avg_sdev_db::<true>(aa, xmer);
                let (mean_lr, std_lr) = self.query_avg_sdev_db::<false>(aa, xmer);
                let zscore_sr = (freq_sr - mean_sr) / std_sr;
                let zscore_lr = (freq_lr - mean_lr) / std_lr;
                let sr_linekey = LineKey::new(zscore_sr);
                let lr_linekey = LineKey::new(zscore_lr);
                let [total_slot, sr_slot, lr_slot] = &mut grid_counts;
                let (total, sr, lr) = self
                    .query_zgrid_db_exact(aa, xmer, sr_linekey, lr_linekey)
                    .unwrap_or_else(|| {
                        self.query_zgrid_db_nearest(aa, xmer, sr_linekey, lr_linekey)
                    });
                *total_slot += total;
                *sr_slot += sr;
                *lr_slot += lr;
            }
            let mut freq_by_grid_sr = 0.0;
            let mut freq_by_grid_lr = 0.0;
            let [total, sr, lr] = grid_counts;
            if total > 0.0 {
                freq_by_grid_sr = sr / total;
                freq_by_grid_lr = lr / total;
            }
            grid[aa].push(FeatureGridEntry {
                sr: freq_by_grid_sr,
                lr: freq_by_grid_lr,
            });
        }
        AAMap(grid.0.map(|buf| leak_vec(buf) as &[_]))
    }
    /// Try and get a subsequence centered at the given `center` index,
    /// starting at `self.max_xmer()` at shrinking until it fits inside
    /// the sequence.
    fn get_seq_centered_at<'a>(
        &self,
        sequence: &'a aa_canonical_str,
        center: usize,
    ) -> &'a aa_canonical_str {
        let space_on_left = center;
        let space_on_right = sequence.len() - 1 - center;
        let min_space = cmp::min(space_on_left, space_on_right);
        let xmer = cmp::min(min_space, self.max_xmer());
        &sequence[center - xmer..center + xmer + 1]
    }
    /// The maximum residue separation that this struct has collected statistics for.
    fn max_xmer(&self) -> usize {
        todo!()
    }
    /// Get `(short-range, long-range)` statistics for the given subsequence.
    fn score_subseq_smart<'a>(&self, subseq: &aa_canonical_str) -> std::vec::Vec<(f64, f64)> {
        debug_assert!(subseq.len() % 2 == 1);
        let midpoint = subseq.len() / 2;
        debug_assert!(midpoint <= self.max_xmer());
        let mut frequency_sum_a = 0.0;
        let mut frequency_sum_b = 0.0;
        let mut denom_total_a = 0.0;
        let mut denom_total_b = 0.0;
        let mid_res_aa = subseq[midpoint];
        let cap = cmp::min(midpoint, self.max_xmer() - 1);
        let mut scores = std::vec::Vec::with_capacity(cap);
        for p in 0..midpoint {
            let c_term_position = midpoint + 1 + p;
            let (freq_a, std_a) =
                self.query_freq_pair_db::<true, true>(mid_res_aa, subseq[c_term_position], p);
            let (freq_b, std_b) =
                self.query_freq_pair_db::<true, false>(mid_res_aa, subseq[c_term_position], p);
            frequency_sum_a += freq_a / std_a;
            frequency_sum_b += freq_b / std_b;
            denom_total_a += 1.0 / std_a;
            denom_total_b += 1.0 / std_b;
            let n_term_position = midpoint - 1 - p;
            let (freq_a, std_a) =
                self.query_freq_pair_db::<false, true>(subseq[n_term_position], mid_res_aa, p);
            let (freq_b, std_b) =
                self.query_freq_pair_db::<false, false>(subseq[n_term_position], mid_res_aa, p);
            frequency_sum_a += freq_a / std_a;
            frequency_sum_b += freq_b / std_b;
            denom_total_a += 1.0 / std_a;
            denom_total_b += 1.0 / std_b;
            if p < cap {
                let final_frequency_a = frequency_sum_a / denom_total_a;
                let final_frequency_b = frequency_sum_b / denom_total_b;
                scores.push((final_frequency_a, final_frequency_b))
            }
        }
        scores
    }
    /// Placeholder method for getting a `(mean, std)` pair from the
    /// "pair-frequency-DB".
    /// 
    /// This method exists as it is right now because of the `GridScore`
    /// class's `PairFreqDB`, which is a dictionary that reports mean
    /// and standard deviation of a biophysical characteristic for pairs
    /// of residues at a given distance.
    /// 
    /// You can find the `GridScore` class in Julie Forman-Kay's lab's
    /// [python package], written by Hao Cai (@haocai1992).
    /// 
    /// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
    fn query_freq_pair_db<const X: bool, const SR: bool>(
        &self,
        aa_a: Aminoacid,
        aa_b: Aminoacid,
        distance: usize,
    ) -> (f64, f64) {
        todo!()
    }
    /// Placeholder method for getting a `(mean, std)` pair from the
    /// "avg-sdev-DB".
    /// 
    /// This method exists as it is right now because of the `GridScore`
    /// class's `AvgSdevDB`, which is a dictionary that reports mean
    /// and standard deviation of a biophysical characteristic for x-mers
    /// containing the given `aa` that are `xmer` long.
    /// 
    /// You can find the `GridScore` class in Julie Forman-Kay's lab's
    /// [python package], written by Hao Cai (@haocai1992).
    /// 
    /// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
    fn query_avg_sdev_db<const SR: bool>(&self, aa: Aminoacid, xmer: usize) -> (f64, f64) {
        todo!()
    }
    /// Placeholder method for getting a `(count, short-range, long-range)`
    /// data tuple from the "z-grid-DB". This method may fail to return
    /// such a tuple if the statistics for the two given [`LineKey`]s are
    /// not present in the database.
    /// 
    /// This method exists as it is right now because of the `GridScore`
    /// class's `ZGridDB`, which is a dictionary that reports the tuple
    /// I just described above (I don't know how to unpack the def more
    /// at this point, but this will all be replaced anyway!).
    /// 
    /// You can find the `GridScore` class in Julie Forman-Kay's lab's
    /// [python package], written by Hao Cai (@haocai1992).
    /// 
    /// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
    fn query_zgrid_db_exact(
        &self,
        aa: Aminoacid,
        xmer: usize,
        sr_linekey: LineKey,
        lr_linekey: LineKey,
    ) -> Option<(f64, f64, f64)> {
        todo!()
    }
    /// Placeholder method for getting a `(count, short-range, long-range)`
    /// data tuple from the "z-grid-DB". This method is expected to be
    /// run on [`LineKey`]s that were not present in the data struct via
    /// [`Self::query_zgrid_db_exact`].
    /// 
    /// This method exists as it is right now because of the `GridScore`
    /// class's `ZGridDB`, which is a dictionary that reports the tuple
    /// I just described above (I don't know how to unpack the def more
    /// at this point, but this will all be replaced anyway!).
    /// 
    /// You can find the `GridScore` class in Julie Forman-Kay's lab's
    /// [python package], written by Hao Cai (@haocai1992).
    /// 
    /// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
    fn query_zgrid_db_nearest(
        &self,
        aa: Aminoacid,
        xmer: usize,
        sr_linekey: LineKey,
        lr_linekey: LineKey,
    ) -> (f64, f64, f64) {
        todo!()
    }
}
impl LineKey {
    /// Make a new [`LineKey`].
    /// 
    /// Ripped from `make_linekey` in Julie Forman-Kay's lab's
    /// standalone LLPhyScore executable in their [python package],
    /// written by Hao Cai (@haocai1992).
    /// 
    /// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
    pub fn new(zscore: f64) -> Self {
        LineKey(((zscore / 2.0).round() * 2.0).clamp(-8.0, 12.0))
    }
}