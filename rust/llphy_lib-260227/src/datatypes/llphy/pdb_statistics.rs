//! Module defining a PDB-statistics data struct
//! (currently [`GridScorer`]) that converts sequences
//! to feature grids.
use crate::{
    datatypes::{AAMap, Aminoacid, aa_canonical_str},
    leak_vec,
};
use bumpalo::{Bump, collections::Vec};
use std::{array, cmp::{self, Ordering}};
/// The minimum residue separation that we have collected residue statistics for.
const MIN_XMER: usize = 1;
/// The maximum residue separation that we have collected residue statistics for.
const MAX_XMER: usize = 40;
/// Struct that uses PDB-statistics to
/// convert sequences to feature grids.
/// 
/// This struct is meant to replicate the `GridScore`
/// class defined in the standalone LLPhyScore executable
/// found in Julie Forman-Kay's lab's [python package],
/// written by Hao Cai (@haocai1992).
/// 
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
pub struct GridScorer<'a> {
    pair_freq_db_x: &'a[AAMap<AAMap<[(f64, f64); 2]>>],
    pair_freq_db_y: &'a[AAMap<AAMap<[(f64, f64); 2]>>],
    avg_sdev_db: AAMap<&'a[[(f64, f64); 2]]>,
    // The two inner `LineKey`-indexed lists are sorted by linekey.
    z_grid_db: AAMap<&'a[&'a[(LineKey, &'a [(LineKey, (f64, f64, f64))])]]>,
}
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
#[derive(Clone, Copy, PartialOrd, PartialEq)]
struct LineKey(f64);
impl GridScorer<'_> {
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
                        self.query_zgrid_db_nearest(aa, xmer, zscore_sr, zscore_lr)
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
    /// starting at `MAX_XMER` at shrinking until it fits inside
    /// the sequence.
    fn get_seq_centered_at<'a>(
        &self,
        sequence: &'a aa_canonical_str,
        center: usize,
    ) -> &'a aa_canonical_str {
        let space_on_left = center;
        let space_on_right = sequence.len() - 1 - center;
        let min_space = cmp::min(space_on_left, space_on_right);
        let xmer = cmp::min(min_space, MAX_XMER);
        &sequence[center - xmer..center + xmer + 1]
    }
    /// Get `(short-range, long-range)` statistics for the given subsequence.
    fn score_subseq_smart<'a>(&self, subseq: &aa_canonical_str) -> impl ExactSizeIterator<Item = (f64, f64)> {
        debug_assert!(subseq.len() % 2 == 1);
        let midpoint = subseq.len() / 2;
        debug_assert!(midpoint <= MAX_XMER);
        let mut frequency_sum_a = 0.0;
        let mut frequency_sum_b = 0.0;
        let mut denom_total_a = 0.0;
        let mut denom_total_b = 0.0;
        let aa = subseq[midpoint];
        let cap = cmp::min(midpoint, MAX_XMER - 1);
        (0..cap).map(move |p| {
            let c_term_position = midpoint + 1 + p;
            let (freq_a, std_a) =
                self.query_freq_pair_db::<true, true>(aa, subseq[c_term_position], p);
            let (freq_b, std_b) =
                self.query_freq_pair_db::<true, false>(aa, subseq[c_term_position], p);
            frequency_sum_a += freq_a / std_a;
            frequency_sum_b += freq_b / std_b;
            denom_total_a += 1.0 / std_a;
            denom_total_b += 1.0 / std_b;
            let n_term_position = midpoint - 1 - p;
            let (freq_a, std_a) =
                self.query_freq_pair_db::<false, true>(subseq[n_term_position], aa, p);
            let (freq_b, std_b) =
                self.query_freq_pair_db::<false, false>(subseq[n_term_position], aa, p);
            frequency_sum_a += freq_a / std_a;
            frequency_sum_b += freq_b / std_b;
            denom_total_a += 1.0 / std_a;
            denom_total_b += 1.0 / std_b;
            let final_frequency_a = frequency_sum_a / denom_total_a;
            let final_frequency_b = frequency_sum_b / denom_total_b;
            (final_frequency_a, final_frequency_b)
        })
    }
    /// Method for getting a `(mean, std)` pair from the "pair-frequency-DB".
    /// 
    /// This method is meant to provide an opaque interface mimicking `GridScore`
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
        let sr_selector = (!SR) as usize;
        if X {
            self.pair_freq_db_x[distance][aa_a][aa_b][sr_selector]
        } else {
            self.pair_freq_db_y[distance][aa_b][aa_a][sr_selector]
        }
    }
    /// Method for getting a `(mean, std)` pair from the "avg-sdev-DB".
    /// 
    /// This method is meant to provide an interface mimicking `GridScore`
    /// class's `AvgSdevDB`, which is a dictionary that reports mean
    /// and standard deviation of a biophysical characteristic for x-mers
    /// containing the given `aa` that are `xmer` long.
    /// 
    /// You can find the `GridScore` class in Julie Forman-Kay's lab's
    /// [python package], written by Hao Cai (@haocai1992).
    /// 
    /// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
    fn query_avg_sdev_db<const SR: bool>(&self, aa: Aminoacid, xmer: usize) -> (f64, f64) {
        let sr_selector = (!SR) as usize;
        self.avg_sdev_db[aa][xmer][sr_selector]
    }
    /// Method for getting a `(count, short-range, long-range)`
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
        let table = self.z_grid_db[aa][xmer];
        let row = table.binary_search_by(|(gridpoint, _)| gridpoint.total_cmp(&sr_linekey)).ok()?;
        let (_, row) = table[row];
        let cell = row.binary_search_by(|(gridpoint, _)| gridpoint.total_cmp(&lr_linekey)).ok()?;
        let (_, tuple) = row[cell];
        Some(tuple)
    }
    /// Method for getting a `(count, short-range, long-range)`
    /// data tuple from the "z-grid-DB". This method is expected to be
    /// run on floats that were not near any gridpoints in the data
    /// struct via [`Self::query_zgrid_db_exact`].
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
        zscore_sr: f64,
        zscore_lr: f64,
    ) -> (f64, f64, f64) {
        let table = self.z_grid_db[aa][xmer];
        let mut best_entry: Option<(f64, usize, usize)> = None;
        for (sr_index, (sr_gridpoint, row)) in table.iter().enumerate() {
            debug_assert!(!row.is_empty());
            let sr_delta = sr_gridpoint.0 - zscore_sr;
            let sr_sqr_delta = sr_delta * sr_delta;
            if let Some((best_sqr_delta, sr_index, lr_index)) = best_entry {
                if best_sqr_delta <= sr_delta && sr_gridpoint.0 > zscore_sr {
                    return table[sr_index].1[lr_index].1
                }
            }
            let (Ok(lr_index) | Err(lr_index)) = row.binary_search_by(|(lr_gridpoint, _)| lr_gridpoint.0.total_cmp(&zscore_lr));
            let mut test_lr_index = |lr_index: usize| {
                if let Some((lr_gridpoint, _)) = row.get(lr_index) {
                    let lr_delta = lr_gridpoint.0 - zscore_lr;
                    let lr_sqr_delta = lr_delta * lr_delta;
                    let sqr_delta = sr_sqr_delta + lr_sqr_delta;
                    if let Some(entry) = best_entry {
                        if sqr_delta < entry.0 {
                            best_entry = Some((sqr_delta, sr_index, lr_index))
                        }
                    }
                }
            };
            test_lr_index(lr_index);
            let Some(lr_index) = lr_index.checked_sub(1) else {
                continue;
            };
            test_lr_index(lr_index);
        } 
        let (_, sr_index, lr_index) = best_entry.unwrap();
        table[sr_index].1[lr_index].1
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
    pub fn total_cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}