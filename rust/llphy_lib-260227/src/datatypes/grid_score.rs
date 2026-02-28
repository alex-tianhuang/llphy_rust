use crate::{
    datatypes::{AAMap, Aminoacid, aa_canonical_str},
    leak_vec,
};
use bumpalo::{Bump, collections::Vec};
use std::{array, cmp};
const MIN_XMER: usize = 1;
pub struct GridScoreOld {}
pub struct ResidueDataOld {
    aa: Aminoacid,
    // n_outside_grid: usize,
    sr: f64,
    lr: f64,
}
pub type ResScoresOld<'a> = AAMap<&'a [ResScoresEntryOld]>;
#[derive(Clone)]
pub struct ResScoresEntryOld {
    pub sr: f64,
    pub lr: f64,
}
impl GridScoreOld {
    pub fn score_sequence_to_res_scores<'a>(
        &self,
        sequence: &aa_canonical_str,
        arena: &'a Bump,
    ) -> ResScoresOld<'a> {
        let mut residue_counts = AAMap([0_usize; 20]);
        for aa in sequence {
            residue_counts[aa] += 1;
        }
        let mut grid = AAMap(array::from_fn(|aaindex| {
            let cap_bound = residue_counts.0[aaindex];
            Vec::with_capacity_in(cap_bound, arena)
        }));
        for ResidueDataOld { aa, sr, lr } in self.score_sequence(sequence) {
            grid[aa].push(ResScoresEntryOld { sr, lr })
        }
        AAMap(grid.0.map(|buf| leak_vec(buf) as &[_]))
    }
    pub fn score_sequence<'a>(&self, sequence: &aa_canonical_str) -> std::vec::Vec<ResidueDataOld> {
        let Some(n_sites) = sequence.len().checked_sub(2 * MIN_XMER) else {
            return std::vec::Vec::new();
        };
        let mut score_data = std::vec::Vec::with_capacity(n_sites);
        for i in 1..=n_sites {
            let mid_res = sequence[i];
            let mut grid_counts = [0.0, 0.0, 0.0];
            // let mut n_outside_grid = 0;
            let m_seq = self.get_middle_seq_centered_at(sequence, i);
            // let mut avg_zscore_sr = 0.0;
            // let mut avg_zscore_lr = 0.0;
            // let mut zscore_n = 0.0;
            debug_assert!(m_seq.len() >= 2 * MIN_XMER + 1);
            let all_xmer_scores = self.score_mseq_smart(m_seq);
            for (xmer, (freq_sr, freq_lr)) in (MIN_XMER..).zip(all_xmer_scores) {
                let (mean_sr, std_sr) = self.query_avg_sdev_db::<true>(mid_res, xmer);
                let (mean_lr, std_lr) = self.query_avg_sdev_db::<false>(mid_res, xmer);
                let zscore_sr = (freq_sr - mean_sr) / std_sr;
                let zscore_lr = (freq_lr - mean_lr) / std_lr;
                let sr_key = make_linekey(zscore_sr);
                let lr_key = make_linekey(zscore_lr);
                // avg_zscore_sr += zscore_sr;
                // avg_zscore_lr += zscore_lr;
                // zscore_n += 1.0;
                let [total_slot, sr_slot, lr_slot] = &mut grid_counts;
                let (total, sr, lr) = self
                    .query_zgrid_db_exact(mid_res, xmer, sr_key, lr_key)
                    .unwrap_or_else(|| {
                        // n_outside_grid += 1;
                        self.query_zgrid_db_nearest(mid_res, xmer, sr_key, lr_key)
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
                // avg_zscore_sr /= zscore_n;
                // avg_zscore_lr /= zscore_n;
            }
            score_data.push(ResidueDataOld {
                aa: mid_res,
                // n_outside_grid,
                sr: freq_by_grid_sr,
                lr: freq_by_grid_lr,
            });
        }
        score_data
    }
    fn get_middle_seq_centered_at<'a>(
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
    fn max_xmer(&self) -> usize {
        todo!()
    }
    fn score_mseq_smart<'a>(&self, m_seq: &aa_canonical_str) -> std::vec::Vec<(f64, f64)> {
        debug_assert!(m_seq.len() % 2 == 1);
        let midpoint = m_seq.len() / 2;
        debug_assert!(midpoint <= self.max_xmer());
        let mut frequency_sum_a = 0.0;
        let mut frequency_sum_b = 0.0;
        let mut denom_total_a = 0.0;
        let mut denom_total_b = 0.0;
        let mid_res_aa = m_seq[midpoint];
        let cap = cmp::min(midpoint, self.max_xmer() - 1);
        let mut scores = std::vec::Vec::with_capacity(cap);
        for p in 0..midpoint {
            let c_term_position = midpoint + 1 + p;
            let (freq_a, std_a) =
                self.query_freq_pair_db::<true, true>(mid_res_aa, m_seq[c_term_position], p);
            let (freq_b, std_b) =
                self.query_freq_pair_db::<true, false>(mid_res_aa, m_seq[c_term_position], p);
            frequency_sum_a += freq_a / std_a;
            frequency_sum_b += freq_b / std_b;
            denom_total_a += 1.0 / std_a;
            denom_total_b += 1.0 / std_b;
            let n_term_position = midpoint - 1 - p;
            let (freq_a, std_a) =
                self.query_freq_pair_db::<false, true>(m_seq[n_term_position], mid_res_aa, p);
            let (freq_b, std_b) =
                self.query_freq_pair_db::<false, false>(m_seq[n_term_position], mid_res_aa, p);
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
    fn query_freq_pair_db<const X: bool, const TAG_A: bool>(
        &self,
        aa_a: Aminoacid,
        aa_b: Aminoacid,
        distance: usize,
    ) -> (f64, f64) {
        todo!()
    }
    fn query_avg_sdev_db<const SR: bool>(&self, mid_res: Aminoacid, xmer: usize) -> (f64, f64) {
        todo!()
    }
    fn query_zgrid_db_exact(
        &self,
        mid_res: Aminoacid,
        xmer: usize,
        sr_linekey: f64,
        lr_linekey: f64,
    ) -> Option<(f64, f64, f64)> {
        todo!()
    }
    fn query_zgrid_db_nearest(
        &self,
        mid_res: Aminoacid,
        xmer: usize,
        sr_linekey: f64,
        lr_linekey: f64,
    ) -> (f64, f64, f64) {
        todo!()
    }
}
fn make_linekey(zscore: f64) -> f64 {
    ((zscore / 2.0).round() * 2.0).clamp(-8.0, 12.0)
}
