//! Module defining [`GridDecoder`].
use crate::datatypes::AAMap;
use crate::features::thresholds::Thresholds;

/// Turns a subarray of a [`GridScore`](crate::features::GridScore)
/// into a sequence-level numeric feature.
pub struct GridDecoder {
    pub sign: i8,
    pub thresholds: Thresholds,
}
impl GridDecoder {
    /// Compute the sequence level feature
    /// from the given grid scores.
    pub fn compute(&self, grid_score: &AAMap<&[f64]>) -> i64 {
        let mut score = 0;
        for (t, &subarr) in self.thresholds.values().zip(grid_score.values()) {
            score += t.score_sites(subarr)
        }
        score * self.sign as i64
    }
}
