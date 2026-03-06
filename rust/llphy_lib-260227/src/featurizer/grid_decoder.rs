//! Module defining [`GridDecoder`].
use borsh::{BorshDeserialize, BorshSerialize};
use crate::datatypes::AAMap;
use crate::featurizer::thresholds::Thresholds;

/// Turns a subarray of a [`GridScore`](crate::featurizer::grid_scorer::GridScore)
/// into a sequence-level numeric feature.
#[derive(BorshSerialize, BorshDeserialize, PartialEq)]
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