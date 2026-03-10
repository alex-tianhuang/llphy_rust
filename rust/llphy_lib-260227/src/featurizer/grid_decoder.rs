//! Module defining [`GridDecoder`].
use borsh::{BorshDeserialize, BorshSerialize};
use crate::datatypes::{AAIndex, AAMap};
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
    /// Compute the residue-level scores using
    /// [`Self::thresholds`]' thresholding rules,
    /// returning an iterator to be consumed by something else.
    pub fn compute_residue_level(&self, sequence: &[AAIndex], grid_score: &AAMap<&[f64]>) -> impl ExactSizeIterator<Item = (usize, i64)> {
        let subsequence: &[AAIndex];
        match sequence.len().checked_sub(2) {
            Some(n_sites) => {
                subsequence = &sequence[1..=n_sites];
            },
            None => {
                subsequence = &[]
            }
        }
        let mut cursors = AAMap(grid_score.0.map(|v|v.len()));
        subsequence.into_iter().enumerate().map(move |(i, &aa)| {
            let cursor = cursors[aa];
            cursors[aa] += 1;
            let score = self.thresholds[aa].score_site(grid_score[aa][cursor]) * self.sign as i64;
            (i, score)
        })
    }
}