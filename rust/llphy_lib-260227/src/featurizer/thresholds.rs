//! Module defining a thresholds container for
//! converting grid-scores to integer values.
//!
//! See [`Thresholds`] and [`ThresholdPair`].
use borsh::{BorshDeserialize, BorshSerialize};

use crate::datatypes::AAMap;
use std::ops::{Deref, DerefMut};

/// A collection of [`ThresholdPair`]s for each aminoacid.
#[derive(BorshSerialize, BorshDeserialize, PartialEq)]
pub struct Thresholds(AAMap<ThresholdPair>);
/// An upper and lower grid-score threshold
/// that converts a grid-score to an integer value.
#[derive(BorshSerialize, BorshDeserialize, PartialEq)]
pub struct ThresholdPair {
    pub upper: f64,
    pub lower: f64,
}
impl Deref for Thresholds {
    type Target = AAMap<ThresholdPair>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for Thresholds {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl Thresholds {
    /// Get a new [`Thresholds`] struct that is filled with `f64::NAN`.
    pub fn new_nan_filled() -> Self {
        Self(AAMap(
            [const {
                ThresholdPair {
                    lower: f64::NAN,
                    upper: f64::NAN,
                }
            }; 20],
        ))
    }
    /// True if there are `f64::NAN`s in any slot of the thresholds.
    pub fn is_nan_free(&self) -> bool {
        !self
            .0
            .0
            .iter()
            .any(|t| t.upper.is_nan() || t.lower.is_nan())
    }
}
impl ThresholdPair {
    /// Get an integer value from a single grid score.
    fn score_site(&self, grid_score: f64) -> i64 {
        let Self { lower, upper } = *self;
        (grid_score > upper) as i64 - (grid_score < lower) as i64
    }
    /// Sum the values of several grid scores.
    pub fn score_sites(&self, grid_scores: &[f64]) -> i64 {
        grid_scores.iter().map(|x| self.score_site(*x)).sum::<i64>()
    }
}
