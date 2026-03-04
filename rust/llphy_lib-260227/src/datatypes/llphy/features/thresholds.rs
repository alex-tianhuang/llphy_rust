use std::ops::{Deref, DerefMut};

use crate::datatypes::AAMap;

pub struct Thresholds(AAMap<ThresholdPair>);
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
    pub fn nan_thresholds() -> AAMap<ThresholdPair> {
        AAMap(
            [const {
                ThresholdPair {
                    lower: f64::NAN,
                    upper: f64::NAN,
                }
            }; 20],
        )
    }
    pub fn is_filled_with_nan(&self) -> bool {
        self.0.0.iter().all(|t| t.upper.is_nan() && t.lower.is_nan())
    }
    pub fn is_nan_free(&self) -> bool {
        !self.0.0.iter().any(|t| t.upper.is_nan() || t.lower.is_nan())
    }
}
pub struct ThresholdPair {
    pub upper: f64,
    pub lower: f64,
}
impl ThresholdPair {
    fn score_site(&self, grid_score: f64) -> i64 {
        let Self { lower, upper } = *self;
        (grid_score > upper) as i64 - (grid_score < lower) as i64
    }
    pub fn score_sites(&self, grid_scores: &[f64]) -> i64 {
        grid_scores.iter().map(|x| self.score_site(*x)).sum::<i64>()
    }
}
