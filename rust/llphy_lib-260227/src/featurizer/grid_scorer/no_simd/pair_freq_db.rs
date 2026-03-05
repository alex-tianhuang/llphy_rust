//! Module defining [`PairFreqDBEntry`] and [`PairFreqEntrySum`]
//! without the use of the `#[portable_simd]` feature.
use std::ops::AddAssign;

/// Weights for each `(aa_x, gap_length, xy_orientation, aa_y)` key.
///
/// Dev note
/// --------
/// The reason this contains data for two features instead
/// of separating them out nicely is because a pair of zscores
/// are required to index into [`super::ZGridDB`] and so it is a
/// win for cache locality to get them loaded in one struct.
pub struct PairFreqDBEntry {
    weight_a: f64,
    weight_b: f64,
    total_a: f64,
    total_b: f64,
}
/// An accumulator struct for taking the sum of many [`PairFreqDBEntry`]s.
pub struct PairFreqEntrySum {
    weight_a: f64,
    weight_b: f64,
    total_a: f64,
    total_b: f64,
}
impl PairFreqDBEntry {
    /// Get a new [`PairFreqDBEntry`] that is filled with `f64::NAN`.
    pub const fn new_nan_filled() -> Self {
        Self {
            weight_a: f64::NAN,
            weight_b: f64::NAN,
            total_a: f64::NAN,
            total_b: f64::NAN,
        }
    }
    /// False if there are `f64::NAN`s in any field.
    pub fn is_nan_free(&self) -> bool {
        [self.weight_a, self.weight_b, self.total_a, self.total_b]
            .into_iter()
            .all(|x| !x.is_nan())
    }
    /// Set weights and total for feature `A`.
    pub fn set_a(&mut self, weight: f64, total: f64) {
        self.weight_a = weight;
        self.total_a = total;
    }
    /// Set weights and total for feature `B`.
    pub fn set_b(&mut self, weight: f64, total: f64) {
        self.weight_b = weight;
        self.total_b = total;
    }
}
impl PairFreqEntrySum {
    /// Get a new [`PairFreqEntrySum`] that is filled with `0.0_f64`.
    pub fn new_zeroed() -> Self {
        Self {
            weight_a: 0.0,
            weight_b: 0.0,
            total_a: 0.0,
            total_b: 0.0,
        }
    }
    /// Helper method for [`crate::featurizer::GridScorer::score_sequence`].
    ///
    /// Equivalent to:
    /// ```
    /// [self.weight_a / self.total_a, self.weight_b / self.total_b]
    /// ```
    pub fn as_frequencies(&self) -> [f64; 2] {
        [self.weight_a / self.total_a, self.weight_b / self.total_b]
    }
}
impl AddAssign<&PairFreqDBEntry> for PairFreqEntrySum {
    fn add_assign(&mut self, rhs: &PairFreqDBEntry) {
        self.weight_a += rhs.weight_a;
        self.weight_b += rhs.weight_b;
        self.total_a += rhs.total_a;
        self.total_b += rhs.total_b;
    }
}
