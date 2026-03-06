//! Module defining [`PairFreqDBEntry`] and [`PairFreqEntrySum`]
//! using `#[portable_simd]`.
use std::{ops::AddAssign, simd::{f64x2, f64x4}};
use crate::{derive_borsh_de_from, derive_borsh_se_into};

/// Weights for each `(aa_x, gap_length, xy_orientation, aa_y)` key.
/// 
/// The field names are not visible in the SIMD representation,
/// but it is essentially an array consisting of named floats
/// `[weight_a, weight_b, total_a, total_b]`.
/// 
/// Dev note
/// --------
/// The reason this contains data for two features instead
/// of separating them out nicely is because a pair of zscores
/// are required to index into [`crate::featurizer::grid_scorer::ZGridDB`]
/// and so it is a win for cache locality to get them loaded in one struct.
pub struct PairFreqDBEntry(f64x4);
/// An accumulator struct for taking the sum of many [`PairFreqDBEntry`]s.
pub struct PairFreqEntrySum(f64x4);

impl PairFreqDBEntry {
    /// Get a new [`PairFreqDBEntry`] that is filled with `f64::NAN`.
    pub const fn new_nan_filled() -> Self {
        Self(f64x4::splat(f64::NAN))
    }
    /// False if there are `f64::NAN`s in any field.
    pub fn is_nan_free(&self) -> bool {
        self.0.as_array()
            .into_iter()
            .all(|x| !x.is_nan())
    }
    /// Set weights and total for feature `A`.
    pub fn set_a(&mut self, weight: f64, total: f64) {
        let this = self.0.as_mut_array();
        this[0] = weight;
        this[2] = total;
    }
    /// Set weights and total for feature `B`.
    pub fn set_b(&mut self, weight: f64, total: f64) {
        let this = self.0.as_mut_array();
        this[1] = weight;
        this[3] = total;
    }
}
impl PairFreqEntrySum {
    /// Get a new [`PairFreqEntrySum`] that is filled with `0.0_f64`.
    pub fn new_zeroed() -> Self {
        Self(f64x4::from_array([0.0; 4]))
    }
    /// Helper method for [`crate::featurizer::GridScorer::score_sequence`].
    /// 
    /// Equivalent to:
    /// ```
    /// [self.weight_a / self.total_a, self.weight_b / self.total_b]
    /// ```
    pub fn as_frequencies(&self) -> f64x2 {
        self.0.extract::<0, 2>() / self.0.extract::<2, 2>()
    }
}
impl AddAssign<&PairFreqDBEntry> for PairFreqEntrySum {
    fn add_assign(&mut self, rhs: &PairFreqDBEntry) {
        self.0 += &rhs.0
    }
}
derive_borsh_de_from!(PairFreqDBEntry as [f64; 4], |a| PairFreqDBEntry(f64x4::from_array(a)));
derive_borsh_se_into!(PairFreqDBEntry as [f64; 4], |this| this.0.to_array());