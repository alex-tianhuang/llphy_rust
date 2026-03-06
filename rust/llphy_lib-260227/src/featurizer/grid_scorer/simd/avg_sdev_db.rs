//! Module defining [`AvgSdevDBEntry`],
//! using `#[portable_simd]`.
use crate::{derive_borsh_de_from, derive_borsh_se_into};
use std::simd::f64x2;

/// A struct containing average and inverse standard deviations
/// of two features for a given `(aa, xmer)` key.
///
/// Dev note
/// --------
/// The reason this contains data for two features instead
/// of separating them out nicely is because a pair of zscores
/// are required to index into [`crate::featurizer::grid_scorer::ZGridDB`]
/// and so it is a win for cache locality to get them loaded in one struct.
pub struct AvgSdevDBEntry {
    avgs: f64x2,
    invstds: f64x2,
}
impl AvgSdevDBEntry {
    /// Get a new [`AvgSdevDBEntry`] that is filled with `f64::NAN`.
    pub const fn new_nan_filled() -> Self {
        Self {
            avgs: f64x2::splat(f64::NAN),
            invstds: f64x2::splat(f64::NAN),
        }
    }
    /// False if there are `f64::NAN`s in any field.
    pub fn is_nan_free(&self) -> bool {
        [self.avgs[0], self.avgs[1], self.invstds[0], self.invstds[1]]
            .into_iter()
            .all(|x| !x.is_nan())
    }
    /// Set mean and inverse std for feature `A`.
    pub fn set_a(&mut self, avg: f64, inv_std: f64) {
        self.avgs[0] = avg;
        self.invstds[0] = inv_std;
    }
    /// Set mean and inverse std for feature `B`.
    pub fn set_b(&mut self, avg: f64, inv_std: f64) {
        self.avgs[1] = avg;
        self.invstds[1] = inv_std;
    }
    /// Helper method for [`crate::featurizer::GridScorer::score_sequence`].
    ///
    /// Equivalent to the non-simd:
    /// ```
    /// let freqs: [f64; 2];
    /// let zscore_a = (freqs[0] - self.avg_a) * self.invstd_a;
    /// let zscore_b = (freqs[1] - self.avg_b) * self.invstd_b;
    /// [zscore_a, zscore_b]
    /// ```
    pub fn freqs_to_zscores(&self, freqs: f64x2) -> f64x2 {
        (freqs - self.avgs) * self.invstds
    }
}
derive_borsh_de_from!(AvgSdevDBEntry as [[f64; 2]; 2], |[avgs, invstds]| {
    AvgSdevDBEntry {
        avgs: f64x2::from_array(avgs),
        invstds: f64x2::from_array(invstds),
    }
});
derive_borsh_se_into!(AvgSdevDBEntry as [[f64; 2]; 2], |this| [this.avgs.to_array(), this.invstds.to_array()]);

