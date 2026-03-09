//! Module defining [`AvgSdevDBEntry`]
//! without the use of the `#[portable_simd]` feature.
use crate::{derive_borsh_de_from, derive_borsh_se_into};
/// A struct containing average and inverse standard deviations
/// of two features for a given `(aa, xmer)` key.
///
/// Dev note
/// --------
/// The reason this contains data for two features instead
/// of separating them out nicely is because a pair of zscores
/// are required to index into [`crate::featurizer::grid_scorer::ZGridDB`]
/// and so it is a win for cache locality to get them loaded in one struct.
#[derive(PartialEq)]
pub struct AvgSdevDBEntry {
    avg_a: f64,
    avg_b: f64,
    invstd_a: f64,
    invstd_b: f64,
}
impl AvgSdevDBEntry {
    /// Get a new [`AvgSdevDBEntry`] that is filled with `f64::NAN`.
    pub const fn new_nan_filled() -> Self {
        Self {
            avg_a: f64::NAN,
            avg_b: f64::NAN,
            invstd_a: f64::NAN,
            invstd_b: f64::NAN,
        }
    }
    /// False if there are `f64::NAN`s in any field.
    pub fn is_nan_free(&self) -> bool {
        [self.avg_a, self.avg_b, self.invstd_a, self.invstd_b]
            .into_iter()
            .all(|x| !x.is_nan())
    }
    /// Set mean and inverse std for feature `A`.
    pub fn set_a(&mut self, avg: f64, inv_std: f64) {
        self.avg_a = avg;
        self.invstd_a = inv_std;
    }
    /// Set mean and inverse std for feature `B`.
    pub fn set_b(&mut self, avg: f64, inv_std: f64) {
        self.avg_b = avg;
        self.invstd_b = inv_std;
    }
    /// Helper method for [`crate::featurizer::GridScorer::score_sequence`].
    ///
    /// Equivalent to:
    /// ```
    /// let freqs: [f64; 2];
    /// let zscore_a = (freqs[0] - self.avg_a) * self.invstd_a;
    /// let zscore_b = (freqs[1] - self.avg_b) * self.invstd_b;
    /// [zscore_a, zscore_b]
    /// ```
    pub fn freqs_to_zscores(&self, freqs: [f64; 2]) -> [f64; 2] {
        let zscore_a = (freqs[0] - self.avg_a) * self.invstd_a;
        let zscore_b = (freqs[1] - self.avg_b) * self.invstd_b;
        [zscore_a, zscore_b]
    }
}
derive_borsh_de_from!(AvgSdevDBEntry as [[f64; 2]; 2], |[
    [avg_a, avg_b],
    [invstd_a, invstd_b],
]| AvgSdevDBEntry {
    avg_a,
    avg_b,
    invstd_a,
    invstd_b
});
derive_borsh_se_into!(AvgSdevDBEntry as [[f64; 2]; 2], |this| [
    [this.avg_a, this.avg_b],
    [this.invstd_a, this.invstd_b]
]);
