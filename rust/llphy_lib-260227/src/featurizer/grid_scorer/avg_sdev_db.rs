//! Module defining [`AvgSdevDB`].
use crate::datatypes::{AAMap, MAX_XMER};
use crate::featurizer::grid_scorer::xmer::XmerIndexableArray;
use std::ops::{Deref, DerefMut};
use std::simd::{f64x4, f64x2};

/// A struct containing [`AvgSdevDBEntry`] for each
/// `(aa, xmer)` key.
pub struct AvgSdevDB(AAMap<XmerIndexableArray<AvgSdevDBEntry>>);
/// A struct containing average and inverse standard deviations
/// of two features for a given `(aa, xmer)` key.
///
/// The field names are not visible in the SIMD representation,
/// but it is essentially an array consisting of named floats
/// [`avg_a`], [`avg_b`], [`invstd_a`], and [`invstd_b`].
/// 
/// [`avg_a`]: Self::avg_a
/// [`avg_b`]: Self::avg_b
/// [`invstd_a`]: Self::invstd_a
/// [`invstd_b`]: Self::invstd_b
/// 
/// Dev note
/// --------
/// The reason this contains data for two features instead
/// of separating them out nicely is because a pair of zscores
/// are required to index into [`super::ZGridDB`] and so it is a
/// win for cache locality to get them loaded in one struct.
pub struct AvgSdevDBEntry(f64x4);
impl Deref for AvgSdevDB {
    type Target = AAMap<XmerIndexableArray<AvgSdevDBEntry>>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for AvgSdevDB {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl AvgSdevDB {
    /// Get a new [`AvgSdevDB`] struct that is filled with `f64::NAN`.
    pub fn new_nan_filled() -> Self {
        AvgSdevDB(AAMap(
            [const { XmerIndexableArray::new([const { AvgSdevDBEntry::new_nan_filled() }; MAX_XMER]) };
                20],
        ))
    }
    /// True if there are `f64::NAN`s in any entry.
    pub fn is_nan_free(&self) -> bool {
        !self.0.values().flatten().any(|e| !e.is_nan_free())
    }
}
/// Shorthand for associating each index of the
/// array of [`AvgSdevDBEntry`] with a named field.
macro_rules! impl_getters {
    ($([$index:literal, $field:ident]),*) => {
        impl AvgSdevDBEntry {
            $(pub fn $field(&self) -> f64 {
                self.0.as_array()[$index]
            })*
        }
    }; 
}
impl_getters!([0, avg_a], [1, avg_b], [2, invstd_a], [3, invstd_b]);
impl AvgSdevDBEntry {
    /// Get a new [`AvgSdevDBEntry`] that is filled with `f64::NAN`.
    const fn new_nan_filled() -> Self {
        Self(f64x4::from_array([f64::NAN; 4]))
    }
    /// False if there are `f64::NAN`s in any field.
    pub fn is_nan_free(&self) -> bool {
        self.0.as_array()
            .into_iter()
            .all(|x| !x.is_nan())
    }
    /// Set mean and inverse std for feature `A`.
    pub fn set_a(&mut self, avg: f64, inv_std: f64) {
        let this = self.0.as_mut_array();
        this[0] = avg;
        this[2] = inv_std;
    }
    /// Set mean and inverse std for feature `B`.
    pub fn set_b(&mut self, avg: f64, inv_std: f64) {
        let this = self.0.as_mut_array();
        this[1] = avg;
        this[3] = inv_std;
    }
    /// Equivalent to:
    /// ```
    /// let freqs: f64x2;
    /// let zscore_a = (freqs.as_array()[0] - self.avg_a()) / self.std_a();
    /// let zscore_b = (freqs.as_array()[1] - self.avg_b()) / self.std_b();
    /// [zscore_a, zscore_b]
    /// ```
    pub fn freqs_to_zscores(&self, freqs: f64x2) -> f64x2 {
        (freqs - self.0.extract::<0, 2>()) * self.0.extract::<2, 2>()
    }
}
