//! Module defining [`AvgSdevDB`].
use crate::datatypes::{AAMap, MAX_XMER, llphy::features::grid_scorer::xmer::XmerIndexableArray};
use std::ops::{Deref, DerefMut};

/// A struct containing [`AvgSdevDBEntry`] for each
/// `(aa, xmer)` key.
pub struct AvgSdevDB(AAMap<XmerIndexableArray<AvgSdevDBEntry>>);
/// A struct containing average and standard deviations
/// of two features for a given `(aa, xmer)` key.
///
/// Dev note
/// --------
/// The reason this contains data for two features instead
/// of separating them out nicely is because a pair of zscores
/// are required to index into [`super::ZGridDB`] and so it is a
/// win for cache locality to get them loaded in one struct.
pub struct AvgSdevDBEntry {
    pub avg_a: f64,
    pub std_a: f64,
    pub avg_b: f64,
    pub std_b: f64,
}
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
    #[cfg(debug_assertions)]
    /// True if there are `f64::NAN`s in any entry.
    pub fn is_nan_free(&self) -> bool {
        !self.0.values().flatten().any(|e| e.is_nan_free())
    }
}
impl AvgSdevDBEntry {
    /// Get a new [`AvgSdevDBEntry`] that is filled with `f64::NAN`.
    const fn new_nan_filled() -> Self {
        Self {
            avg_a: f64::NAN,
            std_a: f64::NAN,
            avg_b: f64::NAN,
            std_b: f64::NAN,
        }
    }
    #[cfg(debug_assertions)]
    /// False if there are `f64::NAN`s in any field.
    pub fn is_nan_free(&self) -> bool {
        [self.avg_a, self.std_a, self.avg_b, self.std_b]
            .into_iter()
            .all(|x| !x.is_nan())
    }
}
