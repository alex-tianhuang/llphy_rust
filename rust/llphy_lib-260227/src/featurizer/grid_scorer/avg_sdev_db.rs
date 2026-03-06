//! Module defining [`AvgSdevDB`],
//! a wrapper around a an `(aa, xmer)`-indexable
//! collection of [`AvgSdevDBEntry`]s.
use std::ops::{Deref, DerefMut};
use borsh::{BorshDeserialize, BorshSerialize};

use crate::{datatypes::{AAMap, MAX_XMER}, featurizer::grid_scorer::{XmerIndexableArray, AvgSdevDBEntry}};


/// A struct containing [`AvgSdevDBEntry`] for each
/// `(aa, xmer)` key.
#[derive(BorshDeserialize, BorshSerialize, PartialEq)]
pub struct AvgSdevDB(AAMap<XmerIndexableArray<AvgSdevDBEntry>>);
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