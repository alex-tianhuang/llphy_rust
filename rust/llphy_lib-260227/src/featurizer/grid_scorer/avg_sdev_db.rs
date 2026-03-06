//! Module defining [`AvgSdevDB`],
//! a wrapper around a an `(aa, xmer)`-indexable
//! collection of [`AvgSdevDBEntry`]s.
use std::{mem::MaybeUninit, ops::{Deref, DerefMut}};
use anyhow::Error;
use borsh::{BorshDeserialize, BorshSerialize};

use crate::{datatypes::{AAMap, MAX_XMER}, featurizer::grid_scorer::{XmerIndexableArray, AvgSdevDBEntry}};


/// A struct containing [`AvgSdevDBEntry`] for each
/// `(aa, xmer)` key.
#[derive(BorshSerialize, PartialEq)]
#[repr(transparent)]
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
    /// Given space for a new [`AvgSdevDB`] table,
    /// initialize all its entries to be `f64::NAN`s.
    pub fn init_with_nans(this: &mut MaybeUninit<Self>) -> &mut Self {
        // SAFETY: AvgSdevDB has #[repr(transparent)] AAMap(XmerIndexableArray)
        //         AAMap has #[repr(transparent)] over [_; 20]
        //         XmerIndexableArray has #[repr(transparent)] over [_; MAX_XMER]
        //         no other initialization guarantees
        let table = unsafe {&mut *this.as_mut_ptr().cast::<[[MaybeUninit<AvgSdevDBEntry>;MAX_XMER]; 20]>()};
        for arr in table {
            for slot in arr {
                slot.write(AvgSdevDBEntry::new_nan_filled());
            }
        }
        // SAFETY: just initialized all slots
        unsafe { this.assume_init_mut() }
    }
    /// True if there are `f64::NAN`s in any entry.
    pub fn is_nan_free(&self) -> bool {
        !self.0.values().flatten().any(|e| !e.is_nan_free())
    }
    /// Moral equivalent of implementing deserialization on [`AvgSdevDB`],
    /// but does it by reference so that don't have to put it on the stack.
    pub fn deserialize_into(this: &mut MaybeUninit<Self>, buf: &mut &[u8]) -> Result<(), Error> {
        let this = Self::init_with_nans(this);
        for arr in this.values_mut() {
            for slot in arr {
                let data = AvgSdevDBEntry::deserialize(buf)?;
                *slot = data;
            }
        }
        debug_assert!(this.is_nan_free());
        Ok(())
    }
}