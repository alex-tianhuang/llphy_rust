//! Module defining common implementation details
//! for [`ZGridDB`], a struct used to query
//! for values in a grid of z-scores.
//!
//! Separated from simd/non-simd specific implementations.

use anyhow::Error;
use borsh::BorshSerialize;
use bumpalo::Bump;
use std::{mem::MaybeUninit, ops::{Deref, DerefMut}};

use crate::{
    datatypes::AAMap,
    featurizer::grid_scorer::{XmerIndexableArray, ZGridSubtable},
};

/// A struct of tables indexable by `(aa, xmer)` keys,
/// where each subtable is 2D `zscore`-indexable.
/// See [`ZGridSubtable`].
///
/// The subtable is expected to be indexed by `zscore_a`
/// and then `zscore_b`, NOT `zscore_b` then `zscore_a`.
#[derive(PartialEq)]
#[repr(transparent)]
pub struct ZGridDB<'a>(AAMap<XmerIndexableArray<ZGridSubtable<'a>>>);
/// The maximum number of Z-score cells in `ZGridSubtable`s
/// derived from Cai's old `LLPhyScore` executable.
/// 
/// Unless the `gridscorer.bin` files get corrupted or
/// somebody brave attempts to add new biophysical features,
/// this number should never be exceeded when deserializing
/// a `ZGridDB`.
pub const KNOWN_MAX_DATA_LEN: usize = 575;
impl<'a> Deref for ZGridDB<'a> {
    type Target = AAMap<XmerIndexableArray<ZGridSubtable<'a>>>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl<'a> DerefMut for ZGridDB<'a> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl<'a> ZGridDB<'a> {
    /// Convert this uninitialized [`ZGridDB`]
    /// to an uninitialized version of its inner field.
    pub fn as_uninit_inner(this: &mut MaybeUninit<Self>) -> &mut AAMap<XmerIndexableArray<MaybeUninit<ZGridSubtable<'a>>>> {
        // SAFETY: ZGridDB has #[repr(transparent)] over its inner field.
        unsafe { 
            &mut *this.as_mut_ptr().cast::<AAMap<XmerIndexableArray<MaybeUninit<ZGridSubtable<'_>>>>>()
        }
    }
    /// Moral equivalent of implementing deserialization on [`ZGridDB`],
    /// but does it by reference so that don't have to put it on the stack.
    pub fn deserialize_into(this: &mut MaybeUninit<Self>, buf: &mut &[u8], arena: &'a Bump) -> Result<(), Error> {
        let this = Self::as_uninit_inner(this);
        for arr in this.values_mut() {
            for slot in arr {
                slot.write(ZGridSubtable::deserialize(buf, arena)?);
            }
        }
        Ok(())
    }
}

impl BorshSerialize for ZGridDB<'_> {
    fn serialize<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        for arr in self.values() {
            for subtable in arr {
                subtable.serialize(writer)?;
            }
        }
        Ok(())
    }
}
/// Find the gridpoint that minimizes the
/// sum of squared differences to the given zscores.
///
/// Helper function for [`ZGridSubtable::lookup`].
///
/// Arguments
/// ---------
/// - `data` is a `row_len x column_len` grid of entries of type `T`,
///   which may or may not be occupied with data
/// - `to_occupied` attempts to convert `T` into an occupied entry
///   of type `U`
/// - `coord_a` and `coord_b` are the two coordinates for which
///   we are trying to find an entry with the minimum distance to.
pub fn lookup_thorough<T, U>(
    data: &[T],
    row_len: usize,
    to_occupied: fn(&T) -> Option<&U>,
    coord_a: f64,
    coord_b: f64,
) -> &U {
    let column_len = data.len() / row_len;
    let mut best_entry: Option<(f64, &U)> = None;
    let mut down_idx_a = (coord_a as usize + 1).min(column_len);
    let mut up_idx_a = down_idx_a;
    let mut search_up = coord_a.round() as usize >= up_idx_a && up_idx_a < column_len;
    for _ in 0..column_len {
        let idx_a = if search_up { up_idx_a } else { down_idx_a - 1 };
        let sqr_delta_a = {
            let delta_a = idx_a as f64 - coord_a;
            delta_a * delta_a
        };
        if let Some((best_sqr_delta, entry)) = best_entry {
            if best_sqr_delta <= sqr_delta_a {
                return entry;
            }
        }
        let row = &data[idx_a * row_len..(idx_a + 1) * row_len];
        debug_assert!(!row.is_empty());
        let idx_is_occupied = |idx| to_occupied(&row[idx]).is_some();
        // call this only if `idx_is_occupied(idx_b)` is true
        let mut check_idx_b = |idx_b: usize| {
            let delta_b = idx_b as f64 - coord_b;
            let sqr_delta_b = delta_b * delta_b;
            let sqr_delta = sqr_delta_a + sqr_delta_b;
            let new_best: bool;
            if let Some(entry) = best_entry {
                new_best = sqr_delta < entry.0;
            } else {
                new_best = true;
            }
            if new_best {
                best_entry = Some((sqr_delta, to_occupied(&row[idx_b]).unwrap()))
            }
        };
        let breakpoint_b = (coord_b as usize + 1).min(row_len);
        for idx_b in breakpoint_b..row_len {
            if idx_is_occupied(idx_b) {
                check_idx_b(idx_b);
                break;
            }
        }
        for idx_b in (0..breakpoint_b).rev() {
            if idx_is_occupied(idx_b) {
                check_idx_b(idx_b);
                break;
            }
        }
        if search_up {
            up_idx_a += 1;
            search_up = down_idx_a == 0;
        } else {
            down_idx_a -= 1;
            search_up = up_idx_a < column_len;
        }
    }
    best_entry.unwrap().1
}
