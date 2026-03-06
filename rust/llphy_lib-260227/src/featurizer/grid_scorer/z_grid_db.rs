//! Module defining common implementation details
//! for [`ZGridDB`], a struct used to query
//! for values in a grid of z-scores.
//!
//! Separated from simd/non-simd specific implementations.

use anyhow::Error;
use borsh::BorshSerialize;
use bumpalo::Bump;
use std::ops::{Deref, DerefMut};

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
    /// Wrap the inner field with a `ZGridDB`.
    pub fn new(inner: AAMap<XmerIndexableArray<ZGridSubtable<'a>>>) -> Self {
        Self(inner)
    }
    /// Moral equivalent of implementing deserialization on [`ZGridDB`],
    /// but with a memory arena to put dynamically allocated subtables into.
    pub fn deserialize(buf: &mut &[u8], arena: &'a Bump) -> Result<Self, Error> {
        let mut arr = [const { None }; 20];
        for slot in arr.iter_mut() {
            let mut arr = [const { None }; _];
            for slot in arr.iter_mut() {
                let subtable = ZGridSubtable::deserialize(buf, arena)?;
                *slot = Some(subtable)
            }
            *slot = Some(XmerIndexableArray::new(arr.map(Option::unwrap)))
        }
        Ok(ZGridDB(AAMap(arr.map(Option::unwrap))))
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
    let start_index_a = (coord_a as usize + 1).min(column_len);
    // I assume this gets optimized out
    let to_sqr_delta_a = |x: usize| {
        let delta_a = x as f64 - coord_a;
        delta_a * delta_a
    };
    let next_search_down =
        { |idx: usize| idx.checked_sub(1).map(|idx| (idx, to_sqr_delta_a(idx))) };
    let next_search_up = |idx: usize| (idx < column_len).then(|| (idx, to_sqr_delta_a(idx)));
    let mut down_search = next_search_down(start_index_a);
    let mut up_search = next_search_up(start_index_a);
    for _ in 0..column_len {
        let (idx_a, sqr_delta_a) = match (down_search, up_search) {
            (Some((down_ptr, down_sqr_delta_a)), Some((up_ptr, up_sqr_delta_a))) => {
                if down_sqr_delta_a < up_sqr_delta_a {
                    down_search = next_search_down(down_ptr);
                    (down_ptr, down_sqr_delta_a)
                } else {
                    up_search = next_search_up(up_ptr + 1);
                    (up_ptr, up_sqr_delta_a)
                }
            }
            (None, Some((up_ptr, up_sqr_delta_a))) => {
                up_search = next_search_up(up_ptr + 1);
                (up_ptr, up_sqr_delta_a)
            }
            (Some((down_ptr, down_sqr_delta_a)), None) => {
                down_search = next_search_down(down_ptr);
                (down_ptr, down_sqr_delta_a)
            }
            (None, None) => unreachable!(),
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
        // I assume this gets optimized out
        let to_abs_delta_b = |x: usize| {
            let delta_b = x as f64 - coord_b;
            delta_b.abs()
        };
        let next_search_down =
            |idx: usize| idx.checked_sub(1).map(|idx| (idx, to_abs_delta_b(idx)));
        let next_search_up = |idx: usize| (idx < row.len()).then(|| (idx, to_abs_delta_b(idx)));
        let start_index_b = (coord_b as usize + 1).min(row_len);
        let mut down_search = next_search_down(start_index_b);
        let mut up_search = next_search_up(start_index_b);
        loop {
            match (down_search, up_search) {
                (Some((down_ptr, down_abs_delta_b)), Some((up_ptr, up_abs_delta_b))) => {
                    if down_abs_delta_b < up_abs_delta_b {
                        down_search = next_search_down(down_ptr);
                        if idx_is_occupied(down_ptr) {
                            check_idx_b(down_ptr);
                            down_search = None
                        }
                    } else {
                        up_search = next_search_up(up_ptr + 1);
                        if idx_is_occupied(up_ptr) {
                            check_idx_b(up_ptr);
                            up_search = None
                        }
                    }
                    continue;
                }
                (None, Some((up_ptr, _))) => {
                    for idx_b in up_ptr..row.len() {
                        if idx_is_occupied(idx_b) {
                            check_idx_b(idx_b);
                            break;
                        }
                    }
                }
                (Some((down_ptr, _)), None) => {
                    for idx_b in (0..=down_ptr).rev() {
                        if idx_is_occupied(idx_b) {
                            check_idx_b(idx_b);
                            break;
                        }
                    }
                }
                (None, None) => (),
            };
            break;
        }
    }
    best_entry.unwrap().1
}
