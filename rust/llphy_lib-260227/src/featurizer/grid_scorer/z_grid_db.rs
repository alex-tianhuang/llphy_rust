//! Module defining [`ZGridDB`].
use crate::datatypes::AAMap;
use crate::featurizer::grid_scorer::xmer::XmerIndexableArray;
use std::ops::{AddAssign, Deref, DerefMut};
use std::simd::{f64x2, i64x4, StdFloat, num::SimdFloat};

/// A struct of tables indexable by `(aa, xmer)` keys,
/// where each subtable is 2D `zscore`-indexable.
/// See [`ZGridSubtable`] and [`ZGridDBEntry`].
///
/// The subtable is expected to be indexed by `zscore_a`
/// and then `zscore_b`, NOT `zscore_b` then `zscore_a`.
pub struct ZGridDB<'a>(AAMap<XmerIndexableArray<ZGridSubtable<'a>>>);
/// A `(zscore_a, zscore_b)`-indexable collection of weights
/// for features `a` and `b`.
///
/// The range of possible zscore pairs is assumed to look like a
/// dense square, such that the entries can be arranged in a matrix
/// without wasting too much space. In practice, about half
/// the space is filled with empty entries, but it hopefully
/// gives a big performance boost over binary search.
pub struct ZGridSubtable<'a> {
    /// Corresponds to double the minimum z-score
    /// for features `A` then feature `B`
    /// (as index 0 and 1 respectively).
    dbl_z_offsets: f64x2,
    /// The length of a row after indexing by `zscore_a`,
    /// or equivalently the number of half-integers
    /// in the range of valid `zscore_b`s.
    row_len: usize,
    data: &'a [ZGridDBEntry],
}
/// A possibly occupied slot for weights for features `a` and `b`,
/// for some `(aa, xmer, zscore_a, zscore_b)` tuples.
///
/// The field names are not visible in the SIMD representation,
/// but 3/4 of the struct is an array consisting of named integers
/// [`weight_total`], [`weight_a`], and [`weight_b`].
///
/// The remaining integer array slot is for [`Self::is_occupied`],
/// which is a boolean that checks if this field contains
/// valid weight data.
///
/// [`weight_total`]: Self::weight_total
/// [`weight_a`]: Self::weight_a
/// [`weight_b`]: Self::weight_b
#[derive(Clone, Copy)]
pub struct ZGridDBEntry(i64x4);
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
}
impl<'a> ZGridSubtable<'a> {
    /// Make a new [`ZGridSubtable`] with no content.
    pub const fn placeholder() -> Self {
        Self {
            dbl_z_offsets: f64x2::from_array([f64::NAN; 2]),
            row_len: 0,
            data: &[],
        }
    }
    /// Make a new [`ZGridSubtable`] with the given data.
    pub fn new(dbl_z_offsets: f64x2, row_len: usize, data: &'a [ZGridDBEntry]) -> Self {
        Self {
            dbl_z_offsets,
            row_len,
            data,
        }
    }
    /// Get the entry best associated with the given `zscores`.
    ///
    /// It is assumed that the two floats represent the zscore
    /// for feature `A` and `B` respectively.
    pub fn lookup<const D: bool>(&self, zscores: f64x2) -> &ZGridDBEntry {
        let dbl_zscores = zscores * f64x2::splat(2.0);
        if let Some(entry) = self.lookup_quick(dbl_zscores) {
            entry
        } else {
            self.lookup_thorough(dbl_zscores)
        }
    }
    /// The number of half-integers in the range of valid `zscore_a`s.
    fn column_len(&self) -> usize {
        self.data.len() / self.row_len
    }
    /// Snap doubled zscores to the grid and
    /// see if an entry exists and is occupied there.
    fn lookup_quick(&self, dbl_zscores: f64x2) -> Option<&ZGridDBEntry> {
        // Bounds derived from Cai's `make_linekey`
        // function from the original `LLPhyScore`.
        let clamped_zscores = dbl_zscores.round().simd_clamp(f64x2::splat(-16.0), f64x2::splat(24.0));
        let [idx_a, idx_b] = (clamped_zscores - self.dbl_z_offsets).to_array();
        if idx_a < 0.0 || idx_b < 0.0 {
            return None;
        }
        let idx_a = idx_a as usize;
        let idx_b = idx_b as usize;
        let row = self
            .data
            .get(self.row_len * idx_a..self.row_len * (idx_a + 1))?;
        let entry = row.get(idx_b)?;
        entry.is_occupied().then_some(entry)
    }
    /// Find the gridpoint that minimizes the sum of squared
    /// differences to the given zscores.
    fn lookup_thorough(&self, dbl_zscores: f64x2) -> &ZGridDBEntry {
        // z-star = (z * 2 - self.dbl_z_offsets);
        let [z_star_a, z_star_b] = (dbl_zscores - self.dbl_z_offsets).to_array();
        let mut best_entry: Option<(f64, &ZGridDBEntry)> = None;
        let start_index_a = (z_star_a as usize + 1).min(self.column_len());
        // I assume this gets optimized out
        let to_sqr_delta_a = |x: usize| {
            let delta_a = x as f64 - z_star_a;
            delta_a * delta_a
        };
        let next_search_down =
            { |idx: usize| idx.checked_sub(1).map(|idx| (idx, to_sqr_delta_a(idx))) };
        let next_search_up =
            |idx: usize| (idx < self.column_len()).then(|| (idx, to_sqr_delta_a(idx)));
        let mut down_search = next_search_down(start_index_a);
        let mut up_search = next_search_up(start_index_a);
        debug_assert!(self.column_len() > 0);
        for _ in 0..self.column_len() {
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
            let row = &self.data[idx_a * self.row_len..(idx_a + 1) * self.row_len];
            debug_assert!(!row.is_empty());
            let mut check_idx_b = |idx_b: usize| {
                let delta_b = idx_b as f64 - z_star_b;
                let sqr_delta_b = delta_b * delta_b;
                let sqr_delta = sqr_delta_a + sqr_delta_b;
                let new_best: bool;
                if let Some(entry) = best_entry {
                    new_best = sqr_delta < entry.0;
                } else {
                    new_best = true;
                }
                if new_best {
                    best_entry = Some((sqr_delta, &row[idx_b]))
                }
            };
            // I assume this gets optimized out
            let to_abs_delta_b = |x: usize| {
                let delta_b = x as f64 - z_star_b;
                delta_b.abs()
            };
            let next_search_down =
                |idx: usize| idx.checked_sub(1).map(|idx| (idx, to_abs_delta_b(idx)));
            let next_search_up = |idx: usize| (idx < row.len()).then(|| (idx, to_abs_delta_b(idx)));
            let start_index_b = (z_star_b as usize + 1).min(self.row_len);
            let mut down_search = next_search_down(start_index_b);
            let mut up_search = next_search_up(start_index_b);
            loop {
                match (down_search, up_search) {
                    (Some((down_ptr, down_abs_delta_b)), Some((up_ptr, up_abs_delta_b))) => {
                        if down_abs_delta_b < up_abs_delta_b {
                            down_search = next_search_down(down_ptr);
                            if row[down_ptr].is_occupied() {
                                check_idx_b(down_ptr);
                                down_search = None
                            }
                        } else {
                            up_search = next_search_up(up_ptr + 1);
                            if row[up_ptr].is_occupied() {
                                check_idx_b(up_ptr);
                                up_search = None
                            }
                        }
                        continue;
                    }
                    (None, Some((up_ptr, _))) => {
                        for idx_b in up_ptr..row.len() {
                            if row[idx_b].is_occupied() {
                                check_idx_b(idx_b);
                                break;
                            }
                        }
                    }
                    (Some((down_ptr, _)), None) => {
                        for idx_b in (0..=down_ptr).rev() {
                            if row[idx_b].is_occupied() {
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
}
/// Shorthand for associating each index of the
/// array of [`ZGridDBEntry`] with a named field.
macro_rules! impl_getters {
    ($([$index:literal, $field:ident]),*) => {
        impl ZGridDBEntry {
            $(pub fn $field(&self) -> i64 {
                self.0.as_array()[$index]
            })*
        }
    };
}
impl_getters!([1, weight_total], [2, weight_a], [3, weight_b]);
impl ZGridDBEntry {
    /// True if this memory location is occupied
    /// with valid/initialized weight data.
    pub fn is_occupied(&self) -> bool {
        self.0.as_array()[0] != 0
    }
    /// Get a new [`ZGridDBEntry`] with all zeroes.
    pub fn new_zeroed() -> Self {
        Self(i64x4::from_array([0; 4]))
    }
    /// Get a new, occupied [`ZGridDBEntry`] from its three fields.
    pub fn new_occupied(weight_total: i64, weight_a: i64, weight_b: i64) -> Self {
        Self(i64x4::from_array([1, weight_total, weight_a, weight_b]))
    }
    /// Utility method for [`super::GridScorer::score_sequence`].
    pub fn freq_a(&self) -> f64 {
        if self.weight_total() == 0 {
            0.0
        } else {
            self.weight_a() as f64 / self.weight_total() as f64
        }
    }
    /// Utility method for [`super::GridScorer::score_sequence`].
    pub fn freq_b(&self) -> f64 {
        if self.weight_total() == 0 {
            0.0
        } else {
            self.weight_b() as f64 / self.weight_total() as f64
        }
    }
}
impl AddAssign<&Self> for ZGridDBEntry {
    fn add_assign(&mut self, rhs: &Self) {
        self.0 += &rhs.0;
    }
}
