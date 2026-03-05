//! Module defining [`ZGridDB`] and [`ZGridEntrySum`]
//! without the use of the `#[portable_simd]` feature.
use crate::featurizer::grid_scorer::xmer::XmerIndexableArray;
use crate::{datatypes::AAMap, featurizer::grid_scorer::z_grid_db::lookup_thorough};
use std::ops::{AddAssign, Deref, DerefMut};

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
    dbl_z_offsets: [f64; 2],
    /// The length of a row after indexing by `zscore_a`,
    /// or equivalently the number of half-integers
    /// in the range of valid `zscore_b`s.
    row_len: usize,
    data: &'a [Option<ZGridDBEntry>],
}
/// An occupied slot for weights for features `a` and `b`,
/// for some `(aa, xmer, zscore_a, zscore_b)` tuples.
#[derive(Clone, Copy)]
pub struct ZGridDBEntry {
    weight_total: i64,
    weight_a: i64,
    weight_b: i64,
}
/// An accumulator for computing the weighted average
/// of a bunch of [`ZGridDBEntry`]s.
pub struct ZGridEntrySum {
    weight_total: i64,
    weight_a: i64,
    weight_b: i64,
}
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
            dbl_z_offsets: [f64::NAN; 2],
            row_len: 0,
            data: &[],
        }
    }
    /// Make a new [`ZGridSubtable`] with the given data.
    pub fn new(dbl_z_offsets: [f64; 2], row_len: usize, data: &'a [Option<ZGridDBEntry>]) -> Self {
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
    pub fn lookup(&self, zscores: [f64; 2]) -> &ZGridDBEntry {
        let dbl_zscores = zscores.map(|x| x * 2.0);
        if let Some(entry) = self.lookup_quick(dbl_zscores) {
            entry
        } else {
            self.lookup_thorough(dbl_zscores)
        }
    }
    /// Snap doubled zscores to the grid and
    /// see if an entry exists and is occupied there.
    fn lookup_quick(&self, dbl_zscores: [f64; 2]) -> Option<&ZGridDBEntry> {
        // Bounds derived from Cai's `make_linekey`
        // function from the original `LLPhyScore`.
        let clamped_zscores = dbl_zscores.map(|x| x.round().clamp(-16.0, 24.0));
        let idx_a = clamped_zscores[0] - self.dbl_z_offsets[0];
        let idx_b = clamped_zscores[1] - self.dbl_z_offsets[1];
        if idx_a < 0.0 || idx_b < 0.0 {
            return None;
        }
        let idx_a = idx_a as usize;
        let idx_b = idx_b as usize;
        let row = self
            .data
            .get(self.row_len * idx_a..self.row_len * (idx_a + 1))?;
        let entry = row.get(idx_b)?;
        entry.as_ref()
    }
    /// Find the gridpoint that minimizes the sum of squared
    /// differences to the given zscores.
    fn lookup_thorough(&self, dbl_zscores: [f64; 2]) -> &ZGridDBEntry {
        let coord_a = dbl_zscores[0] - self.dbl_z_offsets[0];
        let coord_b = dbl_zscores[1] - self.dbl_z_offsets[1];
        lookup_thorough(self.data, self.row_len, Option::as_ref, coord_a, coord_b)
    }
}
impl ZGridDBEntry {
    /// Get a new, occupied [`ZGridDBEntry`] from its three fields.
    /// 
    /// This returns an option so that the [`crate::load_pkg_data::load_grid_scorer`]
    /// function looks a little nicer between simd and no-simd.
    pub fn new_occupied(weight_total: i64, weight_a: i64, weight_b: i64) -> Option<Self> {
        Some(Self {
            weight_total,
            weight_a,
            weight_b,
        })
    }
    /// Make a new unoccupied [`ZGridDBEntry`].
    /// 
    /// This returns an option so that the [`crate::load_pkg_data::load_grid_scorer`]
    /// function looks a little nicer between simd and no-simd.
    pub const fn unoccupied() -> Option<Self> {
        None
    }
}
impl ZGridEntrySum {
    /// Get a new [`ZGridEntrySum`] with all zeroes.
    pub fn new_zeroed() -> Self {
        Self {
            weight_total: 0,
            weight_a: 0,
            weight_b: 0,
        }
    }
    /// After adding many [`ZGridDBEntry`]s to this sum,
    /// report the sum of weights `A` and `B` divided by
    /// the total weight of all observations taken.
    ///
    /// Utility method for [`crate::featurizer::GridScorer::score_sequence`].
    pub fn as_frequencies(&self) -> [f64; 2] {
        let Self {
            weight_total,
            weight_a,
            weight_b,
        } = *self;
        if weight_total > 0 {
            [
                weight_a as f64 / weight_total as f64,
                weight_b as f64 / weight_total as f64,
            ]
        } else {
            [0.0; 2]
        }
    }
}
impl AddAssign<&ZGridDBEntry> for ZGridEntrySum {
    fn add_assign(&mut self, rhs: &ZGridDBEntry) {
        self.weight_total += rhs.weight_total;
        self.weight_a += rhs.weight_a;
        self.weight_b += rhs.weight_b;
    }
}
