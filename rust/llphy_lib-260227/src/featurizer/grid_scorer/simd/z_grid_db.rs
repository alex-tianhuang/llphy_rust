//! Module defining substructures of [`crate::featurizer::grid_scorer::ZGridDB`]
//! using `#[portable_simd]`.
use crate::featurizer::grid_scorer::z_grid_db::lookup_thorough;
use anyhow::Error;
use borsh::{BorshDeserialize, BorshSerialize};
use bumpalo::Bump;
use std::ops::AddAssign;
use std::simd::{
    StdFloat,
    cmp::SimdPartialOrd,
    f64x2, i64x4,
    num::{SimdFloat, SimdInt},
};

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
/// but the first 3/4 of the struct is an array consisting of
/// named integers `[weight_a, weight_b, weight_total]`.
///
/// The remaining integer array slot (at index 3) is for
/// [`Self::is_occupied`], which is a boolean that checks
/// if this field contains valid weight data.
#[derive(Clone, Copy)]
pub struct ZGridDBEntry(i64x4);
/// An accumulator for computing the weighted average
/// of a bunch of [`ZGridDBEntry`]s.
///
/// Like [`ZGridDBEntry`], it is basically just three named
/// integers `[weight_a, weight_b, weight_total]`.
///
/// The `i64` at index 3 is chosen to be ignored because
/// it does not map to numeric data in [`ZGridDBEntry`].
pub struct ZGridEntrySum(i64x4);
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
    pub fn new(dbl_z_offsets: [f64; 2], row_len: usize, data: &'a [ZGridDBEntry]) -> Self {
        Self {
            dbl_z_offsets: f64x2::from_array(dbl_z_offsets),
            row_len,
            data,
        }
    }
    /// Get the entry best associated with the given `zscores`.
    ///
    /// It is assumed that the two floats represent the zscore
    /// for feature `A` and `B` respectively.
    pub fn lookup(&self, zscores: f64x2) -> &ZGridDBEntry {
        let dbl_zscores = zscores * f64x2::splat(2.0);
        if let Some(entry) = self.lookup_quick(dbl_zscores) {
            entry
        } else {
            self.lookup_thorough(dbl_zscores)
        }
    }
    /// Snap doubled zscores to the grid and
    /// see if an entry exists and is occupied there.
    fn lookup_quick(&self, dbl_zscores: f64x2) -> Option<&ZGridDBEntry> {
        // Bounds derived from Cai's `make_linekey`
        // function from the original `LLPhyScore`.
        let clamped_zscores = dbl_zscores
            .round()
            .simd_clamp(f64x2::splat(-16.0), f64x2::splat(24.0));
        let indexes = clamped_zscores - self.dbl_z_offsets;
        if indexes.simd_lt(f64x2::splat(0.0)).any() {
            return None;
        }
        let [idx_a, idx_b] = indexes.cast::<u64>().to_array();
        let idx_a = idx_a as usize;
        let idx_b = idx_b as usize;
        let row = self
            .data
            .get(self.row_len * idx_a..self.row_len * (idx_a + 1))?;
        let entry = row.get(idx_b)?;
        entry.to_occupied()
    }
    /// Find the gridpoint that minimizes the sum of squared
    /// differences to the given zscores.
    fn lookup_thorough(&self, dbl_zscores: f64x2) -> &ZGridDBEntry {
        let [coord_a, coord_b] = (dbl_zscores - self.dbl_z_offsets).to_array();
        lookup_thorough(
            self.data,
            self.row_len,
            ZGridDBEntry::to_occupied,
            coord_a,
            coord_b,
        )
    }
    /// Moral equivalent of implementing deserialize on [`ZGridSubtable`],
    /// but uses a memory arena to hold the dynamically sized subtable entries.
    pub fn deserialize(buf: &mut &[u8], arena: &'a Bump) -> Result<Self, Error> {
        let dbl_z_offsets = <[f64; 2]>::deserialize(buf)?;
        let row_len = <usize>::deserialize(buf)?;
        let data_len = <usize>::deserialize(buf)?;
        // TODO: add size check or convince myself there is not an attack vector here
        let data = arena.alloc_slice_try_fill_with(data_len, |_| {
            let array_repr = <[i64; 4]>::deserialize(buf)?;
            <Result<_, Error>>::Ok(ZGridDBEntry(i64x4::from_array(array_repr)))
        })?;
        Ok(ZGridSubtable {
            dbl_z_offsets: f64x2::from_array(dbl_z_offsets),
            row_len,
            data,
        })
    }
}
impl BorshSerialize for ZGridSubtable<'_> {
    fn serialize<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        self.dbl_z_offsets.as_array().serialize(writer)?;
        self.row_len.serialize(writer)?;
        self.data.len().serialize(writer)?;
        for cell in self.data {
            cell.0.as_array().serialize(writer)?;
        }
        Ok(())
    }
}
impl ZGridDBEntry {
    /// Get a new, occupied [`ZGridDBEntry`] from its three fields.
    pub fn new_occupied(weight_total: i64, weight_a: i64, weight_b: i64) -> Self {
        Self(i64x4::from_array([weight_a, weight_b, weight_total, 1]))
    }
    /// Make a new unoccupied [`ZGridDBEntry`].
    pub const fn unoccupied() -> Self {
        Self(i64x4::splat(0))
    }
    /// Return this entry if it is occupied.
    fn to_occupied(&self) -> Option<&Self> {
        self.is_occupied().then_some(self)
    }
    /// True if this memory location is occupied
    /// with valid/initialized weight data.
    fn is_occupied(&self) -> bool {
        self.0[3] != 0
    }
}

impl ZGridEntrySum {
    /// Get a new [`ZGridEntrySum`] with all zeroes.
    pub fn new_zeroed() -> Self {
        Self(i64x4::splat(0))
    }
    /// The sum of total weights of all entries observed.
    fn weight_total(&self) -> i64 {
        self.0[2]
    }
    /// After adding many [`ZGridDBEntry`]s to this sum,
    /// report the sum of weights `A` and `B` divided by
    /// the total weight of all observations taken.
    ///
    /// Utility method for [`crate::featurizer::GridScorer::score_sequence`].
    pub fn as_frequencies(&self) -> [f64; 2] {
        if self.weight_total() == 0 {
            [0.0; 2]
        } else {
            (self.0.extract::<0, 2>().cast::<f64>() / f64x2::splat(self.weight_total() as f64))
                .to_array()
        }
    }
}
impl AddAssign<&ZGridDBEntry> for ZGridEntrySum {
    fn add_assign(&mut self, rhs: &ZGridDBEntry) {
        self.0 += &rhs.0;
    }
}
