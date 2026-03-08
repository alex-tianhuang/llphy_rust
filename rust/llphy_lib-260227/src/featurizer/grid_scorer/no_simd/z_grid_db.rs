//! Module defining substructures of [`crate::featurizer::grid_scorer::ZGridDB`]
//! without using `#[portable_simd]`.
use crate::featurizer::grid_scorer::z_grid_db::{KNOWN_MAX_DATA_LEN, lookup_thorough};
use anyhow::Error;
use borsh::{BorshDeserialize, BorshSerialize};
use bumpalo::Bump;
use bytesize::ByteSize;
use std::ops::AddAssign;

/// A `(zscore_a, zscore_b)`-indexable collection of weights
/// for features `a` and `b`.
///
/// The range of possible zscore pairs is assumed to look like a
/// dense square, such that the entries can be arranged in a matrix
/// without wasting too much space. In practice, about half
/// the space is filled with empty entries, but it hopefully
/// gives a big performance boost over binary search.
#[derive(PartialEq)]
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
#[derive(Clone, Copy, PartialEq)]
pub struct ZGridDBEntry {
    weight_a: i64,
    weight_b: i64,
    weight_total: i64,
}
/// An accumulator for computing the weighted average
/// of a bunch of [`ZGridDBEntry`]s.
pub struct ZGridEntrySum {
    weight_a: i64,
    weight_b: i64,
    weight_total: i64,
}
impl<'a> ZGridSubtable<'a> {
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
        let idx_a = dbl_zscores[0].round() - self.dbl_z_offsets[0];
        let idx_b = dbl_zscores[1].round() - self.dbl_z_offsets[1];
        if idx_a < 0.0 || idx_b < 0.0 {
            return None;
        }
        let idx_a = idx_a as usize;
        let idx_b = idx_b as usize;
        if idx_a * self.row_len >= self.data.len() || idx_b >= self.row_len {
            return None;
        }
        let entry = &self.data[idx_a * self.row_len + idx_b];
        entry.as_ref()
    }
    /// Find the gridpoint that minimizes the sum of squared
    /// differences to the given zscores.
    fn lookup_thorough(&self, dbl_zscores: [f64; 2]) -> &ZGridDBEntry {
        let coord_a = dbl_zscores[0] - self.dbl_z_offsets[0];
        let coord_b = dbl_zscores[1] - self.dbl_z_offsets[1];
        lookup_thorough(self.data, self.row_len, Option::as_ref, coord_a, coord_b)
    }
    /// Moral equivalent of implementing deserialize on [`ZGridSubtable`],
    /// but uses a memory arena to hold the dynamically sized subtable entries.
    pub fn deserialize(buf: &mut &[u8], arena: &'a Bump) -> Result<Self, Error> {
        let dbl_z_offsets = <[f64; 2]>::deserialize(buf)?;
        let row_len = <usize>::deserialize(buf)?;
        let data_len = <usize>::deserialize(buf)?;
        if data_len > KNOWN_MAX_DATA_LEN {
            let cell_size = size_of::<Option<ZGridDBEntry>>() as u64;
            let asked_for = ByteSize::b(cell_size * data_len as u64);
            let expected = ByteSize::b(cell_size * KNOWN_MAX_DATA_LEN as u64);
            return Err(Error::msg(format!(
                "[ZGridSubtable::deserialize] expected at most {:?} to be allocated, but asking for {:?}",
                expected, asked_for
            )));
        }
        let data = arena.alloc_slice_try_fill_with(data_len, |_| {
            let [weight_a, weight_b, weight_total, flag] = <[i64; 4]>::deserialize(buf)?;
            if flag != 0 {
                <Result<_, Error>>::Ok(Some(ZGridDBEntry {
                    weight_a,
                    weight_b,
                    weight_total,
                }))
            } else {
                <Result<_, Error>>::Ok(None)
            }
        })?;
        Ok(ZGridSubtable::new(dbl_z_offsets, row_len, data))
    }
}
impl BorshSerialize for ZGridSubtable<'_> {
    fn serialize<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        self.dbl_z_offsets.serialize(writer)?;
        self.row_len.serialize(writer)?;
        self.data.len().serialize(writer)?;
        for cell in self.data {
            let array_repr = match cell {
                Some(entry) => [entry.weight_a, entry.weight_b, entry.weight_total, 1],
                None => [0; 4],
            };
            array_repr.serialize(writer)?;
        }
        Ok(())
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
        self.weight_a += rhs.weight_a;
        self.weight_b += rhs.weight_b;
        self.weight_total += rhs.weight_total;
    }
}
