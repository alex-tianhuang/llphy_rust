//! Module defining [`PairFreqDB`].
use crate::datatypes::{AAMap, MAX_XMER};
use std::ops::{AddAssign, Deref, DerefMut};
use std::simd::{f64x4, f64x2};

/// The number of `gap_length`-indexable tables in a [`PairFreqDB`].
///
/// Dev note
/// --------
/// This was probably a bug in the original `LLPhyScore`
/// to include the entry for `0` gap length,
/// but it's being officially adopted by this program.
const PAIR_FREQ_DB_LEN: usize = MAX_XMER + 1;
/// A struct containing [`PairFreqDBEntry`] for each
/// `(aa_x, gap_length, xy_orientation, aa_y)` key.
///
/// - `xy_orientation` could theoretically be represented by a `bool`
///   as it is describing whether `aa_y` is N-terminal or C-terminal
///   of `aa_x`. In practice I just do both so it doesn't matter.
/// - `gap_length` can be any number in `0..=MAX_XMER`.
///   The inner array is one longer than an `xmer`-indexable array.
pub struct PairFreqDB(AAMap<[PairFreqSubtable; PAIR_FREQ_DB_LEN]>);
/// A substructure of [`PairFreqDB`] that organizes entries based
/// on the `xy_orientation` (whether residue `y` is N-terminal or
/// C-terminal of `x`).
pub struct PairFreqSubtable {
    // If `aa_y` is N-terminal to `aa_x`
    // index into this field with `aa_y`.
    pub n_terminal_mapping: AAMap<PairFreqDBEntry>,
    // If `aa_y` is C-terminal to `aa_x`,
    // index into this field with `aa_y`.
    pub c_terminal_mapping: AAMap<PairFreqDBEntry>,
}
/// Weights for each `(aa_x, gap_length, xy_orientation, aa_y)` key.
/// 
/// The field names are not visible in the SIMD representation,
/// but it is essentially an array consisting of named floats
/// [`weight_a`], [`weight_b`], [`total_a`], and [`total_b`].
///
/// [`weight_a`]: Self::weight_a
/// [`weight_b`]: Self::weight_b
/// [`total_a`]: Self::total_a
/// [`total_b`]: Self::total_b
/// 
/// Dev note
/// --------
/// The reason this contains data for two features instead
/// of separating them out nicely is because a pair of zscores
/// are required to index into [`super::ZGridDB`] and so it is a
/// win for cache locality to get them loaded in one struct.
pub struct PairFreqDBEntry(f64x4);
impl Deref for PairFreqDB {
    type Target = AAMap<[PairFreqSubtable; PAIR_FREQ_DB_LEN]>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for PairFreqDB {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl PairFreqDB {
    /// Get a new [`PairFreqDB`] struct that is filled with `f64::NAN`.
    pub fn new_nan_filled() -> Self {
        PairFreqDB(AAMap(
            [const {
                [const {
                    PairFreqSubtable {
                        c_terminal_mapping: AAMap([const { PairFreqDBEntry::new_nan_filled() }; 20]),
                        n_terminal_mapping: AAMap([const { PairFreqDBEntry::new_nan_filled() }; 20]),
                    }
                }; PAIR_FREQ_DB_LEN]
            }; 20],
        ))
    }
    /// False if there are `f64::NAN`s in any entry.
    pub fn is_nan_free(&self) -> bool {
        !self.0.values().flatten().any(|subtable| {
            [subtable.c_terminal_mapping.values(), subtable.n_terminal_mapping.values()]
                .into_iter()
                .flatten()
                .any(|e| !e.is_nan_free())
        })
    }
}

/// Shorthand for associating each index of the
/// array of [`PairFreqDBEntry`] with a named field.
macro_rules! impl_getters {
    ($([$index:literal, $field:ident]),*) => {
        impl PairFreqDBEntry {
            $(pub fn $field(&self) -> f64 {
                self.0.as_array()[$index]
            })*
        }
    };
}
impl_getters!([0, weight_a], [1, weight_b], [2, total_a], [3, total_b]);
impl PairFreqDBEntry {
    /// Get a new [`PairFreqDBEntry`] that is filled with `f64::NAN`.
    const fn new_nan_filled() -> Self {
        Self(f64x4::from_array([f64::NAN; 4]))
    }
    /// False if there are `f64::NAN`s in any field.
    pub fn is_nan_free(&self) -> bool {
        self.0.as_array()
            .into_iter()
            .all(|x| !x.is_nan())
    }
    /// Get a new [`PairFreqDBEntry`] that is filled with `0.0_f64`.
    pub fn new_zeroed() -> Self {
        Self(f64x4::from_array([0.0; 4]))
    }
    /// Set weights and total for feature `A`.
    pub fn set_a(&mut self, weight: f64, total: f64) {
        let this = self.0.as_mut_array();
        this[0] = weight;
        this[2] = total;
    }
    /// Set weights and total for feature `B`.
    pub fn set_b(&mut self, weight: f64, total: f64) {
        let this = self.0.as_mut_array();
        this[1] = weight;
        this[3] = total;
    }
    /// Equivalent to `f64x2::from_array([self.weight_a() / self.total_a(), self.weight_b() / self.total_b()])`.
    pub fn as_frequencies(&self) -> f64x2 {
        self.0.extract::<0, 2>() / self.0.extract::<2, 2>()
    }
}
impl AddAssign<&Self> for PairFreqDBEntry {
    fn add_assign(&mut self, rhs: &Self) {
        self.0 += &rhs.0
    }
}