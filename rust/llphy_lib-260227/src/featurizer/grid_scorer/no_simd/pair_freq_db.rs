//! Module defining [`PairFreqDB`] and [`PairFreqEntrySum`]
//! without the use of the `#[portable_simd]` feature.
use crate::datatypes::{AAMap, MAX_XMER};
use std::ops::{AddAssign, Deref, DerefMut};

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
/// Dev note
/// --------
/// The reason this contains data for two features instead
/// of separating them out nicely is because a pair of zscores
/// are required to index into [`super::ZGridDB`] and so it is a
/// win for cache locality to get them loaded in one struct.
pub struct PairFreqDBEntry {
    weight_a: f64,
    weight_b: f64,
    total_a: f64,
    total_b: f64
}
/// An accumulator struct for taking the sum of many [`PairFreqDBEntry`]s.
pub struct PairFreqEntrySum {
    weight_a: f64,
    weight_b: f64,
    total_a: f64,
    total_b: f64
}
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

impl PairFreqDBEntry {
    /// Shorthand for iterating over fields.
    fn iter(&self) -> impl Iterator<Item = &f64> {
        // SAFETY: this struct is `#[repr(C)]` so I know exactly what it looks like.
        unsafe{ &*(self as *const Self as *const [f64; 4]) }.iter()
    }
    /// Get a new [`PairFreqDBEntry`] that is filled with `f64::NAN`.
    const fn new_nan_filled() -> Self {
        Self {
            weight_a: f64::NAN,
            weight_b: f64::NAN,
            total_a: f64::NAN,
            total_b: f64::NAN,
        }
    }
    /// False if there are `f64::NAN`s in any field.
    pub fn is_nan_free(&self) -> bool {
        self.iter()
            .all(|x| !x.is_nan())
    }
    /// Set weights and total for feature `A`.
    pub fn set_a(&mut self, weight: f64, total: f64) {
        self.weight_a = weight;
        self.total_a = total;
    }
    /// Set weights and total for feature `B`.
    pub fn set_b(&mut self, weight: f64, total: f64) {
        self.weight_b = weight;
        self.total_b = total;
    }
    
}
impl PairFreqEntrySum {
    /// Get a new [`PairFreqEntrySum`] that is filled with `0.0_f64`.
    pub fn new_zeroed() -> Self {
        Self {
            weight_a: 0.0,
            weight_b: 0.0,
            total_a: 0.0,
            total_b: 0.0,
        }
    }
    /// Helper method for [`crate::featurizer::GridScorer::score_sequence`].
    /// 
    /// Equivalent to:
    /// ```
    /// [self.weight_a / self.total_a, self.weight_b / self.total_b]
    /// ```
    pub fn as_frequencies(&self) -> [f64; 2] {
        [self.weight_a / self.total_a, self.weight_b / self.total_b]
    }
}
impl AddAssign<&PairFreqDBEntry> for PairFreqEntrySum {
    fn add_assign(&mut self, rhs: &PairFreqDBEntry) {
        self.weight_a += rhs.weight_a;
        self.weight_b += rhs.weight_b;
        self.total_a += rhs.total_a;
        self.total_b += rhs.total_b;
    }
}