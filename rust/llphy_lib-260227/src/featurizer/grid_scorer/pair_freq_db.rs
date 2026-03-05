//! Module defining [`PairFreqDB`],
//! a `(aa_x, gap_length, xy_orientation, aa_y)`-indexable
//! collection of [`PairFreqDBEntry`]s.
//! 
//! (see [`PairFreqDB`] for what `xy_orientation` means).
use std::ops::{Deref, DerefMut};
use crate::{datatypes::{AAMap, MAX_XMER}, featurizer::grid_scorer::PairFreqDBEntry};


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
