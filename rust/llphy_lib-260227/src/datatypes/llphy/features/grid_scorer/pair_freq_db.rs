//! Module defining [`PairFreqDB`].
use crate::datatypes::{AAMap, MAX_XMER};
use std::ops::{Deref, DerefMut};

/// The number of `gap_length`-indexable tables in a [`PairFreqDB`].
///
/// Dev note
/// --------
/// This was probably a bug in the original `LLPhyScore`
/// to include the entry for `0` gap length,
/// but it's being officially adopted by this program.
const PAIR_FREQ_DB_LEN: usize = MAX_XMER + 1;
/// A struct containing [`PairFreqDBEntry`] for each
/// `(aa_x, gap_length, aa_y)` key.
///
/// `gap_length` can be any number in `0..=MAX_XMER`.
/// The inner array is one longer than an `xmer`-indexable array.
pub struct PairFreqDB(AAMap<[AAMap<PairFreqDBEntry>; PAIR_FREQ_DB_LEN]>);
/// Weights for each `(aa_x, gap_length, aa_y)` key.
/// 
/// Dev note
/// --------
/// The reason this contains data for two features instead
/// of separating them out nicely is because a pair of zscores
/// are required to index into [`super::ZGridDB`] and so it is a
/// win for cache locality to get them loaded in one struct.
pub struct PairFreqDBEntry {
    pub weight_a: f64,
    pub total_a: f64,
    pub weight_b: f64,
    pub total_b: f64,
}
impl Deref for PairFreqDB {
    type Target = AAMap<[AAMap<PairFreqDBEntry>; PAIR_FREQ_DB_LEN]>;
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
                [const { AAMap([const { PairFreqDBEntry::new_nan_filled() }; 20]) };
                    PAIR_FREQ_DB_LEN]
            }; 20],
        ))
    }
    #[cfg(debug_assertions)]
    /// False if there are `f64::NAN`s in any entry.
    pub fn is_nan_free(&self) -> bool {
        !self
            .0
            .values()
            .flatten()
            .any(|subtable| subtable.values().any(|e| e.is_nan_free()))
    }
}
impl PairFreqDBEntry {
    /// Get a new [`PairFreqDBEntry`] that is filled with `f64::NAN`.
    const fn new_nan_filled() -> Self {
        Self {
            weight_a: f64::NAN,
            total_a: f64::NAN,
            weight_b: f64::NAN,
            total_b: f64::NAN,
        }
    }
    #[cfg(debug_assertions)]
    /// False if there are `f64::NAN`s in any field.
    pub fn is_nan_free(&self) -> bool {
        [self.weight_a, self.weight_b, self.total_a, self.total_b]
            .into_iter()
            .all(|x| !x.is_nan())
    }
}
