//! Module defining [`PairFreqDB`],
//! a `(aa_x, gap_length, xy_orientation, aa_y)`-indexable
//! collection of [`PairFreqDBEntry`]s.
//!
//! (see [`PairFreqDB`] for what `xy_orientation` means).
use anyhow::Error;
use borsh::{BorshDeserialize, BorshSerialize};
use std::{
    mem::MaybeUninit,
    ops::{Deref, DerefMut},
    ptr::addr_of_mut,
};

use crate::{
    datatypes::{AAMap, MAX_XMER},
    datatypes::llphy::grid_scorer::PairFreqDBEntry,
};

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
#[derive(BorshSerialize, PartialEq)]
#[repr(transparent)]
pub struct PairFreqDB(AAMap<[PairFreqSubtable; PAIR_FREQ_DB_LEN]>);
/// A substructure of [`PairFreqDB`] that organizes entries based
/// on the `xy_orientation` (whether residue `y` is N-terminal or
/// C-terminal of `x`).
#[derive(BorshSerialize, PartialEq)]
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
    /// Given space for a new [`PairFreqDB`] table,
    /// initialize all its entries to be `f64::NAN`s.
    pub fn init_with_nans(this: &mut MaybeUninit<Self>) -> &mut Self {
        // SAFETY: PairFreqDB has #[repr(transparent)] AAMap([_; PAIR_FREQ_DB_LEN])
        //         AAMap has #[repr(transparent)] over [_; 20]
        //         no other initialization guarantees
        let table = unsafe {
            &mut *this
                .as_mut_ptr()
                .cast::<[[MaybeUninit<PairFreqSubtable>; PAIR_FREQ_DB_LEN]; 20]>()
        };
        for arr in table {
            for subtable in arr {
                let subtable = subtable.as_mut_ptr();
                let n_terminal_mapping: *mut AAMap<PairFreqDBEntry>;
                let c_terminal_mapping: *mut AAMap<PairFreqDBEntry>;

                // SAFETY: `subtable` points to a valid memory region with size/align
                //         appropriate `PairFreqSubtable` (see all the #[repr(transparent)]s)
                unsafe {
                    n_terminal_mapping = addr_of_mut!((*subtable).n_terminal_mapping);
                    c_terminal_mapping = addr_of_mut!((*subtable).n_terminal_mapping);
                }
                for mapping in [n_terminal_mapping, c_terminal_mapping] {
                    // SAFETY: AAMap has #[repr(trasparent)] over [_; 20]
                    let mapping =
                        unsafe { &mut *mapping.cast::<[MaybeUninit<PairFreqDBEntry>; 20]>() };
                    for slot in mapping {
                        slot.write(PairFreqDBEntry::new_nan_filled());
                    }
                }
            }
        }
        // SAFETY: just initialized all slots
        unsafe { this.assume_init_mut() }
    }
    /// False if there are `f64::NAN`s in any entry.
    pub fn is_nan_free(&self) -> bool {
        !self.0.values().flatten().any(|subtable| {
            [
                subtable.c_terminal_mapping.values(),
                subtable.n_terminal_mapping.values(),
            ]
            .into_iter()
            .flatten()
            .any(|e| !e.is_nan_free())
        })
    }
    /// Moral equivalent of implementing deserialization on [`PairFreqDB`],
    /// but does it by reference so that don't have to put it on the stack.
    pub fn deserialize_into(this: &mut MaybeUninit<Self>, buf: &mut &[u8]) -> Result<(), Error> {
        let this = Self::init_with_nans(this);
        for arr in this.0.values_mut() {
            for subtable in arr {
                for mapping in [&mut subtable.n_terminal_mapping, &mut subtable.c_terminal_mapping] {
                    for slot in mapping.values_mut() {
                        let data = PairFreqDBEntry::deserialize(buf)?;
                        *slot = data;
                    }
                }
            }
        }
        debug_assert!(this.is_nan_free());
        Ok(())
    }
}
