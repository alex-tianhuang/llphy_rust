//! Module defining [`GridScorer`] and [`GridScore`].
//!
//! Structs that turn a sequence into a residue-level feature grids.
use crate::datatypes::AAMap;
use anyhow::Error;
use borsh::BorshSerialize;
use bumpalo::Bump;
use std::{mem::MaybeUninit, ptr::addr_of_mut};
pub use xmer::XmerIndexableArray;
mod avg_sdev_db;
pub use avg_sdev_db::AvgSdevDB;
mod pair_freq_db;
mod xmer;
mod z_grid_db;
pub use pair_freq_db::PairFreqDB;
pub use z_grid_db::ZGridDB;
#[cfg(feature = "simd")]
mod simd;
#[cfg(feature = "simd")]
pub use simd::{AvgSdevDBEntry, PairFreqDBEntry, PairFreqEntrySum, ZGridEntrySum, ZGridSubtable};
#[cfg(not(feature = "simd"))]
mod no_simd;
#[cfg(not(feature = "simd"))]
pub use no_simd::{
    AvgSdevDBEntry, PairFreqDBEntry, PairFreqEntrySum, ZGridEntrySum, ZGridSubtable,
};
/// A struct that contains all the necessary data to
/// make biophysical feature grids ([`GridScore`]s)
/// from sequences.
/// 
/// Dev note
/// --------
/// The computation logic is in the [`crate::featurizer`] module.
/// This struct is basically just a bunch of lookup tables.
#[derive(PartialEq, BorshSerialize)]
pub struct GridScorer<'a> {
    pub pair_freqs: PairFreqDB,
    pub avg_sdevs: AvgSdevDB,
    pub z_grid: ZGridDB<'a>,
}
/// A collection of biophysical feature scores for
/// each residue in a sequence.
pub struct GridScore<'a> {
    pub feature_a_scores: AAMap<&'a [f64]>,
    pub feature_b_scores: AAMap<&'a [f64]>,
}
impl<'a> GridScorer<'a> {
    /// Moral equivalent of implementing deserialization on [`GridScorer`],
    /// but does it in a memory arena and returns a reference to it.
    ///
    /// Dev note
    /// --------
    /// Memory arena is for:
    /// 1. `GridScorer` is 1MB (huge) and causes segfaults if I put it on the stack.
    /// 2. `ZGridDB` has a dynamically allocated size but don't want to use global allocator for that
    pub fn deserialize(buf: &mut &[u8], arena: &'a Bump) -> Result<&'a Self, Error> {
        let grid_scorer = arena.alloc_with(<MaybeUninit<GridScorer>>::uninit);
        let target = grid_scorer.as_mut_ptr();
        // SAFETY: target is a valid memory addr for a `GridScorer`,
        //         therefore so is the calculated pointer to `pair_freqs`.
        let pair_freqs =
            unsafe { &mut *addr_of_mut!((*target).pair_freqs).cast::<MaybeUninit<PairFreqDB>>() };
        PairFreqDB::deserialize_into(pair_freqs, buf)?;
        // SAFETY: target is a valid memory addr for a `GridScorer`,
        //         therefore so is the calculated pointer to `avg_sdevs`.
        let avg_sdevs =
            unsafe { &mut *addr_of_mut!((*target).avg_sdevs).cast::<MaybeUninit<AvgSdevDB>>() };
        AvgSdevDB::deserialize_into(avg_sdevs, buf)?;
        // SAFETY: target is a valid memory addr for a `GridScorer`,
        //         therefore so is the calculated pointer to `z_grid`.
        let z_grid = unsafe { &mut *addr_of_mut!((*target).z_grid).cast::<MaybeUninit<ZGridDB>>() };
        ZGridDB::deserialize_into(z_grid, buf, arena)?;
        // SAFETY: just initialized all fields above
        Ok(unsafe { grid_scorer.assume_init_ref() })
    }
}

