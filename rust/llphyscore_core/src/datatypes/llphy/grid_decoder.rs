//! Module defining [`GridDecoder`],
//! a struct that takes a [`GridScore`], a residue-level feature grid
//! and turns it into an integer feature value using some thresholding rules.
//! 
//! [`GridScore`]: crate::datatypes::llphy::GridScore
use anyhow::Error;
use borsh::{BorshDeserialize, BorshSerialize};
use bumpalo::Bump;
use crate::datatypes::AAMap;
use std::{mem::{MaybeUninit}, ops::{Deref, DerefMut}, ptr::addr_of_mut};

/// Two [`GridDecoder`]s with a deserialize implementation.
#[derive(BorshSerialize, PartialEq)]
pub struct GridDecoderPair {
    pub decoder_a: GridDecoder,
    pub decoder_b: GridDecoder
}
/// Takes a [`GridScore`], a residue-level feature grid
/// and turns it into a integer-valued feature using some thresholding rules.
/// 
/// [`GridScore`]: crate::datatypes::llphy::GridScore
#[derive(BorshSerialize, BorshDeserialize, PartialEq)]
pub struct GridDecoder {
    pub sign: i8,
    pub thresholds: Thresholds,
}
impl GridDecoderPair {
    /// Moral equivalent of implementing deserialization on [`GridDecoderPair`],
    /// but does it in a memory arena and returns a reference to it.
    pub fn deserialize<'a>(buf: &mut &[u8], arena: &'a Bump) -> Result<&'a Self, Error> {
        let this = arena.alloc_with(<MaybeUninit<Self>>::uninit);
        let target = this.as_mut_ptr();
        // SAFETY: target is a valid memory addr for a `GridDecoderPair`,
        //         therefore so is the calculated pointer to `decoder_a`.
        let decoder_a = unsafe { &mut *addr_of_mut!((*target).decoder_a).cast::<MaybeUninit<GridDecoder>>()};
        // SAFETY: target is a valid memory addr for a `GridDecoderPair`,
        //         therefore so is the calculated pointer to `decoder_b`.
        let decoder_b = unsafe { &mut *addr_of_mut!((*target).decoder_b).cast::<MaybeUninit<GridDecoder>>()};
        for target in [decoder_a, decoder_b] {
            target.write(GridDecoder::deserialize(buf)?);
        }
        // SAFETY: just initialized all fields
        Ok(unsafe { this.assume_init_ref() })
    }
}
impl GridDecoder {
    /// Compute a sequence-level feature
    /// from the given grid scores.
    pub fn decode(&self, grid_score: &AAMap<&[f64]>) -> i64 {
        let mut score = 0;
        for (t, &subarr) in self.thresholds.values().zip(grid_score.values()) {
            score += t.score_sites(subarr)
        }
        score * self.sign as i64
    }
}
/// A collection of [`ThresholdPair`]s for each aminoacid.
#[derive(BorshSerialize, BorshDeserialize, PartialEq)]
pub struct Thresholds(AAMap<ThresholdPair>);
/// An upper and lower grid-score threshold
/// that converts a grid-score to an integer value.
#[derive(BorshSerialize, BorshDeserialize, PartialEq)]
pub struct ThresholdPair {
    pub upper: f64,
    pub lower: f64,
}
impl Deref for Thresholds {
    type Target = AAMap<ThresholdPair>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for Thresholds {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl Thresholds {
    /// Get a new [`Thresholds`] struct that is filled with `f64::NAN`.
    pub fn new_nan_filled() -> Self {
        Self(AAMap(
            [const {
                ThresholdPair {
                    lower: f64::NAN,
                    upper: f64::NAN,
                }
            }; 20],
        ))
    }
    /// True if there are `f64::NAN`s in any slot of the thresholds.
    pub fn is_nan_free(&self) -> bool {
        !self
            .0
            .0
            .iter()
            .any(|t| t.upper.is_nan() || t.lower.is_nan())
    }
}
impl ThresholdPair {
    /// Get an integer value from a single grid score.
    fn score_site(&self, grid_score: f64) -> i64 {
        let Self { lower, upper } = *self;
        (grid_score > upper) as i64 - (grid_score < lower) as i64
    }
    /// Sum the values of several grid scores.
    pub fn score_sites(&self, grid_scores: &[f64]) -> i64 {
        grid_scores.iter().map(|x| self.score_site(*x)).sum::<i64>()
    }
}
