//! Module of datatypes and associated serialization/deserialization types.
//!
//! Also contains computation logic for [`GridScorer::score_sequence`] and [`GridDecoder::decode`].
//!
//! [`sequences`] for sequence-related datatypes.
//! [`llphy`] for LLPhyScore-specific datatypes.
mod sequences;
pub use sequences::{AAIndex, AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str};
mod llphy;
pub use llphy::{
    FEATURE_NAMES, FeatureMatrix, GridDecoder, GridScore, GridScorer, MAX_XMER, ModelTrainingBase,
    PAIR_NAMES_AND_FEATURE_NAMES, PostProcessedFeatureMatrix, ReferenceFeatureMatrix, ScoreType,
};
