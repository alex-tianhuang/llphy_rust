//! Module of datatypes and associated serialization/deserialization types.
//! 
//! Also contains computation logic for [`GridScorer::score_sequence`] and [`GridDecoder::decode`].
//! 
//! [`sequences`] for sequence-related datatypes.
//! [`llphy`] for LLPhyScore-specific datatypes.
mod sequences;
pub use sequences::{AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str, AAIndex};
mod llphy;
pub use llphy::{
    FeatureMatrix, MAX_XMER, PostProcessedFeatureMatrix, ReferenceFeatureMatrix,
    FEATURE_NAMES, GridScore, GridDecoder, GridScorer, PAIR_NAMES_AND_FEATURE_NAMES
};
