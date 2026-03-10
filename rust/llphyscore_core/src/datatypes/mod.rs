//! Module of datatypes, their serialization/deserialization methods,
//! and some basic computations/conversions that can be performed on them.
//! 
//! [`sequences`] for sequence-related datatypes.
//! [`llphy`] for LLPhyScore-specific datatypes.
mod sequences;
pub use sequences::{AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str, AAIndex};
mod llphy;
pub use llphy::{
    FeatureMatrix, MAX_XMER, PostProcessedFeatureMatrix, ReferenceFeatureMatrix,
    DEFAULT_FEATURES, GridScore, GridDecoder, GridScorer, PAIR_NAMES_AND_FEATURE_NAMES
};
