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
    AvgSdevDB, AvgSdevDBEntry, FEATURE_NAMES, FeatureMatrix, GridDecoder, GridDecoderPair,
    GridScore, GridScorer, MAX_XMER, ModelTrainingBase, PAIR_NAMES_AND_FEATURE_NAMES, PairFreqDB,
    PairFreqDBEntry, PostProcessedFeatureMatrix, PostProcessor, ReferenceFeatureMatrix, ScoreType,
    ZGridDB, ZGridDBEntry, ZGridSubtable, XmerIndexableArray, XmerSize, Thresholds
};
