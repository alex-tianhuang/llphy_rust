//! Datatypes for LLPhyScore.
//!
//! [`features`] for containers that hold feature values.
//! [`grid_scorer`] for biophysical feature grids ([`GridScore`]s) and the models that make them ([`GridScorer`]s).
//! [`grid_decoder`] for [`GridDecoder`] / [`GridDecoderPair`], which turns biophysical feature grids into integer-valued features.
//! [`consts`] for some constants.
pub use consts::{FEATURE_NAMES, MAX_XMER, PAIR_NAMES_AND_FEATURE_NAMES};
pub use features::{FeatureMatrix, PostProcessedFeatureMatrix, ReferenceFeatureMatrix};
pub use grid_decoder::{GridDecoder, GridDecoderPair};
pub use grid_scorer::{
    AvgSdevDB, AvgSdevDBEntry, GridScore, GridScorer, PairFreqDB, PairFreqDBEntry,
    XmerIndexableArray, XmerSize, ZGridDB, ZGridDBEntry, ZGridSubtable,
};
pub use model_train_base::ModelTrainingBase;
pub use post_processor::{PostProcessor, ScoreType};
mod consts;
mod features;
mod grid_decoder;
mod grid_scorer;
mod model_train_base;
mod post_processor;
