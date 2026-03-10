//! Datatypes for LLPhyScore.
//! 
//! [`features`] for containers that hold feature values.
//! [`grid_scorer`] for biophysical feature grids ([`GridScore`]s) and the models that make them ([`GridScorer`]s).
//! [`grid_decoder`] for [`GridDecoder`], which turns biophysical feature grids into sequence-level features.
//! [`consts`] for some constants.
pub use features::{FeatureMatrix, PostProcessedFeatureMatrix, ReferenceFeatureMatrix};
pub use grid_scorer::{GridScorer, GridScore};
pub use grid_decoder::GridDecoder;
pub use consts::{DEFAULT_FEATURES, MAX_XMER, PAIR_NAMES_AND_FEATURE_NAMES};
mod features;
mod grid_scorer;
mod grid_decoder;
mod consts;
// pub use model_train_base::ModelTrainingBase;
// pub use score_type::ScoreType;
