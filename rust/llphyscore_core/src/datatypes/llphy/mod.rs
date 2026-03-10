//! Datatypes for LLPhyScore.
//!
//! [`features`] for containers that hold feature values.
//! [`grid_scorer`] for biophysical feature grids ([`GridScore`]s) and the models that make them ([`GridScorer`]s).
//! [`grid_decoder`] for [`GridDecoder`], which turns biophysical feature grids into sequence-level features.
//! [`consts`] for some constants.
pub use consts::{DEFAULT_FEATURES, MAX_XMER, PAIR_NAMES_AND_FEATURE_NAMES};
pub use features::{FeatureMatrix, PostProcessedFeatureMatrix, ReferenceFeatureMatrix};
pub use grid_decoder::GridDecoder;
pub use grid_scorer::{GridScore, GridScorer};
pub use model_train_base::ModelTrainingBase;
pub use score_type::ScoreType;
mod consts;
mod features;
mod grid_decoder;
mod grid_scorer;
mod model_train_base;
mod score_type;
