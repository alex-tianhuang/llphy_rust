mod model_train_base;
mod score_type;
mod features;
pub use model_train_base::ModelTrainingBase;
pub use score_type::ScoreType;
/// The maximum residue separation that we have collected residue statistics for.
pub(crate) const MAX_XMER: usize = 40;