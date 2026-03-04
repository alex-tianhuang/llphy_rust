//! Module defining [`load_post_processor`].
use anyhow::Error;
use bumpalo::Bump;

use crate::{
    datatypes::{ModelTrainingBase, ScoreType},
    load_pkg_data::load_reference_scores,
    post_processor::PostProcessor,
};

/// Load a [`PostProcessor`] that will produce the requested
/// `score_type`, based on the feature distributions for human sequences
/// according to the model trained on `model_training_base`.
pub fn load_post_processor(
    score_type: ScoreType,
    model_train_base: ModelTrainingBase,
    arena: &Bump,
) -> Result<PostProcessor<'_>, Error> {
    match score_type {
        ScoreType::Raw => Ok(PostProcessor::Raw),
        ScoreType::ZScore => load_reference_scores(model_train_base, arena)
            .map(|ref_scores| PostProcessor::new_zscore(ref_scores, arena)),
        ScoreType::Percentile => {
            load_reference_scores(model_train_base, arena).map(PostProcessor::new_percentile)
        }
    }
}
