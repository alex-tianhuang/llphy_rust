mod sequences;
pub(crate) use sequences::{AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str};
mod llphy;
pub(crate) use llphy::{
    FeatureMatrix, MAX_XMER, ModelTrainingBase, PostProcessedFeatureMatrix, ReferenceFeatureMatrix,
    ScoreType, DEFAULT_FEATURES
};
