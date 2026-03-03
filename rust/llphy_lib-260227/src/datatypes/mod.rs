mod sequences;
pub(crate) use sequences::{AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str};
mod llphy;
pub(crate) use llphy::{
    FeatureGrid, FeatureGridEntry, GridScorer, LLPhyFeature, LLPhySigns, LLPhyThresholds, LineKey,
    MAX_XMER, ModelTrainingType, PostProcessor, ScoreType,
};
