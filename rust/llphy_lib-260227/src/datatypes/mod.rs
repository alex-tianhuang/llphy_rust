mod sequences;
pub(crate) use sequences::{AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str};
mod llphy;
pub(crate) use llphy::{
    GridScorer, LLPhyFeature, PostProcessor, FeatureGridEntry, FeatureGrid, ScoreType, MAX_XMER, LineKey, ModelTrainingType, LLPhySigns, LLPhyThresholds
};
