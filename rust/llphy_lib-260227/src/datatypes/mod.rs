mod sequences;
pub(crate) use sequences::{AAIndex, AAMap, AMINOACIDS, Aminoacid, FastaEntry, aa_canonical_str};
mod llphy;
pub(crate) use llphy::{
    GridScorer, LLPhyFeatureOld, PostProcessor, ResScoresEntryOld, ResScoresOld, ScoreType,
};
