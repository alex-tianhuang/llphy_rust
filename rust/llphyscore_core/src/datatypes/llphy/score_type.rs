//! Module defining the [`ScoreType`] enum.
use clap::ValueEnum;
/// The type of score to return.
#[derive(Clone, ValueEnum)]
#[value(verbatim_doc_comment)]
pub enum ScoreType {
    /// Report raw `llphyscore` values and features.
    Raw,
    /// Report z-scores in comparison to a reference proteome.
    ZScore,
    /// Report values in as percentiles compared to a reference proteome.
    Percentile,
}