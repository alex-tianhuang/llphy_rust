//! LLPhyScore-specific datatypes.
mod pdb_statistics;
use std::fmt::Display;

use anyhow::Error;
pub(crate) use pdb_statistics::{FeatureGrid, FeatureGridEntry, GridScorer, LineKey};
mod model;
pub(crate) use model::{LLPhyFeature, LLPhySigns, LLPhyThresholds};
mod post_processor;
use clap::ValueEnum;
pub(crate) use post_processor::PostProcessor;
use pyo3::{FromPyObject, PyErr};

/// The maximum residue separation that we have collected residue statistics for.
pub(crate) const MAX_XMER: usize = 40;
/// Whether the phase separation model was trained on a negative
/// dataset consisting of the human proteome, the PDB, or both.
#[derive(Copy, Clone, ValueEnum)]
#[value(verbatim_doc_comment)]
pub enum ModelTrainingType {
    /// Model was trained on a human proteome negative dataset.
    Human,
    #[value(name = "PDB")]
    /// Model was trained on the PDB for its negative dataset.
    PDB,
    #[value(name = "human+PDB")]
    /// Model was trained on the PDB + human proteome negative dataset.
    HumanPDB,
}
impl<'a, 'py> FromPyObject<'a, 'py> for ModelTrainingType {
    type Error = PyErr;
    fn extract(obj: pyo3::Borrowed<'a, 'py, pyo3::PyAny>) -> Result<Self, Self::Error> {
        let s = obj.extract()?;
        Ok(ModelTrainingType::from_str(s, false).map_err(Error::msg)?)
    }
}
impl Display for ModelTrainingType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Self::Human => "human",
            Self::PDB => "PDB",
            Self::HumanPDB => "human+PDB",
        };
        f.write_str(s)
    }
}

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
impl<'a, 'py> FromPyObject<'a, 'py> for ScoreType {
    type Error = PyErr;
    fn extract(obj: pyo3::Borrowed<'a, 'py, pyo3::PyAny>) -> Result<Self, Self::Error> {
        let s = obj.extract()?;
        Ok(ScoreType::from_str(s, false).map_err(Error::msg)?)
    }
}
