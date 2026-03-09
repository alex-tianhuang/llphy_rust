//! Module defining the [`ScoreType`] enum.
use anyhow::Error;
use clap::ValueEnum;
use pyo3::{PyErr, FromPyObject};
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
