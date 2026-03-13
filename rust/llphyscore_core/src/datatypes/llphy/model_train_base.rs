//! Module defining the [`ModelTrainingBase`] enum.
use std::fmt::Display;
use clap::ValueEnum;
/// Whether the phase separation model was trained on a negative
/// dataset consisting of the human proteome, the PDB, or both.
#[derive(Copy, Clone, ValueEnum)]
#[value(verbatim_doc_comment)]
pub enum ModelTrainingBase {
    /// Model was trained on a human proteome negative dataset.
    Human,
    #[value(name = "PDB")]
    /// Model was trained on the PDB for its negative dataset.
    PDB,
    #[value(name = "human+PDB")]
    /// Model was trained on the PDB + human proteome negative dataset.
    HumanPDB,
}
impl Display for ModelTrainingBase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Self::Human => "human",
            Self::PDB => "PDB",
            Self::HumanPDB => "human+PDB",
        };
        f.write_str(s)
    }
}
#[cfg(feature = "pyo3")]
impl<'a, 'py> pyo3::FromPyObject<'a, 'py> for ModelTrainingBase {
    type Error = pyo3::PyErr;
    fn extract(obj: pyo3::Borrowed<'a, 'py, pyo3::PyAny>) -> Result<Self, Self::Error> {
        let s = obj.extract()?;
        Ok(ModelTrainingBase::from_str(s, false).map_err(anyhow::Error::msg)?)
    }
}