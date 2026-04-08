//! Crate defining a python library for computing LLPhyScore features.
//! (currently just [`LLPhyScoreCalculator`]).
use pyo3::{Bound, PyResult, pymodule, types::{PyModule, PyModuleMethods}, wrap_pyfunction};
use crate::calculator::LLPhyScoreCalculator;
mod calculator;
mod utils;
mod transform_scores;
use transform_scores::transform_scores as transform_scores_fn;
const PKG_DATA_ROOT: &'static str = env!("PKG_DATA_ROOT");
#[pymodule]
fn llphyscore(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<LLPhyScoreCalculator>()?;
    m.add_function(wrap_pyfunction!(transform_scores_fn, m)?)?;
    Ok(())
}