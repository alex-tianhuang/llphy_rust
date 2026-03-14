use pyo3::{Bound, PyResult, pymodule, types::{PyModule, PyModuleMethods}};

use crate::calculator::LLPhyScoreCalculator;

mod calculator;
mod utils;
const PKG_DATA_ROOT: &'static str = env!("PKG_DATA_ROOT");
#[pymodule]
fn llphyscore(m: Bound<PyModule>) -> PyResult<()> {
    m.add_class::<LLPhyScoreCalculator>()
}