//! Module full of unpickling functions imported
//! from my private `_python` implementations.
//! 
//! It seems strange to have rust call python
//! in a binding for a python library but I
//! have no way of deserializing pickle files
//! with `numpy` objects in them.

use anyhow::Error;
use bumpalo::Bump;
use pyo3::{Python, intern, pyfunction, types::PyAnyMethods};

use crate::{datatypes::MAX_XMER, load_legacy};

#[pyfunction]
pub fn load_named_pdb_statistics() -> Result<(), Error> {
    let arena = Bump::new();
    load_legacy::load_legacy_pdb_statistics(&arena)?;
    Ok(())
}