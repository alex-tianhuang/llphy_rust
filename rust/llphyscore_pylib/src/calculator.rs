//! Module defining [`LLPhyScoreCalculator`]
//! and its main method [`LLPhyScoreCalculator::calculate`].
use crate::{PKG_DATA_ROOT, utils::{load_grid_decoders_static, serialize_row}};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use llphyscore_core::{
    datatypes::{
        GridDecoderPair, ModelTrainingBase, PostProcessedFeatureMatrix,
        ScoreType, aa_canonical_str,
    },
    featurizer::Featurizer,
    load_pkg_data::load_post_processor,
    utils::leak_vec,
};
use pyo3::{
    Bound, PyAny, PyResult, Python,
    exceptions::PyTypeError,
    pyclass, pymethods,
    types::{PyAnyMethods, PyDict, PyDictMethods, PyList, PyListMethods, PyString, PyTypeMethods},
};
use std::{borrow::Cow, path::Path};

/// A calculator class of LLPhyScore features.
/// 
/// Dev note
/// --------
/// Describes the computation to be done but all the
/// resources for doing it (e.g. model data, grid scorers, ...)
/// are loaded lazily in [`LLPhyScoreCalculator::calculate`].
#[pyclass]
#[pyo3(frozen)]
pub struct LLPhyScoreCalculator {
    /// The data that the LLPhyScore model was trained on.
    /// 
    /// See [`ModelTrainingBase`] variants.
    #[pyo3(get)]
    model_training_base: ModelTrainingBase,
    /// Singleton of all [`GridDecoderPair`]s
    /// for the given `model_training_base`.
    grid_decoders: &'static [(&'static str, GridDecoderPair)],
}
/// A helper struct for [`LLPhyScoreCalculator::calculate`]
/// that supports passing in lists and dictionaries.
enum CalculatorInput<'a, 'py> {
    List(&'a [&'a aa_canonical_str]),
    Dict {
        keys: Vec<'a, Bound<'py, PyAny>>,
        values: &'a [&'a aa_canonical_str],
    },
}
#[pymethods]
impl LLPhyScoreCalculator {
    /// Make a new `LLPhyScoreCalculator` from the given
    /// keyword arguments.
    /// 
    /// Argument
    /// --------
    /// model_training_base : "human" | "human+PDB" | "PDB", default = "human+PDB"
    ///                     The dataset that the LLPhyScore model was trained on.
    ///                     Essentially, there are three models trained on different
    ///                     negative datasets which this argument selects for.
    #[new]
    #[pyo3(signature = (*, model_training_base = ModelTrainingBase::HumanPDB))]
    pub fn new(
        py: Python,
        model_training_base: ModelTrainingBase,
    ) -> PyResult<Self> {
        let grid_decoders = load_grid_decoders_static(py, model_training_base)?;
        Ok(Self {
            model_training_base,
            grid_decoders,
        })
    }
    /// Compute LLPhyScore features from the given `sequences`.
    /// 
    /// Arguments
    /// ---------
    /// sequences : list[str] | dict[..., str]
    ///            See section below.
    /// score_type : "raw" | "z-score" | "percentile", default = "percentile"
    ///            The type of score to be reported after computing the raw scores.
    ///            Z-scores and percentiles will be computed in reference to the
    ///            human IDRome.
    /// disable_pbar : bool, default = false
    ///              By default, some progress bars are shown while doing the computation
    ///              of each feature pair. This option turns it off.
    /// 
    /// Sequences
    /// ---------
    /// The `sequences` argument must either be:
    /// 1. A list of aminoacid strings (capitalized, no whitespace), or
    /// 2. A dictionary of arbitrary keys to aminoacid strings (capitalized, no whitespace)
    /// 
    /// In case (1), a list of {feat_name: feat_value} pairs will be returned
    /// in the same order that the sequences were inputted in.
    /// 
    /// In case (2), a dictionary of the form {key: {feat_name: feat_value}}
    /// will be returned, with each key mapping to the result that its
    /// corresponding input sequence generated.
    /// 
    /// If any strings expected to contain only aminoacids have unexpected
    /// characters, the function will error before running any computations.
    #[pyo3(signature = (sequences, *, score_type = ScoreType::Percentile, disable_pbar = false))]
    pub fn calculate<'py>(
        &self,
        py: Python<'py>,
        sequences: Bound<'py, PyAny>,
        score_type: ScoreType,
        disable_pbar: bool
    ) -> PyResult<Bound<'py, PyAny>> {
        let arena = Bump::new();
        let input = CalculatorInput::extract(&sequences, &arena)?;
        let interrupter = || py.check_signals().map_err(Error::new);
        let featurizer = Featurizer::from_parts(Path::new(PKG_DATA_ROOT), self.grid_decoders)
            .with_pbar(!disable_pbar)
            .with_interrupter(&interrupter);
        let post_processor = load_post_processor(
            Path::new(PKG_DATA_ROOT),
            score_type,
            self.model_training_base,
            &arena,
        )?;
        let matrix = featurizer.featurize(input.sequences().iter().cloned(), &arena)?;
        let matrix = post_processor.post_process(matrix, &arena)?;
        Ok(input.serialize(matrix, py)?)
    }
    /// Path to the `llphyscore` database that this calculator uses.
    #[getter]
    pub fn database_path<'py>(&self, py: Python<'py>) -> Bound<'py, PyString> {
        PyString::intern(py, PKG_DATA_ROOT)
    }
}
impl<'a, 'py> CalculatorInput<'a, 'py> {
    /// Moral equivalent of implementing [`pyo3::FromPyObject`] but
    /// strings go into a memory arena to avoid N+1 deallocations.
    pub fn extract(obj: &'a Bound<'py, PyAny>, arena: &'a Bump) -> PyResult<Self> {
        if let Ok(obj) = obj.cast::<PyList>() {
            return Self::extract_list(obj, arena);
        };
        if let Ok(obj) = obj.cast::<PyDict>() {
            return Self::extract_dict(obj, arena);
        };
        Err(PyTypeError::new_err(format!(
            "expected list or dict, got {}",
            obj.get_type().qualname()?
        )))
    }
    /// Helper for [`Self::extract`] that handles lists.
    fn extract_list(obj: &Bound<PyList>, arena: &'a Bump) -> PyResult<Self> {
        let items = arena.alloc_slice_try_fill_iter(obj.iter().map(
            |obj| -> PyResult<&aa_canonical_str> {
                let s = obj.extract::<Cow<str>>()?;
                let s = arena.alloc_str(&s);
                Ok(aa_canonical_str::from_bytes(s.as_bytes())?)
            },
        ))?;
        return Ok(CalculatorInput::List(items));
    }
    /// Helper for [`Self::extract`] that handles dictionaries.
    fn extract_dict(obj: &'a Bound<'py, PyDict>, arena: &'a Bump) -> PyResult<Self> {
        let mut keys = Vec::with_capacity_in(obj.len(), arena);
        let mut values = Vec::with_capacity_in(obj.len(), arena);
        for (key, value) in obj.iter() {
            let s = value.extract::<Cow<str>>()?;
            let s = arena.alloc_str(&s);
            let value = aa_canonical_str::from_bytes(s.as_bytes())?;
            keys.push(key);
            values.push(value)
        }
        return Ok(CalculatorInput::Dict {
            keys,
            values: leak_vec(values),
        });
    }
    /// Get a slice of sequences regardless of the input type.
    pub fn sequences(&self) -> &[&aa_canonical_str] {
        match *self {
            CalculatorInput::List(items) => items,
            CalculatorInput::Dict { values, .. } => values,
        }
    }
    /// Serialize the given matrix to the output format described
    /// in [`LLPhyScoreCalculator::calculate`].
    pub fn serialize(
        self,
        matrix: PostProcessedFeatureMatrix<'_>,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let mut feature_names = std::vec::Vec::with_capacity(matrix.row_len());
        for feat_name in matrix.feature_names() {
            feature_names.push(PyString::new(py, *feat_name))
        }
        feature_names.push(PyString::new(
            py,
            &format!("{}-feature sum", feature_names.len()),
        ));
        match self {
            Self::List(_) => {
                let items = PyList::new(
                    py,
                    matrix.rows().map(|row| {
                        serialize_row(row, &feature_names, py)
                            .expect("failed to __setitem__ on dictionary")
                    }),
                )?;
                Ok(items.into_any())
            },
            Self::Dict { keys, .. } => {
                let items = PyList::new(
                    py,
                    keys.into_iter().zip(matrix.rows()).map(|(key, row)| {
                        (key, serialize_row(row, &feature_names, py)
                            .expect("failed to __setitem__ on dictionary"))
                    }),
                )?;
                let dict = PyDict::from_sequence(&items.into_sequence().into_any())?;
                Ok(dict.into_any())
            }
        }
    }
}
