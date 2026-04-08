use crate::{PKG_DATA_ROOT, utils::serialize_row};
use bumpalo::{Bump, collections::Vec};
use llphyscore_core::{
    datatypes::{FeatureMatrix, ModelTrainingBase, PostProcessedFeatureMatrix, ScoreType},
    load_pkg_data::load_post_processor,
    utils::leak_vec,
};
use pyo3::{
    Bound, PyAny, PyErr, PyResult, Python,
    exceptions::{PyTypeError, PyValueError},
    pyfunction,
    types::{PyAnyMethods, PyDict, PyDictMethods, PyList, PyListMethods, PyString, PyTypeMethods},
};
use std::{borrow::Cow, path::Path};

/// Take raw LLPhyScore values yielded from [`LLPhyScoreCalculator`]
/// and convert them into percentiles or z-scores.
///
/// Arguments
/// ---------
/// data : list[dict[str, int]] | dict[..., dict[str, int]]
///      A list whose items or a dictionary whose values are
///      dictionaries of raw (integer) feature values. If they
///      are floats this function will raise a type error.
/// score_type : "z-score" | "percentile"
///            The type of score to be reported after computing the raw scores.
///            Z-scores and percentiles will be computed in reference to the
///            human IDRome.
/// model_training_base : "human" | "human+PDB" | "PDB", default = "human+PDB"
///                     The dataset that the LLPhyScore model was trained on.
///                     Essentially, there are three models trained on different
///                     negative datasets which this argument selects for.
#[pyfunction]
#[pyo3(signature = (data, *, score_type, model_training_base = ModelTrainingBase::HumanPDB))]
pub fn transform_scores<'py>(
    py: Python<'py>,
    data: Bound<'py, PyAny>,
    score_type: ScoreType,
    model_training_base: ModelTrainingBase,
) -> PyResult<Bound<'py, PyAny>> {
    if matches!(score_type, ScoreType::Raw) {
        return Err(PyValueError::new_err(
            r#"`score_type` must be "z-score" or "percentile""#,
        ));
    }
    let arena = Bump::new();
    let input = TransformScoresInput::extract(&data, &arena)?;
    let post_processor = load_post_processor(
        Path::new(PKG_DATA_ROOT),
        score_type,
        model_training_base,
        &arena,
    )?;
    let (matrix, serializer) = input.into_parts()?;
    let matrix = post_processor.post_process(matrix, &arena)?;
    serializer.serialize(matrix, py)
}

/// A helper struct for [`transform_scores`]
/// that supports passing in lists and dictionaries.
enum TransformScoresInput<'a, 'py> {
    List {
        // Due to current limitations, "8-feature sum" must be the last feature here.
        feature_names: &'a [&'a str],
        values: &'a mut [i64],
    },
    Dict {
        // Due to current limitations, "8-feature sum" must be the last feature here.
        feature_names: &'a [&'a str],
        keys: &'a [Bound<'py, PyAny>],
        values: &'a mut [i64],
    },
}
/// A helper struct for [`transform_scores`]
/// that supports serializing results as a list or dictionary.
enum Serializer<'a, 'py> {
    List,
    Dict { keys: &'a [Bound<'py, PyAny>] },
}
/// An error message to display to the user that "8-feature sum" is a necessary feature of the `transform_score` function.
const ERROR_NEED_8FEATURESUM: &str = r#""8-feature sum" was not found as a feature in this dictionary. Due to current limitations, you must include the "8-feature sum" values when calling `transform_score`."#;
impl<'a, 'py> TransformScoresInput<'a, 'py> {
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
    fn extract_list(obj: &Bound<'py, PyList>, arena: &'a Bump) -> PyResult<Self> {
        if obj.len() == 0 {
            return Ok(TransformScoresInput::List {
                feature_names: &mut [],
                values: &mut [],
            });
        }
        let first_item = obj
            .get_item(0)
            .and_then(|obj| obj.cast_into::<PyDict>().map_err(PyErr::from))?;
        let feature_names = Self::extract_feature_names(&first_item, arena)?;
        let num_rows = obj.len();
        let mut values = Vec::with_capacity_in(num_rows * feature_names.len(), &arena);
        for obj in obj.iter() {
            let obs_len = obj.len()?;
            let expected_len = feature_names.len();
            if obs_len != expected_len {
                return Err(PyValueError::new_err(format!(
                    "expected list of feature dictionaries with identical keys, but found two dictionaries with {} keys and another with {} keys",
                    obs_len, expected_len
                )));
            }
            for fname in feature_names {
                let feat_value =
                    PyAnyMethods::get_item(&obj, fname).and_then(|obj| obj.extract::<i64>())?;
                values.push(feat_value)
            }
        }
        return Ok(Self::List {
            feature_names,
            values: leak_vec(values),
        });
    }
    /// Helper for [`Self::extract`] that handles dictionaries.
    fn extract_dict(obj: &Bound<'py, PyDict>, arena: &'a Bump) -> PyResult<Self> {
        if obj.len() == 0 {
            return Ok(TransformScoresInput::Dict {
                feature_names: &mut [],
                keys: &[],
                values: &mut [],
            });
        }
        let num_rows = obj.len();
        let mut iter = obj.iter();
        let (first_key, first_value) = iter
            .next()
            .expect("checked length > 0 but no items yielded");
        let feature_names = Self::extract_feature_names(first_value.cast::<PyDict>()?, arena)?;

        let mut keys = Vec::with_capacity_in(num_rows, arena);
        let mut values = Vec::with_capacity_in(num_rows * feature_names.len(), arena);
        keys.push(first_key);
        for fname in feature_names {
            let feat_value =
                PyAnyMethods::get_item(&first_value, fname).and_then(|obj| obj.extract::<i64>())?;
            values.push(feat_value)
        }
        for (key, value) in obj.iter() {
            keys.push(key);
            for fname in feature_names {
                let feat_value =
                    PyAnyMethods::get_item(&value, fname).and_then(|obj| obj.extract::<i64>())?;
                values.push(feat_value)
            }
        }
        return Ok(Self::Dict {
            feature_names,
            keys: leak_vec(keys),
            values: leak_vec(values),
        });
    }
    /// Helper for [`Self::extract`] that gets feature names from a dictionary.
    ///
    /// As a hotfix, requires that "8-feature sum" is present in this dictionary.
    /// If it is not, the function will fail (because post-processor is expecting it anyway).
    fn extract_feature_names(obj: &Bound<'py, PyDict>, arena: &'a Bump) -> PyResult<&'a [&'a str]> {
        let feature_names =
            arena.alloc_slice_try_fill_iter(obj.iter().map(|(key, _)| -> PyResult<&'a str> {
                let s = key.extract::<Cow<str>>()?;
                let s = arena.alloc_str(&s);
                Ok(s)
            }))?;
        let n = feature_names
            .iter()
            .position(|s| *s == "8-feature sum")
            .ok_or_else(|| PyValueError::new_err(ERROR_NEED_8FEATURESUM))?;
        feature_names[n..].rotate_left(1);
        Ok(feature_names)
    }
    /// Convert this input into a feature matrix and serializer.
    ///
    /// As a hotfix, requires that "8-feature sum" is present in the
    /// list of feature names, and then truncates the list to not include it.
    ///
    /// If "8-feature sum" is not the last feature name,
    /// this function returns an error.
    pub fn into_parts(self) -> PyResult<(FeatureMatrix<'a, i64>, Serializer<'a, 'py>)> {
        match self {
            Self::List {
                feature_names,
                values,
            } => {
                let [feature_names @ .., last] = feature_names else {
                    return Err(PyValueError::new_err(ERROR_NEED_8FEATURESUM));
                };
                if *last != "8-feature sum" {
                    return Err(PyValueError::new_err(ERROR_NEED_8FEATURESUM));
                }
                Ok((
                    FeatureMatrix {
                        feature_names,
                        data: values,
                    },
                    Serializer::List,
                ))
            }
            Self::Dict {
                feature_names,
                keys,
                values,
            } => {
                let [feature_names @ .., last] = feature_names else {
                    return Err(PyValueError::new_err(ERROR_NEED_8FEATURESUM));
                };
                if *last != "8-feature sum" {
                    return Err(PyValueError::new_err(ERROR_NEED_8FEATURESUM));
                }
                Ok((
                    FeatureMatrix {
                        feature_names,
                        data: values,
                    },
                    Serializer::Dict { keys },
                ))
            }
        }
    }
}
impl<'a, 'py> Serializer<'a, 'py> {
    /// Serialize a feature matrix as a list or dictionary.
    pub fn serialize(
        self,
        matrix: PostProcessedFeatureMatrix<'a>,
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
            Self::List => {
                let items = PyList::new(
                    py,
                    matrix.rows().map(|row| {
                        serialize_row(row, &feature_names, py)
                            .expect("failed to __setitem__ on dictionary")
                    }),
                )?;
                Ok(items.into_any())
            }
            Self::Dict { keys, .. } => {
                let items = PyList::new(
                    py,
                    keys.into_iter().zip(matrix.rows()).map(|(key, row)| {
                        (
                            key,
                            serialize_row(row, &feature_names, py)
                                .expect("failed to __setitem__ on dictionary"),
                        )
                    }),
                )?;
                let dict = PyDict::from_sequence(&items.into_sequence().into_any())?;
                Ok(dict.into_any())
            }
        }
    }
}
