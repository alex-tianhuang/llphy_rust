//! Module defining [`load_grid_decoders`].
use crate::{
    datatypes::{Aminoacid, ModelTrainingBase},
    features::{GridDecoder, Thresholds},
    load_pkg_data::read_archive_file,
};
use anyhow::Error;
use bumpalo::Bump;
use pyo3::{
    Python,
    types::{PyAnyMethods, PyDict, PyDictMethods, PyList, PyListMethods, PyTuple, PyTupleMethods},
};
use std::path::PathBuf;

/// Signs of each feature pair,
/// ripped from the original python library this package is from.
/// 
/// Pairs are organized like [`sign_for_feat_a`, `sign_for_feat_b`].
const FEATURE_SIGNS: &'static [(&'static str, [i8; 2])] = &[
    ("S2.SUMPI", [1, 1]),
    ("S3.WATER.V2", [1, -1]),
    ("S4.SSPRED", [-1, 1]),
    ("S5.DISO", [1, 1]),
    ("S6.CHARGE.V2", [1, -1]),
    ("S7.ELECHB.V2", [1, 1]),
    ("S8.CationPi.V2", [-1, -1]),
    ("S9.LARKS.V2", [1, -1]),
];

/// Get all [`GridDecoder`]s for the given [`ModelTrainingType`].
///
/// Dev note
/// --------
/// Uses `Python` because one of the pickled objects is a numpy
/// array which `serde_pickle` does not natively deserialize.
///
/// Probably one day I will remove this useless `numpy` dependency,
/// but not today.
pub fn load_grid_decoders<'a>(
    model_train_base: ModelTrainingBase,
    arena: &'a Bump,
    py: Python,
) -> Result<&'a [(&'a str, [GridDecoder; 2])], Error> {
    let filename = PathBuf::from(format!(
        "trained_weights.8FEATURES.{}.pkl",
        model_train_base
    ));
    let model_pickled = read_archive_file(&filename)?;
    let locals = PyDict::new(py);
    locals.set_item("pickled", model_pickled)?;
    py.run(c"import pickle\nunpickled = [(name, [(aa, block.tolist()) for aa, block in entry.items()]) for name, entry in pickle.loads(pickled).items()]", None, Some(&locals))?;
    let unpickled = locals.get_item("unpickled")?.unwrap();
    let unpickled_list = unpickled.cast_into::<PyList>().unwrap();
    Ok(
        &*arena.alloc_slice_try_fill_iter(unpickled_list.iter().map(|obj| {
            let item = obj.cast_into::<PyTuple>().unwrap();
            let [pair_name, sublist] = item.as_slice() else {
                unreachable!()
            };
            let pair_name = pair_name.extract::<&str>().unwrap();
            let (pair_name_static, signs) = *FEATURE_SIGNS
                .iter()
                .find(|(known_name, _)| *known_name == pair_name)
                .ok_or_else(|| {
                    Error::msg(format!(
                        "unrecognized pair name while unloading model: {}",
                        pair_name
                    ))
                })?;
            let mut thresholds_a = Thresholds::new_nan_filled();
            let mut thresholds_b = Thresholds::new_nan_filled();
            for item in sublist.cast::<PyList>().unwrap().iter() {
                let (aa, [upper_bound_a, lower_bound_a, upper_bound_b, lower_bound_b]) =
                    item.extract::<(Aminoacid, [f64; 4])>()?;
                thresholds_a[aa].upper = upper_bound_a;
                thresholds_a[aa].lower = lower_bound_a;
                thresholds_b[aa].upper = upper_bound_b;
                thresholds_b[aa].lower = lower_bound_b;
            }
            debug_assert!(thresholds_a.is_nan_free());
            debug_assert!(thresholds_b.is_nan_free());
            <Result<_, Error>>::Ok((
                pair_name_static,
                [
                    GridDecoder {
                        thresholds: thresholds_a,
                        sign: signs[0],
                    },
                    GridDecoder {
                        thresholds: thresholds_b,
                        sign: signs[1],
                    },
                ],
            ))
        }))?,
    )
}
