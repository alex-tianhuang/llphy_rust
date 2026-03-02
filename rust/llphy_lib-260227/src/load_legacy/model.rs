use anyhow::Error;
use bumpalo::Bump;
use pyo3::{
    Python,
    types::{PyAnyMethods, PyDict, PyList, PyListMethods, PyTuple},
};
use std::{array, path::PathBuf};

use crate::{
    datatypes::{AAMap, Aminoacid, LLPhyFeature, LLPhySigns, LLPhyThresholds, ModelTrainingType},
    load_legacy::read_archive_file,
};

pub fn load_legacy_model_and_feature_order<'a>(
    model_train_base: ModelTrainingType,
    arena: &'a Bump,
    py: Python,
) -> Result<&'a mut [(&'a str, LLPhyFeature)], Error> {
    let filename = PathBuf::from(format!(
        "trained_weights.8FEATURES.{}.pkl",
        model_train_base
    ));
    let model_pickled = read_archive_file(&filename)?;
    let locals = PyDict::new(py);
    locals.set_item("pickled", model_pickled)?;
    py.run(c"import pickle\nunpickled = [(name, [(aa, block.tolist()) for aa, block in entry.items()]) for name, entry in pickle.loads(pickled).items()]", None, Some(&locals))?;
    let unpickled = locals.get_item("unpickled").unwrap();
    let unpickled_list = unpickled.cast_into::<PyList>().unwrap();
    arena.alloc_slice_try_fill_iter(unpickled_list.iter().map(|obj| {
        let item = obj.cast_into::<PyTuple>().unwrap();
        let grid_name = item.get_item(0).unwrap();
        let grid_name = grid_name.extract::<&str>().unwrap();
        let (_, signs) = *FEATURE_SIGNS
            .iter()
            .find(|(known_name, _)| *known_name == grid_name)
            .ok_or_else(|| {
                Error::msg(format!(
                    "unrecognized grid name while unloading model: {}",
                    grid_name
                ))
            })?;
        let sublist = item.get_item(1)?.cast_into::<PyList>().unwrap();
        let mut thresholds = AAMap(array::from_fn(|_| LLPhyThresholds {
            sr: [f64::NAN; 2],
            lr: [f64::NAN; 2],
        }));
        for item in sublist.iter() {
            let (aa, [sr_upper, sr_lower, lr_upper, lr_lower]) =
                item.extract::<(Aminoacid, [f64; 4])>()?;
            let entry = &mut thresholds[aa];
            entry.sr = [sr_upper, sr_lower];
            entry.lr = [lr_upper, lr_lower];
        }
        debug_assert!(thresholds.values().all(|entry| {
            [entry.sr, entry.lr]
                .into_iter()
                .flatten()
                .all(|slot| !slot.is_nan())
        }));
        <Result<_, Error>>::Ok((
            arena.alloc_str(grid_name) as &str,
            LLPhyFeature { thresholds, signs },
        ))
    }))
}
const FEATURE_SIGNS: &'static [(&'static str, LLPhySigns)] = &[
    ("S2.SUMPI", LLPhySigns { sr: 1, lr: 1 }),
    ("S3.WATER.V2", LLPhySigns { sr: 1, lr: -1 }),
    ("S4.SSPRED", LLPhySigns { sr: -1, lr: 1 }),
    ("S5.DISO", LLPhySigns { sr: 1, lr: 1 }),
    ("S6.CHARGE.V2", LLPhySigns { sr: 1, lr: -1 }),
    ("S7.ELECHB.V2", LLPhySigns { sr: 1, lr: 1 }),
    ("S8.CationPi.V2", LLPhySigns { sr: -1, lr: -1 }),
    ("S9.LARKS.V2", LLPhySigns { sr: 1, lr: -1 }),
];
