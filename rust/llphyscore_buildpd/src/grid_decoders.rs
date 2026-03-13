//! Module defining [`build_griddecoders`].
use anyhow::Error;
use borsh::{BorshDeserialize, BorshSerialize};
use bumpalo::Bump;
use llphyscore_core::datatypes::{Aminoacid, GridDecoder, GridDecoderPair, ModelTrainingBase, Thresholds};
use pyo3::{
    Python,
    types::{PyAnyMethods, PyDict, PyDictMethods, PyList, PyListMethods, PyTuple, PyTupleMethods},
};
use std::{fs::File, io::Read, path::{Path}};

/// Load all [`GridDecoderPair`]s from Cai's (@haocai1992) old data files,
/// serialize them, verify roundtrip integrity, and then write them
/// to an appropriate output directory.
///
/// The old data files are read from a database at `scoredb_root` and the
/// new `gridscorer.bin` files are dumped to a database at `pkg_data_root`.
pub fn build_griddecoders(scoredb_root: &Path, pkg_data_root: &Path) -> Result<(), Error> {
    Python::attach(|py| {
        let mut arena = Bump::new();
        let mut bytes = Vec::new();
        for model_train_base in [ModelTrainingBase::Human, ModelTrainingBase::HumanPDB, ModelTrainingBase::PDB] {
            arena.reset();
            let grid_decoders = load_grid_decoders(scoredb_root, model_train_base, &arena, py)?;
            for (pair_name, grid_decoder_pair) in grid_decoders {
                let outdir = pkg_data_root.join("feature_pairs").join(pair_name);
                let outpath = outdir.join(format!("{}.model.bin", model_train_base));
                bytes.clear();
                if outpath.exists() {
                    let mut file = File::open(&outpath)?;
                    file.read_to_end(&mut bytes)?;
                    let round_trip = GridDecoderPair::deserialize(&mut &*bytes)?;
                    if &round_trip != grid_decoder_pair {
                        let inpath = scoredb_root.join(format!("trained_weights.8FEATURES.{}.pkl", model_train_base));
                        return Err(Error::msg(format!(
                            "pre-existing file @ {} is not the same model data as loaded from {}",
                            outpath.display(),
                            inpath.display()
                        )));
                    }
                } else {
                    grid_decoder_pair.serialize(&mut bytes)?;
                    let round_trip = GridDecoderPair::deserialize(&mut &*bytes)?;
                    assert!(&round_trip == grid_decoder_pair);
                    let mut file = File::options().create(true).truncate(true).write(true).open(&outpath)?;
                    grid_decoder_pair.serialize(&mut file)?;
                }
            }
        }
        Ok(())
    })
}

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
fn load_grid_decoders<'a>(
    scoredb_root: &Path,
    model_train_base: ModelTrainingBase,
    arena: &'a Bump,
    py: Python,
) -> Result<&'a [(&'a str, GridDecoderPair)], Error> {
    let filename = scoredb_root.join(format!(
        "trained_weights.8FEATURES.{}.pkl",
        model_train_base
    ));
    let mut bytes = Vec::new();
    File::open(&filename)?.read_to_end(&mut bytes)?;
    let model_pickled = bytes;
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
                    item.extract::<(char, [f64; 4])>()?;
                let aa = Aminoacid::try_from(aa)?;
                thresholds_a[aa].upper = upper_bound_a;
                thresholds_a[aa].lower = lower_bound_a;
                thresholds_b[aa].upper = upper_bound_b;
                thresholds_b[aa].lower = lower_bound_b;
            }
            debug_assert!(thresholds_a.is_nan_free());
            debug_assert!(thresholds_b.is_nan_free());
            <Result<_, Error>>::Ok((
                pair_name_static,
                GridDecoderPair {
                    decoder_a: GridDecoder {
                        thresholds: thresholds_a,
                        sign: signs[0],
                    },
                    decoder_b: GridDecoder {
                        thresholds: thresholds_b,
                        sign: signs[1],
                    },
                },
            ))
        }))?,
    )
}
