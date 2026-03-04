pub use grid_decoder::GridDecoder;
pub use grid_scorer::{GridScore, GridScorer};
pub use thresholds::{ThresholdPair, Thresholds};
mod grid_decoder;
pub(crate) mod grid_scorer;
mod thresholds;

// pub fn load_legacy_model_and_feature_order<'a>(
//     model_train_base: ModelTrainingType,
//     arena: &'a Bump,
//     py: Python,
// ) -> Result<&'a mut [(&'a str, LLPhyFeature)], Error> {
//     let filename = PathBuf::from(format!(
//         "trained_weights.8FEATURES.{}.pkl",
//         model_train_base
//     ));
//     let model_pickled = read_archive_file(&filename)?;
//     let locals = PyDict::new(py);
//     locals.set_item("pickled", model_pickled)?;
//     py.run(c"import pickle\nunpickled = [(name, [(aa, block.tolist()) for aa, block in entry.items()]) for name, entry in pickle.loads(pickled).items()]", None, Some(&locals))?;
//     let unpickled = locals.get_item("unpickled").unwrap();
//     let unpickled_list = unpickled.cast_into::<PyList>().unwrap();
//     arena.alloc_slice_try_fill_iter(unpickled_list.iter().map(|obj| {
//         let item = obj.cast_into::<PyTuple>().unwrap();
//         let grid_name = item.get_item(0).unwrap();
//         let grid_name = grid_name.extract::<&str>().unwrap();
//         let (_, signs) = *FEATURE_SIGNS
//             .iter()
//             .find(|(known_name, _)| *known_name == grid_name)
//             .ok_or_else(|| {
//                 Error::msg(format!(
//                     "unrecognized grid name while unloading model: {}",
//                     grid_name
//                 ))
//             })?;
//         let sublist = item.get_item(1)?.cast_into::<PyList>().unwrap();
//         let mut thresholds = AAMap(array::from_fn(|_| LLPhyThresholds {
//             sr: [f64::NAN; 2],
//             lr: [f64::NAN; 2],
//         }));
//         for item in sublist.iter() {
//             let (aa, [sr_upper, sr_lower, lr_upper, lr_lower]) =
//                 item.extract::<(Aminoacid, [f64; 4])>()?;
//             let entry = &mut thresholds[aa];
//             entry.sr = [sr_upper, sr_lower];
//             entry.lr = [lr_upper, lr_lower];
//         }
//         debug_assert!(thresholds.values().all(|entry| {
//             [entry.sr, entry.lr]
//                 .into_iter()
//                 .flatten()
//                 .all(|slot| !slot.is_nan())
//         }));
//         <Result<_, Error>>::Ok((
//             arena.alloc_str(grid_name) as &str,
//             LLPhyFeature { thresholds, signs },
//         ))
//     }))
// }
