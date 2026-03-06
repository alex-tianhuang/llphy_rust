mod legacy;
// use crate::{
//     datatypes::{ModelTrainingBase, ReferenceFeatureMatrix, ScoreType},
//     featurizer::{GridDecoder, GridScorer},
//     io::read_file_into_global,
//     leak_vec, post_processor::PostProcessor,
// };
// use anyhow::Error;
// use borsh::BorshDeserialize;
// use bumpalo::{Bump, collections::Vec};
// use pyo3::Python;
// use std::{fs, path::Path};
pub use legacy::{load_grid_decoders, load_grid_scorer, load_post_processor};
// /// Replacing [`legacy::load_grid_decoders`].
// pub fn load_grid_decoders<'a>(
//     model_train_base: ModelTrainingBase,
//     arena: &'a Bump,
//     py: Python,
// ) -> Result<&'a [(&'a str, [GridDecoder; 2])], Error> {
//     let root = Path::new(env!("CARGO_MANIFEST_DIR"))
//         .join("pkg_data")
//         .join("feature_pairs");
//     let mut decoders = Vec::new_in(arena);
//     for notification in fs::read_dir(&root)? {
//         let dir_entry = notification?;
//         let pair_name = dir_entry
//             .file_name()
//             .into_string()
//             .map_err(|s| Error::msg(format!("unexpected folder name: {}", s.to_string_lossy())))?;
//         let filepath = root
//             .join(&pair_name)
//             .join(format!("{}.model.bin", model_train_base));
//         let bytes = read_file_into_global(&filepath)?;
//         let decoder_pair = <[GridDecoder; 2]>::deserialize(&mut &*bytes)?;
//         decoders.push((&*arena.alloc_str(&pair_name), decoder_pair));
//     }
//     Ok(leak_vec(decoders))
// }
// /// Replacing [`legacy::load_grid_scorer`].
// pub fn load_grid_scorer<'a>(
//     pair_name: &str,
//     z_grid_db_arena: &'a Bump,
// ) -> Result<GridScorer<'a>, Error> {
//     let root = Path::new(env!("CARGO_MANIFEST_DIR"))
//         .join("pkg_data")
//         .join("feature_pairs");
//     let filepath = root.join(pair_name).join("gridscorer.bin");
//     let bytes = read_file_into_global(&filepath)?;
//     GridScorer::deserialize(&mut &*bytes, z_grid_db_arena)
// }
// /// Replacing [`legacy::load_post_processor`].
// pub fn load_post_processor(
//     score_type: ScoreType,
//     model_train_base: ModelTrainingBase,
//     arena: &Bump,
// ) -> Result<PostProcessor<'_>, Error> {
//     match score_type {
//         ScoreType::Raw => Ok(PostProcessor::Raw),
//         ScoreType::ZScore => load_reference_scores(model_train_base, arena)
//             .map(|ref_scores| PostProcessor::new_zscore(ref_scores, arena)),
//         ScoreType::Percentile => {
//             load_reference_scores(model_train_base, arena).map(PostProcessor::new_percentile)
//         }
//     }
// }
// /// Replacing [`legacy::load_reference_scores`].
// pub fn load_reference_scores(model_train_base: ModelTrainingBase, arena: &Bump) -> Result<ReferenceFeatureMatrix<'_>, Error> {
//     let root = Path::new(env!("CARGO_MANIFEST_DIR"))
//         .join("pkg_data")
//         .join("human_reference_data");
//     let filepath = root.join(format!("{}.distr.bin", model_train_base));
//     let bytes = read_file_into_global(&filepath)?;
//     ReferenceFeatureMatrix::deserialize(&mut &*bytes, arena)
// }
