//! Functions for loading models and reference data.
//!
//! Database file structure
//! -----------------------
//! 
//! This module assumes a specific directory structure
//! for models and reference data, like so:
//! ```text
//! |- `feature_pairs`
//! |  |- `S2.SUMPI`
//! |  |  |- `gridscorer.bin`
//! |  |  |- `human.model.bin`
//! |  |  |- `human+PDB.model.bin`
//! |  |  \- `PDB.model.bin`
//! |  \- `S3.WATER.V2`
//! |  |  |- `gridscorer.bin`
//! |  ... (for each feature pair)
//! |
//! \- `human_reference_data`
//!    |- `human.distr.bin`
//!    |- `human+PDB.distr.bin`
//!    \- `PDB.distr.bin`
//! ```
//! As you can see above, the three types of files are:
//! 1. `gridscorer.bin`, which contains a [`GridScorer`]
//! 2. `...model.bin`, which contains [`GridDecoderPair`]s
//!    that have been trained on human, human+PDB, or PDB negative sets
//!    (i.e. the variants of [`ModelTrainingBase`])
//! 3. `...distr.bin`, which contains a [`ReferenceFeatureMatrix`]
//!    which represents feature scores of human IDRs that have been
//!    obtained by using `human`, `human+PDB`, or `PDB` [`GridDecoder`]s.
use crate::{
    datatypes::{
        GridDecoderPair, GridScorer, ModelTrainingBase, PostProcessor, ReferenceFeatureMatrix, ScoreType
    },
    utils::leak_vec,
};
use anyhow::Error;
use borsh::BorshDeserialize;
use bumpalo::{Bump, collections::Vec};
use std::{
    fs::{self, File},
    io::Read,
    path::Path,
};
/// Load all [`GridDecoderPair`]s in a given database.
/// 
/// These are the `...model.bin` files described in the
/// [module level docs](self).
pub fn load_grid_decoders<'a>(
    pkg_data_root: &Path,
    model_train_base: ModelTrainingBase,
    arena: &'a Bump,
) -> Result<&'a [(&'a str, GridDecoderPair)], Error> {
    let feature_pairs_root = pkg_data_root.join("feature_pairs");
    let mut decoders = Vec::with_capacity_in(fs::read_dir(&feature_pairs_root)?.count(), arena);
    let mut byte_buffer = std::vec::Vec::new();
    for notification in fs::read_dir(&feature_pairs_root)? {
        let dir_entry = notification?;
        let pair_name = dir_entry
            .file_name()
            .into_string()
            .map_err(|s| Error::msg(format!("unexpected folder name: {}", s.to_string_lossy())))?;
        let filepath = feature_pairs_root
            .join(&pair_name)
            .join(format!("{}.model.bin", model_train_base));
        byte_buffer.clear();
        File::open(&filepath)?.read_to_end(&mut byte_buffer)?;
        let decoder_pair = GridDecoderPair::deserialize(&mut &*byte_buffer)?;
        decoders.push((&*arena.alloc_str(&pair_name), decoder_pair));
    }
    Ok(leak_vec(decoders))
}
/// Load a single [`GridScorer`]s in a given database
/// into the given memory arena (by feature pair name).
/// 
/// These are the `gridscorer.bin` files described in the
/// [module level docs](self).
pub fn load_grid_scorer<'a>(pkg_data_root: &Path, pair_name: &str, arena: &'a Bump) -> Result<&'a GridScorer<'a>, Error> {
    let filepath = pkg_data_root
        .join("feature_pairs").join(pair_name).join("gridscorer.bin");
    let mut byte_buffer = std::vec::Vec::new();
    File::open(&filepath)?.read_to_end(&mut byte_buffer)?;
    GridScorer::deserialize(&mut &*byte_buffer, arena)
}
/// If requested (via `score_type`), load reference data
/// and then convert it to the appropriate post-processor.
/// 
/// The data is loaded from the `...distr.bin` files
/// described in the [module level docs](self).
pub fn load_post_processor<'a>(
    pkg_data_root: &Path,
    score_type: ScoreType,
    model_train_base: ModelTrainingBase,
    arena: &'a Bump,
) -> Result<PostProcessor<'a>, Error> {
    match score_type {
        ScoreType::Raw => Ok(PostProcessor::Raw),
        ScoreType::ZScore => load_reference_scores(pkg_data_root, model_train_base, arena)
            .map(|ref_scores| PostProcessor::new_zscore(ref_scores, arena)),
        ScoreType::Percentile => {
            load_reference_scores(pkg_data_root, model_train_base, arena).map(PostProcessor::new_percentile)
        }
    }
}
/// Load a single [`ReferenceFeatureMatrix`]s in a given database
/// into the given memory arena (by [`ModelTrainingBase`]).
/// 
/// These are the `...distr.bin` files described in the
/// [module level docs](self).
pub fn load_reference_scores<'a>(
    pkg_data_root: &Path,
    model_train_base: ModelTrainingBase,
    arena: &'a Bump,
) -> Result<ReferenceFeatureMatrix<'a>, Error> {
    let filepath = pkg_data_root
        .join("human_reference_data").join(format!("{}.distr.bin", model_train_base));
    let mut byte_buffer = std::vec::Vec::new();
    File::open(&filepath)?.read_to_end(&mut byte_buffer)?;
    ReferenceFeatureMatrix::deserialize(&mut &*byte_buffer, arena)
}
