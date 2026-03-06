mod legacy;
use crate::{
    datatypes::ModelTrainingBase,
    featurizer::{GridDecoder, GridScorer},
    io::read_file_into_global,
    leak_vec,
};
use anyhow::Error;
use borsh::BorshDeserialize;
use bumpalo::{Bump, collections::Vec};
pub use legacy::load_post_processor;
use pyo3::Python;
use std::{fs, path::Path};

/// Replacing [`legacy::load_grid_decoders`].
pub fn load_grid_decoders<'a>(
    model_train_base: ModelTrainingBase,
    arena: &'a Bump,
    py: Python,
) -> Result<&'a [(&'a str, [GridDecoder; 2])], Error> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("pkg_data")
        .join("feature_pairs");
    let mut decoders = Vec::new_in(arena);
    for notification in fs::read_dir(&root)? {
        let dir_entry = notification?;
        let pair_name = dir_entry
            .file_name()
            .into_string()
            .map_err(|s| Error::msg(format!("unexpected folder name: {}", s.to_string_lossy())))?;
        let filepath = root
            .join(&pair_name)
            .join(format!("{}.model.bin", model_train_base));
        let bytes = read_file_into_global(&filepath)?;
        let decoder_pair = <[GridDecoder; 2]>::deserialize(&mut &*bytes)?;
        decoders.push((&*arena.alloc_str(&pair_name), decoder_pair));
    }
    Ok(leak_vec(decoders))
}
/// Replacing [`legacy::load_grid_scorer`].
pub fn load_grid_scorer<'a>(
    pair_name: &str,
    z_grid_db_arena: &'a Bump,
) -> Result<GridScorer<'a>, Error> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("pkg_data")
        .join("feature_pairs");
    let filepath = root.join(pair_name).join("gridscorer.bin");
    let bytes = read_file_into_global(&filepath)?;
    GridScorer::deserialize(&mut &*bytes, z_grid_db_arena)
}
