mod legacy;
use std::{fs, path::PathBuf};
use anyhow::Error;
use borsh::BorshDeserialize;
use bumpalo::{Bump, collections::Vec};
pub use legacy::{load_grid_scorer, load_post_processor};
use pyo3::Python;
use crate::{datatypes::ModelTrainingBase, featurizer::GridDecoder, io::read_file_into_global, leak_vec};
/// Path relative to the cargo manifest directory
/// of the legacy data directory.
const DATA_DIR: &'static str = "pkg_data";

/// Replacing [`legacy::load_grid_decoders`].
pub fn load_grid_decoders<'a>(
    model_train_base: ModelTrainingBase,
    arena: &'a Bump,
    py: Python,
) -> Result<&'a [(&'a str, [GridDecoder; 2])], Error> {
    let mut root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    for name in DATA_DIR.split("/") {
        root = root.join(name);
    }
    root = root.join("feature_pairs");
    let mut decoders = Vec::new_in(arena);
    for notification in fs::read_dir(&root)? {
        let dir_entry = notification?;
        let pair_name = dir_entry.file_name().into_string().map_err(|s| Error::msg(format!("unexpected folder name: {}", s.to_string_lossy())))?;
        let filepath = root.join(&pair_name).join(format!("{}.model.bin", model_train_base));
        let bytes = read_file_into_global(&filepath)?;
        let decoder_pair = <[GridDecoder; 2]>::deserialize(&mut &*bytes)?;
        decoders.push((&*arena.alloc_str(&pair_name), decoder_pair));
    };
    Ok(leak_vec(decoders))
}