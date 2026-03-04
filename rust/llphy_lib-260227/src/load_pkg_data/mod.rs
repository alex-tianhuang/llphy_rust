use crate::io::read_file_into_global;
use anyhow::Error;
use std::path::{Path, PathBuf};
mod grid_decoder;
mod grid_scorer;
mod post_processor;
mod reference_scores;
pub use grid_decoder::load_grid_decoders;
pub use grid_scorer::load_grid_scorer;
pub use post_processor::load_post_processor;
use reference_scores::load_reference_scores;
/// Path relative to the cargo manifest directory
/// of the legacy data directory.
const DATA_DIR: &'static str = "pkg_data/legacy";

/// Load a file from the archive using a relative filepath.
///
/// Uses [`read_file_into_global`].
fn read_archive_file(filepath: &Path) -> Result<Vec<u8>, Error> {
    let mut root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    for name in DATA_DIR.split("/") {
        root = root.join(name);
    }
    read_file_into_global(&root.join(filepath))
}
