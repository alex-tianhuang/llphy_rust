use anyhow::Error;
use std::path::Path;
mod model;
mod pdb_statistics;
mod reference_g2w_scores;
pub use model::load_legacy_model_and_feature_order;
pub use pdb_statistics::load_legacy_pdb_statistics;
pub use reference_g2w_scores::load_legacy_reference_g2w_scores;

use crate::io::read_file_into_global;
fn read_archive_file(filepath: &Path) -> Result<Vec<u8>, Error> {
    const DEV_ARCHIVE_ROOT: &'static str = env!("LEGACY_ARCHIVE_DIR");
    let filepath = Path::new(DEV_ARCHIVE_ROOT).join(filepath);
    read_file_into_global(&filepath)
}
pub const LEGACY_FEATURE_NAMES: &'static [&'static str] = &[
    "protein-water",
    "protein-carbon",
    "hydrogen bond (long-range)",
    "pi-pi (long-range)",
    "disorder (long)",
    "K-Beta similarity",
    "disorder (short)",
    "electrostatic (short-range)",
];
