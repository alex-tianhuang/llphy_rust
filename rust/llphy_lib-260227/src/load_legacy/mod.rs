use anyhow::Error;
use std::path::Path;
mod pdb_statistics;
mod reference_g2w_scores;
pub use pdb_statistics::load_legacy_pdb_statistics;
pub use reference_g2w_scores::load_legacy_reference_g2w_scores;

use crate::io::read_file_into_global;
fn read_archive_file(filepath: &Path) -> Result<Vec<u8>, Error> {
    const DEV_ARCHIVE_ROOT: &'static str = env!("LEGACY_ARCHIVE_DIR");
    let filepath = Path::new(DEV_ARCHIVE_ROOT).join(filepath);
    read_file_into_global(&filepath)
}
