use anyhow::Error;
use std::{
    collections::HashMap, fs::File, io::{Cursor, Read}, path::{Path, PathBuf}
};
mod pdb_statistics;
pub use pdb_statistics::load_legacy_pdb_statistics;
fn read_archive_file(filepath: &Path) -> Result<Vec<u8>, Error> {
    const DEV_ARCHIVE_ROOT: &'static str = env!("LEGACY_ARCHIVE_DIR");
    let filepath = Path::new(DEV_ARCHIVE_ROOT).join(filepath);
    let mut file = File::open(filepath)?;
    let filesize = file.metadata()?.len() as usize;
    let mut bytes = Vec::with_capacity(filesize);
    unsafe {
        bytes.set_len(filesize);
    }
    file.read_exact(&mut bytes[..])?;
    Ok(bytes)
}