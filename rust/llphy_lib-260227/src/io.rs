//! Module defining various IO utilities
//! like [`read_file`] and [`read_file_into_global`].
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use std::{
    fs::File,
    io::Read,
    path::Path,
};

/// Read the file at `path` into a byte vector
/// managed by the global allocator.
/// 
/// For non-persistent (and large!) files.
pub fn read_file_into_global(path: &Path) -> Result<std::vec::Vec<u8>, Error>{
    let mut file = File::open(path)?;
    let filesize = file.metadata()?.len() as usize;
    let mut buf = std::vec::Vec::with_capacity(filesize);
    unsafe {
        buf.set_len(filesize);
    }
    file.read_exact(&mut buf)?;
    Ok(buf)
}
/// Read the file at `path` into the given memory arena
/// as a vector of bytes.
pub fn read_file<'a>(path: &Path, arena: &'a Bump) -> Result<Vec<'a, u8>, Error> {
    let mut file = File::open(path)?;
    // Truncation made by casting does not matter.
    // If the filesize is greater than `usize::MAX`
    // you will have different problems.
    let filesize = file.metadata()?.len() as usize;
    let mut buf = Vec::with_capacity_in(filesize, arena);
    // SAFETY: bytes have no initialization guarantees
    unsafe {
        buf.set_len(filesize);
    };
    
    file.read_exact(&mut buf)?;
    Ok(buf)
}