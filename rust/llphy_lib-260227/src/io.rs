//! Module defining [`read_file`].
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use bytesize::ByteSize;
use std::{
    fs::File,
    io::{ErrorKind, Read, Seek},
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
    file.read_exact(&mut buf[..])?;
    ensure_file_at_eof(file, buf.len() as u64)?;
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
    let mut bytes = Vec::with_capacity_in(filesize, arena);
    // SAFETY: bytes have no initialization guarantees
    unsafe {
        bytes.set_len(filesize);
    };
    let buf = &mut *bytes;
    match file.read_exact(buf) {
        Ok(()) => (),
        // This probably doesn't happen but I don't know if
        // metadata is guaranteed to be accurate.
        Err(ref e) if matches!(e.kind(), ErrorKind::UnexpectedEof) => {
            let actual_size = file.seek(std::io::SeekFrom::End(0))?;
            return Err(wrong_file_size(buf.len() as _, actual_size));
        }
        Err(e) => return Err(e.into()),
    };
    ensure_file_at_eof(file, buf.len() as _)?;
    Ok(bytes)
}
/// A sanity check for me to make sure we read the whole file.
fn ensure_file_at_eof(mut file: File, current_position: u64) -> Result<(), Error> {
    let end_position = file.seek(std::io::SeekFrom::End(0))?;
    if end_position == current_position {
        Ok(())
    } else {
        // This probably doesn't happen but I don't know if
        // metadata is guaranteed to be accurate.
        Err(wrong_file_size(current_position, end_position))
    }
}
/// The error returned when the number of bytes read from
/// the file doesn't match the metadata size.
#[cold]
fn wrong_file_size(expected_size: u64, actual_size: u64) -> Error {
    Error::msg(format!(
        "file was expected to be {} because of metadata but it was actually {}",
        ByteSize::b(expected_size),
        ByteSize::b(actual_size)
    ))
}
