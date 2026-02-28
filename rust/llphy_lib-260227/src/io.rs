//! Module defining [`read_fasta`]
use crate::{datatypes::FastaEntry, fasta::parse_fasta_entries};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use bytesize::ByteSize;
use std::{
    fs::File,
    io::{ErrorKind, Read, Seek},
    path::{Path, PathBuf},
    ptr::slice_from_raw_parts_mut,
};

/// Read the file at `path` into a memory arena, and parse it into
/// a vector of [`FastaEntry`] structs.
/// 
/// In addition to failing if the file cannot be read from, it will
/// fail if the file is empty or does not start with a `>` character.
/// 
/// IO
/// --
/// May log warnings to stderr, because it calls [`parse_fasta_entries`]
/// which logs to stderr.
pub fn read_fasta(path: PathBuf, arena: &Bump) -> Result<Vec<'_, FastaEntry<'_>>, Error> {
    let bytes = read_file(&path, arena)?;
    if bytes.is_empty() {
        return Err(Error::msg(format!(
            "expected non-empty fasta file at {}, got empty file",
            path.display()
        )));
    }
    if !bytes.starts_with(&[b'>']) {
        return Err(Error::msg(format!(
            "expected fasta file at {}, got non-`>` character on first line",
            path.display()
        )));
    }
    Ok(parse_fasta_entries(bytes, arena))
}
/// Read the file at `path` into the given memory arena
/// as a vector of bytes.
pub fn read_file<'a, 'b>(path: &'a Path, arena: &'b Bump) -> Result<Vec<'b, u8>, Error> {
    let mut file = File::open(path)?;
    // Truncation made by casting does not matter.
    // If the filesize is greater than `usize::MAX`
    // you will have different problems.
    let filesize = file.metadata()?.len() as usize;
    let mut bytes = Vec::with_capacity_in(filesize, arena);
    // SAFETY: bytes have no initialization guarantees
    let buf = unsafe { &mut *slice_from_raw_parts_mut(bytes.as_mut_ptr(), filesize) };
    match file.read_exact(buf) {
        Ok(()) => (),
        // This probably doesn't happen but I don't know if
        // metadata is guaranteed to be accurate.
        Err(ref e) if matches!(e.kind(), ErrorKind::UnexpectedEof) => {
            let actual_size = file.seek(std::io::SeekFrom::End(0))?;
            return Err(wrong_file_size(buf.len() as _, actual_size))
        },
        Err(e) => return Err(e.into())
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