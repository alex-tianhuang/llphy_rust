use crate::{datatypes::FastaEntry, fasta::parse_fasta_entries};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use bytesize::ByteSize;
use std::{
    fs::File,
    io::{Read, Seek},
    path::{Path, PathBuf},
    ptr::slice_from_raw_parts_mut,
};

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
pub fn read_file<'a, 'b>(path: &'a Path, arena: &'b Bump) -> Result<Vec<'b, u8>, Error> {
    let mut file = File::open(path)?;
    // Truncation made by casting does not matter.
    // If the filesize is greater than `usize::MAX`
    // you will have different problems.
    let filesize = file.metadata()?.len() as usize;
    let mut bytes = Vec::with_capacity_in(filesize, arena);
    let buf = unsafe { &mut *slice_from_raw_parts_mut(bytes.as_mut_ptr(), filesize) };
    file.read_exact(buf)?;
    ensure_file_at_eof(file, buf.len() as _)?;
    Ok(bytes)
}
fn ensure_file_at_eof(mut file: File, current_position: u64) -> Result<(), Error> {
    let end_position = file.seek(std::io::SeekFrom::End(0))?;
    if end_position == current_position {
        Ok(())
    } else {
        Err(Error::msg(format!(
            "file was expected to be {} because of metadata but it was actually {}",
            ByteSize::b(current_position),
            ByteSize::b(end_position)
        )))
    }
}
