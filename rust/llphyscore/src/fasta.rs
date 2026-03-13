use crate::utils::alloc_dyn_writer;
use anyhow::Error;
use bumpalo::{Bump, boxed::Box, collections::Vec};
use llphyscore_core::{
    datatypes::{FastaEntry, aa_canonical_str},
    utils::{leak_vec, read_file},
};
use std::{
    fs::File,
    io::Write,
    path::{Path, PathBuf},
};

/// Helper function for callers of [`read_fasta`].
///
/// Construct a type erased thing that logs errors,
/// which can either be:
/// 1. `std::io::stderr`
/// 2. A `File` at a filepath given by the user.
pub fn alloc_error_handler(
    log_seq_errs_to: Option<PathBuf>,
    arena: &Bump,
) -> Result<Box<'_, dyn Write>, Error> {
    match log_seq_errs_to {
        Some(path) => {
            let writer = File::options().append(true).create(true).open(path)?;
            Ok(alloc_dyn_writer(writer, &arena))
        }
        None => {
            let writer = std::io::stderr().lock();
            Ok(alloc_dyn_writer(writer, &arena))
        }
    }
}

/// Read the file at `path` into a memory arena,
/// and parse it into a vector of [`FastaEntry`] structs.
///
/// Report errors to `err_out`, much like using `stderr`.
///
/// Errors
/// ------
/// This function fails if:
/// 1. The file at `path` cannot be read from.
/// 2. If the file at `path` is empty or does not start with a `>` character.
/// 3. If errors cannot be written to `err_out`.
pub fn read_fasta<'a>(
    path: &Path,
    arena: &'a Bump,
    err_out: &mut dyn Write,
) -> Result<&'a [FastaEntry<'a>], Error> {
    let bytes = read_file(path, arena)?;
    if bytes.is_empty() {
        return Err(Error::msg("got empty file"));
    }
    if !bytes.starts_with(&[b'>']) {
        return Err(Error::msg("got non-`>` character on first line"));
    }
    let mut entries = Vec::new_in(arena);
    for notification in parse_fasta_entries(leak_vec(bytes)) {
        match notification {
            Ok(entry) => {
                entries.push(entry);
            }
            Err(e) => {
                writeln!(err_out, "{}", e)?;
            }
        }
    }
    Ok(leak_vec(entries))
}
/// Iterate over a buffer with FASTA-formatted sequences
/// turned into [`FastaEntry`] structs, also returning parsing
/// errors as they come up.
///
/// This function will treat the first line of this slice as a header
/// (assumes it starts with `>` character). On debug builds this crashes
/// the program if the assumption is wrong.
fn parse_fasta_entries(
    mut slice: &mut [u8],
) -> impl Iterator<Item = Result<FastaEntry<'_>, Error>> {
    debug_assert!(slice.starts_with(b">"));
    std::iter::from_fn(move || {
        if slice.is_empty() {
            return None;
        }
        let working_slice: &mut [u8] = std::mem::replace(&mut slice, &mut []);
        let Some(idx) = end_of_current_header(working_slice) else {
            return Some(Err(format!(
                "could not parse entry (empty sequence): {}",
                String::from_utf8_lossy(working_slice.strip_prefix(b">").unwrap_or(working_slice))
            )));
        };
        let (header_slice, rest) = working_slice.split_at_mut(idx);
        debug_assert!(header_slice.starts_with(b">"));
        let header_slice = &header_slice[1..];
        let working_slice = &mut rest[1..];
        let sequence_slice: &mut [u8];
        match next_start_of_header(working_slice) {
            Some(idx) => {
                let rest: &mut [u8];
                ((sequence_slice, rest)) = working_slice.split_at_mut(idx);
                slice = &mut rest[1..];
            }
            None => {
                sequence_slice = working_slice;
                slice = &mut [];
            }
        }
        let header = match str::from_utf8(header_slice) {
            Ok(header) => header,
            Err(e) => {
                return Some(Err(format!(
                    "could not parse entry ({}): {}",
                    e,
                    String::from_utf8_lossy(header_slice)
                )));
            }
        };
        let r = match aa_canonical_str::join_multiline(sequence_slice) {
            Ok(sequence) => {
                if sequence.is_empty() {
                    Err(format!(
                        "could not parse entry (empty sequence): {}",
                        header
                    ))
                } else {
                    Ok(FastaEntry { header, sequence })
                }
            }
            Err(e) => Err(format!("could not parse entry ({}): {}", e, header)),
        };
        Some(r)
    })
    .map(|r| r.map_err(Error::msg))
}
/// Find the index of the next newline from the beginning of this slice.
fn end_of_current_header(bytes: &[u8]) -> Option<usize> {
    let b = bytes.iter().find(|b| **b == b'\n')?;
    let offset = (b as *const u8 as usize) - (bytes.as_ptr() as usize);
    Some(offset)
}
/// Find the index of the next newline followed by a `>` character
/// from the beginning of this slice.
fn next_start_of_header(bytes: &[u8]) -> Option<usize> {
    let ptr = bytes.windows(2).find(|pair| *pair == b"\n>")?;
    let offset = (ptr.as_ptr() as usize) - (bytes.as_ptr() as usize);
    Some(offset)
}
