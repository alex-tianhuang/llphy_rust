//! Module defining [`read_fasta`].
use crate::{
    datatypes::{Aminoacid, FastaEntry, aa_canonical_str},
    io::read_file,
    leak_vec,
};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use std::path::PathBuf;

/// Read the file at `path` into a memory arena, and parse it into
/// a vector of [`FastaEntry`] structs.
///
/// In addition to failing if the file cannot be read from, it will
/// fail if the file is empty or does not start with a `>` character.
///
/// IO
/// --
/// May log warnings to stderr,
/// because it calls [`parse_fasta_entries`] which logs to stderr.
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
/// Parse a buffer with FASTA-formatted sequences
/// into a vector of [`FastaEntry`] structs.
///
/// IO
/// --
/// This function logs to stderr when it encounters:
/// - Headers that are non-UTF8
/// - Sequences that contain non aminoacid (lowercase or uppercase) characters
/// - Sequences that are empty (after removing all non aminoacid characters)
fn parse_fasta_entries<'a>(bytes: Vec<'a, u8>, arena: &'a Bump) -> Vec<'a, FastaEntry<'a>> {
    let mut slice = leak_vec(bytes);
    let mut entries = Vec::new_in(arena);
    while !slice.is_empty() {
        let Some(idx) = end_of_current_header(slice) else {
            eprintln!(
                "could not parse entry (empty sequence): {}",
                String::from_utf8_lossy(slice.strip_prefix(&[b'>']).unwrap_or(slice))
            );
            return entries;
        };
        let (header_slice, rest) = slice.split_at_mut(idx);
        debug_assert!(!header_slice.is_empty());
        debug_assert_eq!(header_slice[0], b'>');
        let header_slice = &header_slice[1..];
        slice = &mut rest[1..];
        let rest: &mut [u8];
        let sequence_slice: &mut [u8];
        match next_start_of_header(slice) {
            Some(idx) => {
                ((sequence_slice, rest)) = slice.split_at_mut(idx);
                slice = &mut rest[1..];
            }
            None => {
                ((sequence_slice, rest)) = slice.split_at_mut(slice.len());
                slice = rest;
            }
        }
        let header: &str;
        match str::from_utf8(header_slice) {
            Ok(h) => {
                header = h;
            }
            Err(e) => {
                eprintln!(
                    "could not parse entry ({}): {}",
                    e,
                    String::from_utf8_lossy(header_slice)
                );
                continue;
            }
        }
        let sequence: &aa_canonical_str;
        match parse_sequence(sequence_slice) {
            ParseResult::Ok(s) => {
                sequence = s;
            }
            ParseResult::WarnInvalidBytes(s) => {
                eprintln!(
                    "sequence was modified, non aminoacid bytes removed while parsing this entry: {}\nmodified sequence: {}",
                    header,
                    s.as_str(),
                );
                sequence = s;
            }
            ParseResult::Empty => {
                eprintln!("could not parse entry (empty sequence): {}", header);
                continue;
            }
            ParseResult::EmptyAndInvalidBytes => {
                eprintln!("could not parse entry (contained invalid bytes): {}", header);
                continue;
            }
        }
        entries.push(FastaEntry { header, sequence })
    }
    entries
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
    for pair in bytes.windows(2) {
        if pair == &[b'\n', b'>'] {
            let offset = (pair.as_ptr() as usize) - (bytes.as_ptr() as usize);
            return Some(offset);
        }
    }
    None
}
/// All possible cases returned by [`parse_sequence`].
enum ParseResult<'a> {
    /// The sequence is non-empty and does not merit any warning.
    Ok(&'a aa_canonical_str),
    /// The sequence contained non-aminoacid and non-whitespace bytes,
    /// but after removing those it yielded a non-empty sequence.
    WarnInvalidBytes(&'a aa_canonical_str),
    /// The sequence contained no characters or only whitespace.
    Empty,
    /// The sequence contained characters but none of them were
    /// convertible to aminoacids.
    EmptyAndInvalidBytes,
}
/// Parse an [`aa_canonical_str`] from the given slice of bytes.
///
/// Parsing is done loosely, converting lowercase to uppercase and
/// removing unexpected characters rather than failing outright.
/// It will only fail if the sequence does not contain any aminoacids.
fn parse_sequence(slice: &mut [u8]) -> ParseResult<'_> {
    let mut non_aa_start_idx = None;
    let mut warn = false;
    if slice.is_empty() {
        return ParseResult::Empty;
    }
    for i in 0..slice.len() {
        let b = slice[i];
        if Aminoacid::try_from(b).is_ok() {
            continue;
        }
        if Aminoacid::try_from(b.to_ascii_uppercase()).is_ok() {
            slice[i] = b.to_ascii_uppercase();
            continue;
        };
        non_aa_start_idx = Some(i);
        warn = should_warn(b);
        break;
    }
    let Some(start) = non_aa_start_idx else {
        return ParseResult::Ok(aa_canonical_str::from_bytes(slice).unwrap())
    };
    let mut len = start;
    let mut cursor = start;
    loop {
        cursor += 1;
        let Some(&b) = slice.get(cursor) else { break };
        if Aminoacid::try_from(b.to_ascii_uppercase()).is_err() {
            warn |= should_warn(b);
            continue;
        };
        slice[len] = b.to_ascii_uppercase();
        len += 1;
    }
    if len == 0 {
        if warn {
            return ParseResult::EmptyAndInvalidBytes;
        } else {
            return ParseResult::Empty;
        }
    }
    let s = aa_canonical_str::from_bytes(&slice[..len]).unwrap();
    if warn {
        ParseResult::WarnInvalidBytes(s)
    } else {
        ParseResult::Ok(s)
    }
}
/// Checks whether this byte merits a warning for this sequence.
///
/// Right now I ignore whitespace characters.
fn should_warn(b: u8) -> bool {
    !b.is_ascii_whitespace()
}
