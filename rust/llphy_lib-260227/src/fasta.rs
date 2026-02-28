use crate::datatypes::{Aminoacid, FastaEntry, aa_canonical_str};
use bumpalo::{Bump, collections::Vec};
use std::mem::{self, ManuallyDrop};

pub fn parse_fasta_entries<'a>(bytes: Vec<'a, u8>, arena: &'a Bump) -> Vec<'a, FastaEntry<'a>> {
    let mut slice = leak_vec(bytes);
    let mut entries = Vec::new_in(arena);
    while !slice.is_empty() {
        let Some(idx) = end_of_current_header(slice) else {
            eprintln!(
                "could not parse entry (empty sequence): {}",
                String::from_utf8_lossy(slice)
            );
            return entries;
        };
        let (header_slice, rest) = slice.split_at_mut(idx);
        debug_assert!(!header_slice.is_empty());
        debug_assert_eq!(header_slice[0], b'>');
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
        match parse_aastr(sequence_slice) {
            ParseResult::Ok(s) => {
                sequence = s;
            }
            ParseResult::Warn(s) => {
                eprintln!(
                    "non aminoacid bytes removed while parsing this entry: {}",
                    header
                );
                sequence = s;
            }
            ParseResult::Empty => {
                eprintln!("could not parse entry (empty sequence): {}", header);
                continue;
            }
        }
        entries.push(FastaEntry { header, sequence })
    }
    entries
}
fn end_of_current_header(bytes: &[u8]) -> Option<usize> {
    let b = bytes.iter().find(|b| **b == b'\n')?;
    let offset = (b as *const u8 as usize) - (bytes.as_ptr() as usize);
    Some(offset)
}
fn next_start_of_header(bytes: &[u8]) -> Option<usize> {
    for pair in bytes.windows(2) {
        if pair == &[b'\n', b'>'] {
            let offset = (pair.as_ptr() as usize) - (bytes.as_ptr() as usize);
            return Some(offset);
        }
    }
    None
}
pub enum ParseResult<'a> {
    Ok(&'a aa_canonical_str),
    Warn(&'a aa_canonical_str),
    Empty,
}
fn parse_aastr(slice: &mut [u8]) -> ParseResult<'_> {
    let mut non_aa_start_idx = None;
    let mut warn = false;
    if slice.is_empty() {
        return ParseResult::Empty;
    }
    for i in 0..slice.len() {
        let b = slice[i];
        if let Ok(_) = Aminoacid::try_from(b) {
            continue;
        }
        if let Ok(_) = Aminoacid::try_from(b.to_ascii_uppercase()) {
            slice[i] = b.to_ascii_uppercase();
            continue;
        };
        non_aa_start_idx = Some(i + 1);
        warn = should_warn(b);
        break;
    }
    let Some(start) = non_aa_start_idx else {
        // SAFETY: to reach this point all bytes in `slice`
        //         having never reached the break statement
        //         in the loop above
        unsafe { return ParseResult::Ok(aa_canonical_str::from_bytes_unchecked(slice)) }
    };
    let mut shift = 1;
    for i in start..slice.len() {
        let b = slice[i];
        if let Ok(_) = Aminoacid::try_from(b.to_ascii_uppercase()) {
            slice[i - shift] = b.to_ascii_uppercase();
            continue;
        };
        shift += 1;
        warn |= should_warn(b);
    }
    if slice.len() == shift {
        return ParseResult::Empty;
    }

    let len = slice.len() - shift;
    // SAFETY: safe for the same reason `Vec::retain` is safe.
    let s = unsafe { aa_canonical_str::from_bytes_unchecked(&slice[..len]) };
    if warn {
        ParseResult::Warn(s)
    } else {
        ParseResult::Ok(s)
    }
}
fn should_warn(b: u8) -> bool {
    !b.is_ascii_whitespace()
}
/// Leak a [`Vec`] managed by an arena into a slice
/// that lives for as long as the arena does not reset.
pub fn leak_vec<'a, T>(buf: Vec<'a, T>) -> &'a mut [T] {
    let mut buf = ManuallyDrop::new(buf);
    // SAFETY: the arena must de-allocate before this vec does.
    //         Also, do to the use of `ManuallyDrop`, these bytes
    //         are not de-allocated and reused by the arena.
    unsafe { mem::transmute::<&mut [T], &'a mut [T]>(buf.as_mut_slice()) }
}
