/// A fork of my other project's aminoacid/sequence related definitions,
/// except for [`FastaEntry`] which I made for this project.
mod aamap;
mod aastring;
mod aminoacid;
pub use aamap::AAMap;
pub use aastring::aa_canonical_str;
pub use aminoacid::{AAIndex, AMINOACIDS, Aminoacid};

/// Fasta entry, containing a sequence and header.
pub struct FastaEntry<'a> {
    pub header: &'a str,
    pub sequence: &'a aa_canonical_str,
}
