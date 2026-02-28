/// A fork of my other project's aminoacid/sequence related definitions,
/// except for [`FastaEntry`] which I made for this project.
mod aamap;
mod aaset;
mod aastring;
mod aminoacid;
pub(crate) use aamap::const_aamap;
pub use aamap::{AAMap, 
    // AAWeights
};
pub use aaset::AASet;
pub(crate) use aastring::{
    AACanonicalString, AACanonicalStringStrict, NotCanonicalAAStrError, aa_canonical_str,
    generic,
};
pub(crate) use aminoacid::{AAIndex, AMINOACIDS, Aminoacid, NotAminoacidError};

/// Fasta entry, containing a sequence and header.
#[derive(Clone, Copy)]
pub struct FastaEntry<'a> {
    pub header: &'a str,
    pub sequence: &'a aa_canonical_str
}