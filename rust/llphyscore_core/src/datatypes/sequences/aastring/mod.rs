//! Module defining aminoacid string [`aa_canonical_str`].
//! 
//! Dev note
//! --------
//! It seems a bit more complicated than it needs to be
//! because this submodule comes from a crate where the
//! definition of "valid aminoacid character" can vary
//! depending on whether gaps or non-canonical aminoacids
//! are allowed, the most general implementation are in
//! the [`generic`] module.
pub(crate) mod generic;
use crate::datatypes::Aminoacid;
use generic::AALike;

/// Borrowed `str` analogue for [`Aminoacid`] strings.
#[allow(non_camel_case_types)]
pub type aa_canonical_str = generic::aa_str<Aminoacid>;
unsafe impl AALike for Aminoacid {
    const DESCRIBE: &'static str = "single-letter aminoacid character";
}
