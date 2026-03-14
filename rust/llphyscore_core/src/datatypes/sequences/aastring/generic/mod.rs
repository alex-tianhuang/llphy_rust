//! Module defining a generic "aminoacid string" datatype ([`aa_str`]),
//! which is just a wrapper around slices that can be cast safely to strings.
//!
//! The reason this is so complicated just to represent aminoacid strings
//! is because this module comes from a larger codebase where there are
//! actually multiple aminoacid variants and so it is useful to have a
//! generic `aa_str`.
use std::fmt::Debug;
mod borrowed;
pub use borrowed::aa_str;
/// A one-byte datatype that is representable by a subset
/// (single-byte, not arbitrary chars) of displayable characters.
///
/// Safety
/// ------
/// Only implement this trait if you are sure a slice of these bytes
/// can be safely cast as a string slice of single byte characters.
pub unsafe trait AALike: Sized + Debug + Copy + TryFrom<char> + Into<u8> {
    /// A lowercase string describing what bytes of this type represent.
    ///
    /// Make sure it reads well in the sentence:
    /// `expected string of {DESCRIBE}s`
    const DESCRIBE: &'static str;
}
