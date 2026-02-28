mod sequences;
pub(crate) use sequences::{
    AACanonicalString, AACanonicalStringStrict, AAIndex, AAMap, AASet, AMINOACIDS, Aminoacid,
    FastaEntry, NotAminoacidError, NotCanonicalAAStrError, aa_canonical_str,
};
