//! Crate defining functionality used by most
//! of the other crates in this workspace.
#![cfg_attr(feature = "simd", feature(portable_simd))]
pub mod datatypes;
pub mod featurizer;
pub mod utils;
pub mod load_pkg_data;