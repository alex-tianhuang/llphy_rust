//! Some utility functions for the crate
//! that I couldn't bother organizing.

/// Leak a [`Vec`] managed by an arena into a slice
/// that lives for as long as the arena does not reset.
pub(crate) fn leak_vec<'a, T>(buf: bumpalo::collections::Vec<'a, T>) -> &'a mut [T] {
    let mut buf = std::mem::ManuallyDrop::new(buf);
    // SAFETY: the arena must de-allocate before this vec does.
    //         Also, do to the use of `ManuallyDrop`, these bytes
    //         are not de-allocated and reused by the arena.
    unsafe { std::mem::transmute::<&mut [T], &'a mut [T]>(buf.as_mut_slice()) }
}

/// Boilerplate that looks like `#[serde(from)]`.
macro_rules! derive_borsh_de_from {
    ($my_ty:ty as $std_ty:ty, $from_impl:expr) => {
        impl borsh::BorshDeserialize for $my_ty {
            fn deserialize(buf: &mut &[u8]) -> std::io::Result<Self> {
                <$std_ty as borsh::BorshDeserialize>::deserialize(buf).map($from_impl)
            }
            fn deserialize_reader<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
                <$std_ty as borsh::BorshDeserialize>::deserialize_reader(reader).map($from_impl)
            }
        }
    };
}

/// Boilerplate that looks like `#[serde(into)]`.
macro_rules! derive_borsh_se_into {
    ($my_ty:ty as $std_ty:ty, |$this:ident| $to_impl:expr) => {
        impl borsh::BorshSerialize for $my_ty {
            fn serialize<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
                let $this = self;
                let intermediate: $std_ty = $to_impl;
                intermediate.serialize(writer)
            }
        }
    };
}
use {
    bumpalo::Bump,
    std::{iter::Cloned, mem::ManuallyDrop},
};
pub(crate) use {derive_borsh_de_from, derive_borsh_se_into};
