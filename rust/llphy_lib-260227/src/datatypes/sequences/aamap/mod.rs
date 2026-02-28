//! Module defining the [`AAMap`] type.

use crate::datatypes::{AMINOACIDS, Aminoacid};
use std::{
    fmt::Debug,
    ops::{Index, IndexMut},
};
// mod aaweights;
// pub use aaweights::AAWeights;

/// A residue mapping type.
///
/// Uses [`Aminoacid::to_aaindex`] as an indexing function.
/// Eagerly allocates 20 slots. May be better to use
/// `Option<Box<T>>` than `Option<T>` for large,
/// lazily constructed types.
///
/// Usage
/// -----
/// 
/// Construct like so:
/// ```
/// let map = AAMap([0.0_f32; 20]);
/// ```
/// or like so:
/// ```
/// let map = <AAMap<Option<Box<Foo>>>>::default();
/// ```
///
/// And index into it like so:
/// ```
/// let value = map[Aminoacid::R];
/// map[Aminoacid::T] = value;
/// ```
///
/// Iterate over values like so:
/// ```
/// for (aa, value) in map.iter() {
///     /* ... */
/// }
/// ```
#[derive(Clone, Default, PartialEq)]
pub struct AAMap<T>(pub [T; AMINOACIDS.len()]);
impl<T> AAMap<T> {
    /// Iterate over (aa, &value) pairs, like a regular map.
    pub fn iter(&self) -> impl Iterator<Item = (Aminoacid, &T)> {
        AMINOACIDS.into_iter().zip(self.0.iter())
    }
    /// Iterate over (aa, &mut value) pairs, like a regular map.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (Aminoacid, &mut T)> {
        AMINOACIDS.into_iter().zip(self.0.iter_mut())
    }
    /// Iterate over values by reference, like a regular map.
    /// The order of the iterator is still the same as the declaration
    /// order of amino acids.
    pub fn values(&self) -> impl Iterator<Item = &T> {
        self.0.iter()
    }
    /// Iterate over values by mutable reference, like a regular map.
    /// The order of the iterator is still the same as the declaration
    /// order of amino acids.
    pub fn values_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.0.iter_mut()
    }
    /// Map this residue map to a new one by applying `f` to value references.
    pub fn map_refs<V>(&self, f: impl FnMut(&T) -> V) -> AAMap<V> {
        AAMap(self.0.each_ref().map(f))
    }
}
/// Construct an [`AAMap`] from an array of ([`Aminoacid`], value)
/// pairs as a `const` expression.
macro_rules! const_aamap {
    ([(_, $default:expr),$($pair:expr),+]) => {{
        let mut map = [$default; 20];
        $(let (aa, value): (crate::datatypes::Aminoacid, _) = $pair; map[aa.to_aaindex() as usize] = value;)+
        crate::datatypes::AAMap(map)
    }};
}
pub(crate) use const_aamap;
impl<T> Index<Aminoacid> for AAMap<T> {
    type Output = T;
    fn index(&self, index: Aminoacid) -> &Self::Output {
        self.0.index(index.to_aaindex() as usize)
    }
}
impl<T> IndexMut<Aminoacid> for AAMap<T> {
    fn index_mut(&mut self, index: Aminoacid) -> &mut Self::Output {
        self.0.index_mut(index.to_aaindex() as usize)
    }
}
impl<T: Default> FromIterator<(Aminoacid, T)> for AAMap<T> {
    fn from_iter<I: IntoIterator<Item = (Aminoacid, T)>>(iter: I) -> Self {
        let mut return_value = Self::default();
        for (aa, value) in iter {
            return_value[aa] = value
        }
        return return_value;
    }
}
impl<T> IntoIterator for AAMap<T> {
    type IntoIter =
        std::iter::Zip<std::array::IntoIter<Aminoacid, 20>, std::array::IntoIter<T, 20>>;
    type Item = (Aminoacid, T);
    fn into_iter(self) -> Self::IntoIter {
        AMINOACIDS.into_iter().zip(self.0.into_iter())
    }
}
impl<T: Debug> Debug for AAMap<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_map().entries(self.iter()).finish()
    }
}
