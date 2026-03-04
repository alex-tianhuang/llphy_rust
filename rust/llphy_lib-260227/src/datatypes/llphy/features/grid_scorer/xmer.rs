//! Module defining [`XmerIndexableArray`] and [`XmerSize`].
use std::{num::NonZero, ops::Index};
use crate::datatypes::MAX_XMER;

/// A newtype wrapper that indicates something is
/// accessible via an index in `1..=MAX_XMER`.
///
/// See also [`XmerSize`].
pub struct XmerIndexableArray<T>([T; MAX_XMER]);
/// A newtype wrapper that indicates a number is
/// in `1..=MAX_XMER`. For [`XmerIndexableArray`].
#[derive(Clone, Copy)]
pub struct XmerSize(NonZero<usize>);
/// Like `1..=MAX_XMER` but returns [`XmerSize`]s.
pub fn xmer_sizes() -> impl ExactSizeIterator<Item = XmerSize> {
    (0..MAX_XMER).map(|n| unsafe { XmerSize(NonZero::new_unchecked(n + 1)) })
}
impl<T> XmerIndexableArray<T> {
    pub const fn new(arr: [T; MAX_XMER]) -> Self {
        Self(arr)
    }
}
impl<T> Index<XmerSize> for XmerIndexableArray<T> {
    type Output = T;
    fn index(&self, index: XmerSize) -> &Self::Output {
        self.0.index(index.get() - 1)
    }
}
impl<'a, T> IntoIterator for &'a XmerIndexableArray<T> {
    type IntoIter = std::slice::Iter<'a, T>;
    type Item = &'a T;
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}
impl XmerSize {
    /// Make a new [`XmerSize`] from a plain number.
    pub fn new(n: usize) -> Option<Self> {
        if (1..=MAX_XMER).contains(&n) {
            Some(unsafe { XmerSize(NonZero::new_unchecked(n)) })
        } else {
            None
        }
    }
    /// Get a plain number from an [`XmerSize`].
    pub fn get(&self) -> usize {
        self.0.get()
    }
}
