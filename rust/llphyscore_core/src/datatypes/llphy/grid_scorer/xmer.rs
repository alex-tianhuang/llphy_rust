//! Module defining [`XmerIndexableArray`] and [`XmerSize`].
use std::{num::NonZero, ops::{Index, IndexMut}};
use borsh::{BorshDeserialize, BorshSerialize};
use crate::datatypes::MAX_XMER;

/// A newtype wrapper that indicates something is
/// accessible via an index in `1..=MAX_XMER`.
///
/// See also [`XmerSize`].
#[derive(BorshDeserialize, BorshSerialize, PartialEq)]
#[repr(transparent)]
pub struct XmerIndexableArray<T>([T; MAX_XMER]);
/// A newtype wrapper that indicates a number is
/// in `1..=MAX_XMER`. For [`XmerIndexableArray`].
#[derive(Clone, Copy)]
pub struct XmerSize(NonZero<usize>);
impl<T> XmerIndexableArray<T> {
    /// Constructor.
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
impl<T> IndexMut<XmerSize> for XmerIndexableArray<T> {
    fn index_mut(&mut self, index: XmerSize) -> &mut Self::Output {
        self.0.index_mut(index.get() - 1)
    }
}
impl<'a, T> IntoIterator for &'a XmerIndexableArray<T> {
    type IntoIter = std::slice::Iter<'a, T>;
    type Item = &'a T;
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}
impl<'a, T> IntoIterator for &'a mut XmerIndexableArray<T> {
    type IntoIter = std::slice::IterMut<'a, T>;
    type Item = &'a mut T;
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}
impl XmerSize {
    /// Make a new [`XmerSize`] from a plain number.
    pub fn new(n: usize) -> Option<Self> {
        if (1..=MAX_XMER).contains(&n) {
            Some(unsafe { Self::new_unchecked(n) })
        } else {
            None
        }
    }
    /// Make a new [`XmerSize`] from a plain number
    /// without checking bounds.
    pub unsafe fn new_unchecked(n: usize) -> Self {
        unsafe {Self(NonZero::new_unchecked(n))}
    }
    /// Get a plain number from an [`XmerSize`].
    pub fn get(&self) -> usize {
        self.0.get()
    }
}
