//! Module defining [`ZGridDB`].
use std::ops::{AddAssign, Deref, DerefMut};
use crate::datatypes::AAMap;
use crate::featurizer::grid_scorer::xmer::XmerIndexableArray;
use std::simd::i64x4;

/// A struct of tables indexable by `(aa, xmer)` keys,
/// where each subtable is 2D `zscore`-indexable.
/// See [`ZGridSubtable`] and [`ZGridDBEntry`].
///
/// The subtable is expected to be indexed by `zscore_a`
/// and then `zscore_b`, NOT `zscore_b` then `zscore_a`.
pub struct ZGridDB<'a>(AAMap<XmerIndexableArray<ZGridSubtable<'a>>>);
/// A `(zscore_a, zscore_b)`-indexable collection of weights
/// for features `a` and `b`.
///
/// For quick lookup, the two layers of slices are expected to be
/// sorted by `f64`, the `gridpoint`s associated to each index.
pub struct ZGridSubtable<'a>(&'a [(f64, &'a [ZGridDBEntry])]);
/// Weights for features `a` and `b`,
/// for each `(aa, xmer, zscore_a, zscore_b)` tuple.
/// 
/// The field names are not visible in the SIMD representation,
/// but it is essentially an array consisting of named integers
/// [`weight_total`], [`weight_a`], and [`weight_b`].
/// 
/// [`gridpoint`] (due to SIMD needing four slots) must be
/// represented as a integer but is actually an `f64`.
/// 
/// [`gridpoint`]: Self::gridpoint
/// [`weight_total`]: Self::weight_total
/// [`weight_a`]: Self::weight_a
/// [`weight_b`]: Self::weight_b
pub struct ZGridDBEntry(i64x4);
impl<'a> Deref for ZGridDB<'a> {
    type Target = AAMap<XmerIndexableArray<ZGridSubtable<'a>>>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl<'a> DerefMut for ZGridDB<'a> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl<'a> ZGridDB<'a> {
    /// Wrap the inner field with a `ZGridDB`.
    pub fn new(inner: AAMap<XmerIndexableArray<ZGridSubtable<'a>>>) -> Self {
        Self(inner)
    }
}
impl<'a> ZGridSubtable<'a> {
    /// Make a new [`ZGridSubtable`] with no content.
    pub const fn placeholder() -> Self {
        Self(&mut [])
    }
    /// Make a new [`ZGridSubtable`] with the given data.
    pub fn new(data: &'a [(f64, &'a [ZGridDBEntry])]) -> Self {
        Self(data)
    }
    /// Get the entry best associated with the
    /// given `zscore`s for features `a` and `b`.
    pub fn lookup(&self, zscore_a: f64, zscore_b: f64) -> &ZGridDBEntry {
        if let Some(entry) = self.lookup_quick(zscore_a, zscore_b) {
            entry
        } else {
            self.lookup_thorough(zscore_a, zscore_b)
        }
    }
    /// Snap the zscores to a `0.5`-unit-size grid and see
    /// if there exists an exactly matching entry.
    fn lookup_quick(&self, zscore_a: f64, zscore_b: f64) -> Option<&ZGridDBEntry> {
        let key_a = snap_to_grid(zscore_a);
        let key_b = snap_to_grid(zscore_b);
        let i = self
            .0
            .binary_search_by(|probe| probe.0.total_cmp(&key_a))
            .ok()?;
        let (_, ref subarr) = self.0[i];
        let j = subarr
            .binary_search_by(|probe| probe.gridpoint().total_cmp(&key_b))
            .ok()?;
        Some(&subarr[j])
    }
    /// Find the gridpoint that minimizes the sum of squared
    /// differences to the given zscores.
    fn lookup_thorough(&self, zscore_a: f64, zscore_b: f64) -> &ZGridDBEntry {
        let mut best_entry: Option<(f64, &ZGridDBEntry)> = None;
        let (Ok(start_index) | Err(start_index)) = self
            .0
            .binary_search_by(|(gridpoint, _)| gridpoint.total_cmp(&zscore_a));
        // I assume this gets optimized out
        let to_sqr_delta_a = |x: f64| {
            let delta_a = x - zscore_a;
            delta_a * delta_a
        };
        let next_search_down =
            |idx: usize| idx.checked_sub(1).map(|i| (i, to_sqr_delta_a(self.0[i].0)));
        let next_search_up = |idx: usize| {
            Some(idx).and_then(|i| {
                self.0
                    .get(i)
                    .map(|(gridpoint, _)| (i, to_sqr_delta_a(*gridpoint)))
            })
        };
        let mut down_search = next_search_down(start_index);
        let mut up_search = next_search_up(start_index);
        debug_assert!(self.0.len() > 0);
        for _ in 0..self.0.len() {
            let (i, sqr_delta_a) = match (down_search, up_search) {
                (Some((down_ptr, down_sqr_delta_a)), Some((up_ptr, up_sqr_delta_a))) => {
                    if down_sqr_delta_a < up_sqr_delta_a {
                        down_search = next_search_down(down_ptr);
                        (down_ptr, down_sqr_delta_a)
                    } else {
                        up_search = next_search_up(up_ptr + 1);
                        (up_ptr, up_sqr_delta_a)
                    }
                }
                (None, Some((up_ptr, up_sqr_delta_a))) => {
                    up_search = next_search_up(up_ptr + 1);
                    (up_ptr, up_sqr_delta_a)
                }
                (Some((down_ptr, down_sqr_delta_a)), None) => {
                    down_search = next_search_down(down_ptr);
                    (down_ptr, down_sqr_delta_a)
                }
                (None, None) => unreachable!(),
            };
            if let Some((best_sqr_delta, entry)) = best_entry {
                if best_sqr_delta <= sqr_delta_a {
                    return entry;
                }
            }
            let (_, ref subarr) = self.0[i];
            debug_assert!(!subarr.is_empty());
            let (Ok(j) | Err(j)) =
                subarr.binary_search_by(|entry| entry.gridpoint().total_cmp(&zscore_b));
            let mut check_j = |j: usize| {
                if let Some(entry) = subarr.get(j) {
                    let delta_b = entry.gridpoint() - zscore_b;
                    let sqr_delta_b = delta_b * delta_b;
                    let sqr_delta = sqr_delta_a + sqr_delta_b;
                    let new_best: bool;
                    if let Some(entry) = best_entry {
                        new_best = sqr_delta < entry.0;
                    } else {
                        new_best = true;
                    }
                    if new_best {
                        best_entry = Some((sqr_delta, &subarr[j]))
                    }
                }
            };
            check_j(j);
            let Some(j) = j.checked_sub(1) else {
                continue;
            };
            check_j(j);
        }
        best_entry.unwrap().1
    }
}
/// Shorthand for associating each index of the
/// array of [`ZGridDBEntry`] with a named field.
macro_rules! impl_getters {
    ($([$index:literal, $field:ident]),*) => {
        impl ZGridDBEntry {
            $(pub fn $field(&self) -> i64 {
                self.0.as_array()[$index]
            })*
        }
    };
}
impl_getters!([1, weight_total], [2, weight_a], [3, weight_b]);
impl ZGridDBEntry {
    /// Get the z-score value of feature `B`
    /// that this entry is associated with.
    pub fn gridpoint(&self) -> f64 {
        f64::from_bits(self.0.as_array()[0].cast_unsigned())
    }
    /// Get a new [`ZGridDBEntry`] with all `0.0_f64`.
    pub fn new_zeroed() -> Self {
        Self(i64x4::from_array([0; 4]))
    }
    /// Get a new [`ZGridDBEntry`] from its four fields.
    pub fn from_parts(gridpoint: f64, weight_total: i64, weight_a: i64, weight_b: i64) -> Self {
        let gridpoint = f64::to_bits(gridpoint).cast_signed();
        Self(i64x4::from_array([gridpoint, weight_total, weight_a, weight_b]))
    }
    /// Utility method for [`super::GridScorer::score_sequence`].
    pub fn freq_a(&self) -> f64 {
        if self.weight_total() == 0 {
            0.0
        } else {
            self.weight_a() as f64 / self.weight_total() as f64
        }
    }
    /// Utility method for [`super::GridScorer::score_sequence`].
    pub fn freq_b(&self) -> f64 {
        if self.weight_total() == 0 {
            0.0
        } else {
            self.weight_b() as f64 / self.weight_total() as f64
        }
    }
}
impl AddAssign<&Self> for ZGridDBEntry {
    fn add_assign(&mut self, rhs: &Self) {
        self.0 += &rhs.0;
    }
}
/// Snap a `zscore` to the grid (half-integers).
///
/// Ripped from `make_linekey` in Julie Forman-Kay's lab's
/// standalone LLPhyScore executable in their [python package],
/// written by Hao Cai (@haocai1992).
///
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
pub fn snap_to_grid(zscore: f64) -> f64 {
    ((zscore * 2.0).round() / 2.0).clamp(-8.0, 12.0)
}
