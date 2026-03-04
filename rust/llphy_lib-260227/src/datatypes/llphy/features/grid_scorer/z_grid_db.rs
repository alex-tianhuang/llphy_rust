//! Module defining [`ZGridDB`].
use std::ops::{Deref, DerefMut};

use crate::datatypes::{AAMap, llphy::features::grid_scorer::xmer::XmerIndexableArray};

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
pub struct ZGridSubtable<'a>(&'a mut [(f64, &'a mut [(f64, ZGridDBEntry)])]);
/// Weights for features `a` and `b`,
/// for each `(aa, xmer, zscore_a, zscore_b)` tuple.
pub struct ZGridDBEntry {
    pub weight_total: i64,
    pub weight_a: i64,
    pub weight_b: i64,
}
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
impl ZGridSubtable<'_> {
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
            .binary_search_by(|probe| probe.0.total_cmp(&key_b))
            .ok()?;
        Some(&subarr[j].1)
    }
    /// Find the gridpoint that minimizes the sum of squared
    /// differences to the given zscores.
    fn lookup_thorough(&self, zscore_a: f64, zscore_b: f64) -> &ZGridDBEntry {
        let mut best_entry: Option<(f64, &ZGridDBEntry)> = None;
        let (Ok(start_index) | Err(start_index)) =
            self.0.binary_search_by(|(gridpoint, _)| gridpoint.total_cmp(&zscore_a));
        // I assume this gets optimized out
        let to_sqr_delta_a = |x: f64| {
            let delta_a = x - zscore_a;
            delta_a * delta_a
        };
        let next_search_down = |idx: usize| {
            idx.checked_sub(1)
                .map(|i| (i, to_sqr_delta_a(self.0[i].0)))
        };
        let next_search_up = |idx: usize| {
            Some(idx).and_then(|i| self.0.get(i).map(|(gridpoint, _)| (i, to_sqr_delta_a(*gridpoint))))
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
                    return entry
                }
            }
            let (_, ref subarr) = self.0[i];
            debug_assert!(!subarr.is_empty());
            let (Ok(j) | Err(j)) =
                subarr.binary_search_by(|(gridpoint, _)| gridpoint.total_cmp(&zscore_b));
            let mut check_j = |j: usize| {
                if let Some((gridpoint, _)) = subarr.get(j) {
                    let delta_b = *gridpoint - zscore_b;
                    let sqr_delta_b = delta_b * delta_b;
                    let sqr_delta = sqr_delta_a + sqr_delta_b;
                    let new_best: bool;
                    if let Some(entry) = best_entry {
                        new_best = sqr_delta < entry.0;
                    } else {
                        new_best = true;
                    }
                    if new_best {
                        best_entry = Some((sqr_delta, &subarr[j].1))
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
