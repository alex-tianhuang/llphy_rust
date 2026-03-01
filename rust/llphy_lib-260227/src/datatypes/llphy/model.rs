//! Module defining structs to map residue-level feature grids to
//! sequence-level features (currently via [`LLPhyFeature`]).
use crate::datatypes::{AAMap, FeatureGridEntry};

/// Thresholds and signs for mapping short-range
/// and long-range feature-grids to sequence-level feature values.
pub struct LLPhyFeature {
    thresholds: AAMap<LLPhyThresholds>,
    signs: LLPhySigns,
}
/// Bounds for the short-range and long-range feature grids.
/// 
/// Each tuple represents an `(upper, lower)` threshold,
/// where being greater than the `upper` threshold is `+1`
/// to the score, and being smaller than the `lower` threshold
/// is `-1` to the score.
struct LLPhyThresholds {
    /// Short-range thresholds.
    sr: (f64, f64),
    /// Long-range thresholds.
    lr: (f64, f64),
}
/// The sign to apply to short-range or long-range features.
struct LLPhySigns {
    sr: bool,
    lr: bool,
}

impl LLPhyFeature {
    /// Get either the short-range or long-range feature values
    /// given a `grid_row` (feature-grid for one feature).
    pub fn get_g2w_score_for_subfeature<const SR: bool>(
        &self,
        grid_row: &AAMap<&[FeatureGridEntry]>,
    ) -> isize {
        let mut sum_score = 0;
        for (weight, &res_scores) in self.thresholds.values().zip(grid_row.values()) {
            let (upper, lower) = if SR { weight.sr } else { weight.lr };
            for entry in res_scores {
                let grid_value = if SR { entry.sr } else { entry.lr };
                sum_score += (grid_value > upper) as isize;
                sum_score -= (grid_value < lower) as isize;
            }
        }
        let sign = if SR { self.signs.sr } else { self.signs.lr };
        sum_score * sign as isize
    }
}
