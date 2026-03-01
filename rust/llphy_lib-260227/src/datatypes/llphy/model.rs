use crate::datatypes::{AAMap, FeatureGrid, FeatureGridEntry};

pub struct LLPhyFeatureOld {
    weights: AAMap<LLPhyWeightOld>,
    signs: LLPhySignOld,
}
pub struct LLPhyWeightOld {
    sr: (f64, f64),
    lr: (f64, f64),
}
pub struct LLPhySignOld {
    sr: bool,
    lr: bool,
}

impl LLPhyFeatureOld {
    pub fn get_sr_lr_g2w_score(&self, grid_row: &FeatureGrid<'_>) -> (isize, isize) {
        let sr_score = self.get_g2w_score_for_subfeature::<true>(grid_row);
        let lr_score = self.get_g2w_score_for_subfeature::<false>(grid_row);
        (sr_score, lr_score)
    }
    fn get_g2w_score_for_subfeature<const SR: bool>(
        &self,
        grid_row: &AAMap<&[FeatureGridEntry]>,
    ) -> isize {
        let mut sum_score = 0;
        for (weight, &res_scores) in self.weights.values().zip(grid_row.values()) {
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
