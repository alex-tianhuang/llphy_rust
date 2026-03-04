pub use thresholds::Thresholds;
pub use grid_scorer::GridScorer;
mod thresholds;
mod grid_scorer;
pub struct LLPhyFeature<'a> {
    grid_scorer: GridScorer<'a>,
    thresholds: Thresholds,
    sign: i8,
}