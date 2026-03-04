pub use grid_scorer::{GridScorer, GridScore};
use grid_decoder::GridDecoder;
mod grid_decoder;
mod grid_scorer;
mod thresholds;
pub struct LLPhyFeaturePair<'a> {
    grid_scorer: GridScorer<'a>,
    decoder_a: Option<GridDecoder>,
    decoder_b: Option<GridDecoder>,
}
