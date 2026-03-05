//! Module defining substructures of [`super::GridScorer`]
//! using `#[portable_simd]`.
pub use avg_sdev_db::AvgSdevDBEntry;
pub use pair_freq_db::{PairFreqDBEntry, PairFreqEntrySum};
pub use z_grid_db::{ZGridDBEntry, ZGridSubtable, ZGridEntrySum};
mod z_grid_db;
mod avg_sdev_db;
mod pair_freq_db;
mod xmer;
