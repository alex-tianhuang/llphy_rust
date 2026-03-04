pub use avg_sdev_db::AvgSdevDB;
pub use pair_freq_db::PairFreqDB;
pub use xmer::{XmerIndexableArray, XmerSize, xmer_sizes};
pub use z_grid_db::ZGridDB;

use crate::datatypes::AAMap;
mod avg_sdev_db;
mod pair_freq_db;
mod xmer;
mod z_grid_db;
/// A struct that contains all the necessary data to
/// make biophysical feature grids ([`GridScore`]s)
/// from sequences.
pub struct GridScorer<'a> {
    pub pair_freqs_n_term_anchored: PairFreqDB,
    pub pair_freqs_c_term_anchored: PairFreqDB,
    pub avg_sdevs: AvgSdevDB,
    pub z_grid: ZGridDB<'a>,
}
/// A collection of biophysical feature scores for
/// each residue in a sequence (not tracked)
pub struct GridScore<'a> {
    pub feature_a_scores: Option<AAMap<&'a [f64]>>,
    pub feature_b_scores: Option<AAMap<&'a [f64]>>
}
