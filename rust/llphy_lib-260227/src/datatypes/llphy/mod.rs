mod features;
mod model_train_base;
mod score_type;
pub use features::{FeatureMatrix, PostProcessedFeatureMatrix, ReferenceFeatureMatrix};
pub use model_train_base::ModelTrainingBase;
pub use score_type::ScoreType;
/// The maximum residue separation that Cai collected residue statistics for.
pub(crate) const MAX_XMER: usize = 40;
/// The default features calculated by Cai's old LLPhyScore program.
pub(crate) const DEFAULT_FEATURES: &'static [&'static str] = &[
    "protein-water",
    "protein-carbon",
    "hydrogen bond (long-range)",
    "pi-pi (long-range)",
    "disorder (long)",
    "K-Beta similarity",
    "disorder (short)",
    "electrostatic (short-range)",
];
/// The default pair names and feature names that Cai's
/// old LLPhyScore program used to organize computations.
const PAIR_NAMES_AND_FEATURE_NAMES: &'static [(&'static str, [&'static str; 2])] = &[
    ("S2.SUMPI", ["pi-pi (short-range)", "pi-pi (long-range)"]),
    ("S3.WATER.V2", ["protein-water", "protein-carbon"]),
    (
        "S4.SSPRED",
        ["sec. structure (helices)", "sec. structure (strands)"],
    ),
    ("S5.DISO", ["disorder (long)", "disorder (short)"]),
    (
        "S6.CHARGE.V2",
        ["electrostatic (short-range)", "electrostatic (long-range)"],
    ),
    (
        "S7.ELECHB.V2",
        ["hydrogen bond (short-range)", "hydrogen bond (long-range)"],
    ),
    (
        "S8.CationPi.V2",
        ["cation-pi (short-range)", "cation-pi (long-range)"],
    ),
    (
        "S9.LARKS.V2",
        ["K-Beta similarity", "K-Beta non-similarity"],
    ),
];
/// Utility method for looking up in the above data struct
/// by using a `feature_name` (not a `pair_name`).
pub fn find_pair_and_features_from_one_feature_name(feature_name: &str) -> Option<(&'static str, [&'static str; 2])> {
    PAIR_NAMES_AND_FEATURE_NAMES.iter().find(|(_, known_names)| {
        known_names.contains(&feature_name)
    }).copied()
}