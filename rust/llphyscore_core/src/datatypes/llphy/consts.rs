//! Some constants.
/// The maximum residue separation that Cai collected residue statistics for.
pub const MAX_XMER: usize = 40;
/// The default features calculated by Cai's old LLPhyScore program.
pub const DEFAULT_FEATURES: [&'static str; 8] = [
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
pub const PAIR_NAMES_AND_FEATURE_NAMES: [(&'static str, [&'static str; 2]); 8] = [
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