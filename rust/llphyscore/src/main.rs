//! Program that computes phase separation propensity and related
//! biophysical features of sequences from a given input file.
//!
//! Based on the [python library] from Julie Forman-Kay's lab,
//! written by Hao Cai (@haocai1992).
//!
//! [python library]: https://github.com/julie-forman-kay-lab/LLPhyScore
use crate::{fasta::read_fasta, output::write_output};
use anyhow::Error;
use bumpalo::Bump;
use clap::Parser;
use llphyscore_core::{
    datatypes::{ModelTrainingBase, ScoreType},
    featurizer::Featurizer,
    load_pkg_data::load_post_processor,
};
use std::path::{Path, PathBuf};
mod fasta;
mod output;
mod utils;
/// Program that computes phase separation propensity and related
/// biophysical features of sequences from a given input file.
///
/// Based on the [python library] from Julie Forman-Kay's lab,
/// written by Hao Cai (@haocai1992).
///
/// [python library]: https://github.com/julie-forman-kay-lab/LLPhyScore
#[derive(Parser)]
#[clap(verbatim_doc_comment, name = "llphyscore")]
pub struct Args {
    /// Input file, in FASTA format.
    #[arg(short, long)]
    input_file: PathBuf,
    /// Output file, a CSV.
    #[arg(short, long)]
    output_file: Option<PathBuf>,
    /// The type of score to output.
    #[arg(short, long, default_value = "percentile", value_enum)]
    score_type: ScoreType,
    /// The phase separation model to use,
    /// categorized by the negative set it was trained on.
    #[arg(short, long, default_value = "human+PDB", value_enum)]
    model_train_base: ModelTrainingBase,
    /// When given, put all sequence validation errors
    /// into this file instead of clouding `stderr`.
    ///
    /// Creates the file if does not exist,
    /// appends to the file if it does.
    #[arg(long)]
    log_seq_errs_to: Option<PathBuf>,
    /// When set, turn off the progress bar.
    #[arg(long, action)]
    disable_pbar: bool,
    /// When set, sets `--disable-pbar` and `--log-seq-errs-to`
    /// to the equivalent of `/dev/null`.
    #[arg(short, long, action)]
    quiet: bool,
}
const PKG_DATA_ROOT: &'static str = env!("PKG_DATA_ROOT");
fn main() -> Result<(), Error> {
    let Args {
        input_file,
        output_file,
        score_type,
        model_train_base,
        log_seq_errs_to,
        disable_pbar,
        quiet: quiet_override,
    } = Args::parse();
    let arena = Bump::new();
    let sequences = {
        let mut err_out = fasta::alloc_error_handler(log_seq_errs_to, quiet_override, &arena)?;
        read_fasta(&input_file, &arena, &mut *err_out)?
    };
    let featurizer = Featurizer::load_new(Path::new(PKG_DATA_ROOT), model_train_base, &arena)?
        .with_pbar(!(disable_pbar || quiet_override));
    let post_processor = load_post_processor(Path::new(PKG_DATA_ROOT), score_type, model_train_base, &arena)?;
    let matrix = featurizer.featurize(sequences.iter().map(|ent| ent.sequence), &arena)?;
    let matrix = post_processor.post_process(matrix, &arena)?;
    write_output(output_file, sequences, matrix, &arena)
}
