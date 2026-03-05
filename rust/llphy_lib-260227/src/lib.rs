#![cfg_attr(feature = "simd", feature(portable_simd))]
//! A library for computing phase separation propensity
//! and related biophysical features of sequences.
//!
//! Based on the [python library] from Julie Forman-Kay's lab,
//! written by Hao Cai (@haocai1992).
//!
//! [python library]: https://github.com/julie-forman-kay-lab/LLPhyScore
use std::{
    env::args_os,
    mem::{self, ManuallyDrop},
    path::PathBuf,
};
mod io;
mod load_pkg_data;
mod post_processor;
use crate::{
    datatypes::{DEFAULT_FEATURES, ModelTrainingBase, ScoreType},
    fasta::read_fasta,
    featurizer::featurize,
    load_pkg_data::{load_grid_decoders, load_post_processor},
    output::write_output,
};
use anyhow::{Context, Error};
use bumpalo::{Bump, collections::Vec};
use clap::Parser;
mod datatypes;
mod fasta;
mod featurizer;
mod output;
mod pbar;
use pyo3::{
    Bound, FromPyObject, PyResult, Python,
    prelude::pymodule,
    pyfunction,
    types::{PyModule, PyModuleMethods},
    wrap_pyfunction,
};
/// Program that computes phase separation propensity and related
/// biophysical features of sequences from a given input file.
///
/// Based on the [python library] from Julie Forman-Kay's lab,
/// written by Hao Cai (@haocai1992).
///
/// [python library]: https://github.com/julie-forman-kay-lab/LLPhyScore
#[derive(Parser, FromPyObject)]
#[clap(verbatim_doc_comment, name = "llphyscore")]
pub struct Args {
    /// Input file, in FASTA format.
    #[arg(short, long)]
    input_file: PathBuf,
    /// Output file, a CSV.
    #[arg(short, long)]
    #[pyo3(default = None)]
    output_file: Option<PathBuf>,
    /// The type of score to output.
    #[arg(short, long, default_value = "percentile", value_enum)]
    #[pyo3(default = ScoreType::Percentile)]
    score_type: ScoreType,
    /// The phase separation model to use,
    /// categorized by the negative set it was trained on.
    #[arg(short, long, default_value = "human+PDB", value_enum)]
    #[pyo3(default = ModelTrainingBase::HumanPDB)]
    model_train_base: ModelTrainingBase,
    /// When set, suppress the progress bar and warnings
    /// about sequences that contained invalid characters
    /// (which may lead to unexpected sequences).
    #[arg(short, long, action)]
    #[pyo3(default = false)]
    quiet: bool,
}
#[pymodule(name = "_rust")]
fn _module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(main_cli, m)?)?;
    m.add_function(wrap_pyfunction!(run_fasta_scorer, m)?)?;
    Ok(())
}
/// Entry point for the `cli` that computes biophysical features.
/// 
/// Wrapper around the library function [`run_fasta_scorer`].
#[pyfunction]
fn main_cli(py: Python) -> Result<(), Error> {
    // For some reason `args_os()` gets an extra copy
    // of argument 0 so we skip that.
    let args = Args::parse_from(args_os().skip(1));
    run_fasta_scorer(py, args)
}
/// Loads sequences from `args.input_file`,
/// computes biophysical feature vectors and then
/// computes the phase separation propensity,
/// reporting it to `args.output_file`.
///
/// This function is based on the `run_fasta_scorer`
/// function in the standalone LLPhyScore executable
/// found in Julie Forman-Kay's lab's [python package],
/// written by Hao Cai (@haocai1992).
///
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
///
/// IO
/// --
/// This function calls [`read_fasta`] and [`featurize`],
/// which both touch the stderr output. It also prints a couple
/// messages to stderr.
#[pyfunction]
fn run_fasta_scorer(py: Python, args: Args) -> Result<(), Error> {
    /// Print or not depending on the local variable `quiet`.
    macro_rules! eprintln_or_quiet {
        ($quiet:ident, $($tt:tt)*) => {
            if !$quiet {
                std::eprintln!($($tt)*);
            }
        };
    }
    let Args {
        input_file,
        output_file,
        score_type,
        model_train_base,
        quiet
    } = args;
    // Everything de-allocates when this does,
    // so no need to de-allocate anything else!
    let arena = Bump::new();
    eprintln_or_quiet!(quiet, "Reading fasta file...");
    let sequences = leak_vec(maybe_quiet!(quiet, read_fasta(&input_file, &arena)).with_context(|| {
        format!(
            "failed to get sequences from fasta file @ {}",
            input_file.display()
        )
    })?);
    eprintln_or_quiet!(quiet, "Loading feature grid decoders...");
    let grid_decoders = load_grid_decoders(model_train_base, &arena, py)?;
    eprintln_or_quiet!(quiet, "Loading reference feature distribution...");
    let post_processor = load_post_processor(score_type, model_train_base, &arena)?;
    let matrix = maybe_quiet!(quiet, featurize(sequences, DEFAULT_FEATURES, grid_decoders, &arena, py))?;
    let matrix = post_processor.post_process(matrix, &arena)?;
    write_output(output_file, sequences, matrix, &arena)
}

/// Leak a [`Vec`] managed by an arena into a slice
/// that lives for as long as the arena does not reset.
fn leak_vec<'a, T>(buf: Vec<'a, T>) -> &'a mut [T] {
    let mut buf = ManuallyDrop::new(buf);
    // SAFETY: the arena must de-allocate before this vec does.
    //         Also, do to the use of `ManuallyDrop`, these bytes
    //         are not de-allocated and reused by the arena.
    unsafe { mem::transmute::<&mut [T], &'a mut [T]>(buf.as_mut_slice()) }
}
/// Shorthand for setting a function to
/// quiet or not depending on a non-const value.
/// 
/// Equivalent to:
/// ```
/// if quiet {foo::<QUIET: true>(...)} else {foo::<QUIET: false>(...)}
/// ```
macro_rules! maybe_quiet {
    ($quiet:ident, $func:ident($($args:tt)*)) => {
        if $quiet {$func::<true>($($args)*)} else {$func::<false>($($args)*)}
    };
}
use maybe_quiet;