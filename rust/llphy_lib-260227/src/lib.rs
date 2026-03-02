//! A library for computing phase separation propensity
//! and related biophysical features of sequences.
//!
//! Based on the [python library] from Julie Forman-Kay's lab,
//! written by Hao Cai (@haocai1992).
//!
//! [python library]: https://github.com/julie-forman-kay-lab/LLPhyScore
use std::{
    fmt::Write,
    fs::File,
    io::{StdoutLock, stdout},
    mem::{self, ManuallyDrop},
    path::PathBuf,
};
mod io;
use crate::{
    datatypes::{
        AAMap, FastaEntry, FeatureGrid, GridScorer, LLPhyFeature, ModelTrainingType, PostProcessor,
        ScoreType,
    },
    fasta::read_fasta,
    load_legacy::{
        LEGACY_FEATURE_NAMES, load_legacy_model_and_feature_order, load_legacy_pdb_statistics,
        load_legacy_reference_g2w_scores,
    },
};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use clap::Parser;
mod datatypes;
mod fasta;
mod load_legacy;
use indicatif::{self, ProgressBar};
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
    output_file: Option<PathBuf>,
    /// The type of score to output.
    #[arg(short, long, default_value = "percentile", value_enum)]
    #[pyo3(default = ScoreType::Percentile)]
    score_type: ScoreType,
    /// The phase separation model to use,
    /// categorized by the negative set it was trained on.
    #[arg(short, long, default_value = "human+PDB", value_enum)]
    #[pyo3(default = ModelTrainingType::HumanPDB)]
    model_train_base: ModelTrainingType,
}
#[pymodule(name = "_rust")]
fn _module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_fasta_scorer, m)?)?;
    Ok(())
}
/// Loads sequences from `args.input_file`,
/// computes biophysical feature grids and then
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
/// This function calls [`read_fasta`] and [`seqs_to_grids`],
/// which both touch the stderr output. It also prints one message to stderr.
#[pyfunction]
fn run_fasta_scorer(py: Python, args: Args) -> Result<(), Error> {
    let Args {
        input_file,
        output_file,
        score_type,
        model_train_base,
    } = args;
    // Everything de-allocates when this does,
    // so no need to de-allocate anything else!
    let arena = Bump::new();
    let (pdb_statistics_names, pdb_statistics_scorer) = load_named_pdb_statistics(&arena)?;
    let feature_names = leak_vec(load_feature_names(&arena));
    feature_names.sort_by_key(|name| {
        feature_traversal_order(name, pdb_statistics_names)
    });
    let post_processor = load_post_processor(score_type, model_train_base, &arena)?;
    let model = load_llphy_model_and_feature_order(model_train_base, &arena, py)?;
    model.sort_by_key(|(name, _)| {
        pdb_statistics_names.iter().enumerate().find(|(_, (known_name, _))| known_name == name).map(|(idx, _)|idx).unwrap_or(usize::MAX)
    });
    let model = arena.alloc_slice_fill_iter(model.iter().map(|(_, feature)| feature.clone()));
    let sequences = leak_vec(read_fasta(input_file, &arena)?);
    let grids = leak_vec(seqs_to_grids::<true>(
        sequences,
        pdb_statistics_scorer,
        pdb_statistics_names,
        feature_names,
        &arena,
    ));
    eprintln!("CALCULATING SCORES");
    let output_scores = leak_vec(get_g2w_scores(
        grids,
        model,
        pdb_statistics_names,
        feature_names,
        &arena,
    ));
    post_process(post_processor, output_scores);
    write_output(output_file, output_scores, sequences, feature_names, &arena)?;
    Ok(())
}

/// Placeholder for some feature names.
///
/// In the standalone LLPhyScore executable found in
/// Julie Forman-Kay's lab's [python package] written
/// by Hao Cai (@haocai1992), feature names exist
/// as a global hard coded list.
///
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
fn load_feature_names(arena: &Bump) -> Vec<'_, &str> {
    Vec::from_iter_in(LEGACY_FEATURE_NAMES.iter().cloned(), arena)
}
/// Docs out of date at the moment, but this is a helper function for
/// sorting things by output order.
/// 
/// The traversal order is defined by flattening the `pdb_statistics_names`
/// and filtering out names that do not appear in `feature_names`.
fn feature_traversal_order<'a, 'b>(
    feature_name: &str,
    pdb_statistics_names: &[(&str, [&str; 2])],
) -> (usize, u8) {
    pdb_statistics_names
        .iter()
        .enumerate()
        .find(|(_, (_, [sr_name, lr_name]))| feature_name == *sr_name || feature_name == *lr_name)
        .map(|(slot, (sr_name, _))| {
            if *sr_name == feature_name {
                (slot, 0)
            } else {
                (slot, 1)
            }
        })
        .unwrap_or((usize::MAX, 0))
}
/// Load PDB statistics that help compute [`seqs_to_grids`],
/// and the names of the long and short range statistics.
///
/// At the moment, this replaces all of the global state logic
/// involved in setting up the `GridScore` class instances in the
/// standalone LLPhyScore executable found in Julie Forman-Kay's
/// lab's [python package], written by Hao Cai (@haocai1992).
///
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
fn load_named_pdb_statistics<'a>(
    arena: &'a Bump,
) -> Result<(&'a [(&'a str, [&'a str; 2])], &'a [GridScorer<'a>]), Error> {
    load_legacy_pdb_statistics(arena)
}

/// Load a "post-processor" (see [`PostProcessor`]).
///
/// At the moment, this replaces all of the global state logic
/// involved in setting up the `human_g2w_scores` and
/// `g2w_means_stds_human` dictionaries in the standalone
/// LLPhyScore executable found in Julie Forman-Kay's
/// lab's [python package], written by Hao Cai (@haocai1992).
///
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
fn load_post_processor<'a>(
    score_type: ScoreType,
    model_train_base: ModelTrainingType,
    arena: &'a Bump,
) -> Result<PostProcessor<'a>, Error> {
    match score_type {
        ScoreType::Raw => Ok(PostProcessor::new_raw()),
        ScoreType::ZScore => {
            let ref_scores = load_reference_g2w_scores(model_train_base, arena)?;
            let num_features = ref_scores[0].subfeatures.len();
            Ok(PostProcessor::new_zscore(ref_scores, num_features, arena))
        }
        ScoreType::Percentile => {
            let ref_scores = load_reference_g2w_scores(model_train_base, arena)?;
            let num_features = ref_scores[0].subfeatures.len();
            Ok(PostProcessor::new_percentile(
                ref_scores,
                num_features,
                arena,
            ))
        }
    }
}
/// Helper for [`load_post_processor`].
///
/// Load a reference dataset of biophysical feature scores
/// (currently [`G2WScores`]) that corresponds to the human
/// proteome or a PDB-derived dataset.
fn load_reference_g2w_scores(
    model_train_base: ModelTrainingType,
    arena: &Bump,
) -> Result<&[G2WScores<'_>], Error> {
    load_legacy_reference_g2w_scores(model_train_base, arena)
}

/// Load the ML model that converts residue-level feature grids
/// into sequence-level biophysical feature values in [`get_g2w_scores`].
///
/// At the moment, this replaces all of the global state logic
/// involved in setting up the global `final_models` dictionary
/// in the standalone LLPhyScore executable found in Julie
/// Forman-Kay's lab's [python package], written by Hao Cai
/// (@haocai1992).
///
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
fn load_llphy_model_and_feature_order<'a>(
    model_train_base: ModelTrainingType,
    arena: &'a Bump,
    py: Python,
) -> Result<&'a mut [(&'a str, LLPhyFeature)], Error> {
    load_legacy_model_and_feature_order(model_train_base, arena, py)
}

/// Given sequences and PDB statistics (currently `&[GridScorer]`),
/// compute a sequence x feature grid of biophysical scores.
///
/// May waste some time computing short-range / long-range features
/// even though they are not used later. Will fix someday!
///
/// This function is based on the `seqs2grids`
/// function in the standalone LLPhyScore executable
/// found in Julie Forman-Kay's lab's [python package],
/// written by Hao Cai (@haocai1992).
///
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
///
/// IO
/// --
/// If `PBAR` is true, this function
/// prints the grid that is currently being worked on along with a
/// progress bar to whatever the default terminal for [`indicatif`] is.
fn seqs_to_grids<'a, const PBAR: bool>(
    sequences: &[FastaEntry<'_>],
    pdb_statistics_scorer: &[GridScorer],
    pdb_statistics_names: &[(&str, [&str; 2])],
    feature_names_requested: &[&str],
    arena: &'a Bump,
) -> Vec<'a, &'a mut [FeatureGrid<'a>]> {
    let num_grids = pdb_statistics_names
        .iter()
        .filter(|(_, [sr_name, lr_name])| {
            feature_names_requested.contains(sr_name) || feature_names_requested.contains(lr_name)
        })
        .count();
    let mut grids = Vec::from_iter_in(
        (0..sequences.len()).map(|_| {
            leak_vec(Vec::from_iter_in(
                (0..num_grids).map(|_| AAMap::default()),
                arena,
            ))
        }),
        arena,
    );
    let mut idx = 0;
    for ((grid_name, [sr_name, lr_name]), grid_scorer) in
        pdb_statistics_names.iter().zip(pdb_statistics_scorer)
    {
        if feature_names_requested.contains(sr_name) || feature_names_requested.contains(lr_name) {
            if PBAR {
                let pbar = ProgressBar::new(sequences.len() as _);
                pbar.println(
                    bumpalo::format!(in arena, "CONVERTING SEQUENCES TO {} GRIDS", grid_name),
                );
                for (grid_row, entry) in grids.iter_mut().zip(sequences) {
                    grid_row[idx] = grid_scorer.score_sequence(entry.sequence, arena);
                    pbar.inc(1);
                }
                pbar.abandon();
            } else {
                for (grid_row, entry) in grids.iter_mut().zip(sequences) {
                    grid_row[idx] = grid_scorer.score_sequence(entry.sequence, arena);
                }
            }
            idx += 1;
        }
    }
    grids
}

/// Leak a [`Vec`] managed by an arena into a slice
/// that lives for as long as the arena does not reset.
pub fn leak_vec<'a, T>(buf: Vec<'a, T>) -> &'a mut [T] {
    let mut buf = ManuallyDrop::new(buf);
    // SAFETY: the arena must de-allocate before this vec does.
    //         Also, do to the use of `ManuallyDrop`, these bytes
    //         are not de-allocated and reused by the arena.
    unsafe { mem::transmute::<&mut [T], &'a mut [T]>(buf.as_mut_slice()) }
}
/// A struct that holds all the computed biophysical
/// feature scores and their sum for one sequence.
pub struct G2WScores<'a> {
    feature_sum: f64,
    subfeatures: &'a mut [f64],
}
/// Given residue-level feature grids and an LLPhyScore model
/// (currently `&[LLPhyFeatureOld]`), compute the sequence-level
/// biophysical feature grids (currently [`G2WScores`]).
///
/// This function is based on the `get_g2w_scores`
/// function in the standalone LLPhyScore executable
/// found in Julie Forman-Kay's lab's [python package],
/// written by Hao Cai (@haocai1992).
///
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
fn get_g2w_scores<'a>(
    grids: &[&mut [FeatureGrid<'_>]],
    model: &[LLPhyFeature],
    pdb_statistics_names: &[(&str, [&str; 2])],
    feature_names_requested: &[&str],
    arena: &'a Bump,
) -> Vec<'a, G2WScores<'a>> {
    let mut scores = Vec::with_capacity_in(grids.len(), arena);
    for grid in grids {
        let mut subfeatures = Vec::with_capacity_in(model.len() * 2, arena);
        let mut grid_iter = grid.iter();
        for ((_, [sr_name, lr_name]), subfeature) in pdb_statistics_names.iter().zip(model) {
            let mut grid_row = None;
            if feature_names_requested.contains(sr_name) {
                let sr_feat = subfeature.get_g2w_score_for_subfeature::<true>(
                    grid_row.get_or_insert_with(|| grid_iter.next().unwrap()),
                );
                subfeatures.push(sr_feat as f64);
            }
            if feature_names_requested.contains(lr_name) {
                let lr_feat = subfeature.get_g2w_score_for_subfeature::<false>(
                    grid_row.get_or_insert_with(|| grid_iter.next().unwrap()),
                );
                subfeatures.push(lr_feat as f64);
            }
        }
        let feature_sum = subfeatures.iter().sum::<f64>();
        scores.push(G2WScores {
            feature_sum,
            subfeatures: leak_vec(subfeatures),
        })
    }
    scores
}

/// Convert raw biophysical feature scores (currently [`G2WScores`])
/// to a percentile or z-score if requested by the post-processor loaded.
fn post_process(post_processor: PostProcessor<'_>, output_scores: &mut [G2WScores<'_>]) {
    for output_row in output_scores {
        post_processor.transform(output_row);
    }
}
/// Output the final biophysical scores (currently [`G2WScores`])
/// associated with the given `sequences` to a file at the given
/// `path` or stdout.
fn write_output(
    path: Option<PathBuf>,
    output_scores: &[G2WScores<'_>],
    sequences: &[FastaEntry<'_>],
    feature_names: &[&str],
    arena: &Bump,
) -> Result<(), Error> {
    let mut stdout_writer: csv::Writer<StdoutLock>;
    let mut file_writer: csv::Writer<File>;
    let write_row: &mut dyn FnMut(&[&str]) -> Result<(), csv::Error>;
    match path {
        Some(path) => {
            file_writer = csv::Writer::from_path(path)?;
            write_row = arena.alloc(|row| file_writer.write_record(row))
        }
        None => {
            stdout_writer = csv::Writer::from_writer(stdout().lock());
            write_row = arena.alloc(|row| stdout_writer.write_record(row))
        }
    }
    let mut column_buffer = Vec::with_capacity_in(feature_names.len() + 2, arena);
    let values_buffer = leak_vec(Vec::from_iter_in(
        (0..feature_names.len() + 1).map(|_| bumpalo::format!(in arena, "{}", 2.0_f64.sqrt())),
        arena,
    ));
    let feature_sum_column_name = bumpalo::format!(in arena, "{}-feature sum", feature_names.len());
    column_buffer.push("tag");
    column_buffer.extend_from_slice(feature_names);
    column_buffer.push(unsafe { &*(&*feature_sum_column_name as *const str) });
    write_row(&column_buffer)?;
    for (entry, output_row) in sequences.iter().zip(output_scores) {
        let [tag_slot, feature_slots @ .., feature_sum_slot] = &mut *column_buffer else {
            unreachable!()
        };
        let [subfeatures_buffer @ .., feature_sum_buffer] = values_buffer else {
            unreachable!()
        };
        *tag_slot = entry.header;
        for (value_buffer, value) in subfeatures_buffer
            .iter_mut()
            .zip(output_row.subfeatures.iter())
        {
            value_buffer.clear();
            write!(value_buffer, "{}", value).unwrap();
        }
        for (slot, value_str) in feature_slots.iter_mut().zip(subfeatures_buffer.iter()) {
            *slot = unsafe { &*(&**value_str as *const str) };
        }
        feature_sum_buffer.clear();
        write!(feature_sum_buffer, "{}", output_row.feature_sum).unwrap();
        *feature_sum_slot = unsafe { &*(&**feature_sum_buffer as *const str) };
        write_row(&column_buffer)?;
    }
    Ok(())
}
