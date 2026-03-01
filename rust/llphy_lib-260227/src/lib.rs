//! A library for computing phase separation propensity
//! and related biophysical features of sequences.
//! 
//! Based on the [python library] from Julie Forman-Kay's lab,
//! written by Hao Cai (@haocai1992).
//! 
//! [python library]: https://github.com/julie-forman-kay-lab/LLPhyScore
use std::{
    fmt::Write, fs::File, io::{StdoutLock, stdout}, mem::{self, ManuallyDrop}, path::PathBuf, ptr
};
mod io;
use crate::{
    datatypes::{
        FastaEntry, GridScorer, LLPhyFeatureOld, PostProcessor, ResScoresOld, ScoreType,
    },
    io::read_fasta,
};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use clap::Parser;
mod datatypes;
mod fasta;

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
}
/// Main method for the binary.
/// 
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
pub fn bin_main(args: Args) -> Result<(), Error> {
    let Args {
        input_file,
        output_file,
        score_type,
    } = args;
    // Everything de-allocates when this does,
    // so no need to de-allocate anything else!
    let arena = Bump::new();
    let feature_names = leak_vec(load_feature_names(&arena));
    let pdb_statistics = leak_vec(load_pdb_statistics(feature_names, &arena));
    let post_processor = load_post_processor(score_type, feature_names, &arena);
    let model = leak_vec(load_llphy_model(feature_names, &arena));
    let sequences = leak_vec(read_fasta(input_file, &arena)?);
    let grids = leak_vec(seqs_to_grids(sequences, pdb_statistics, &arena));
    let output_scores = leak_vec(get_g2w_scores(grids, model, &arena));
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
    todo!()
}
/// Load PDB statistics that help compute [`seqs_to_grids`].
/// `feature_names` should be used to determine the order of
/// features returned. 
/// 
/// At the moment, this replaces all of the global state logic
/// involved in setting up the `GridScore` class instances in the
/// standalone LLPhyScore executable found in Julie Forman-Kay's
/// lab's [python package], written by Hao Cai (@haocai1992).
/// 
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
fn load_pdb_statistics<'a>(feature_names: &[&str], arena: &'a Bump) -> Vec<'a, GridScorer> {
    todo!()
}

/// Load a "post-processor" (see [`PostProcessor`]).
/// `feature_names` should be used to determine the order of
/// features returned. 
/// 
/// At the moment, this replaces all of the global state logic
/// involved in setting up the `human_g2w_scores` and
/// `g2w_means_stds_human` dictionaries in the standalone
/// LLPhyScore executable found in Julie Forman-Kay's
/// lab's [python package], written by Hao Cai (@haocai1992).
/// 
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
fn load_post_processor<'a>(score_type: ScoreType, feature_names: &[&str], arena: &'a Bump) -> PostProcessor<'a> {
    match score_type {
        ScoreType::Raw => PostProcessor::new_raw(),
        ScoreType::ZScore => {
            let ref_scores = load_reference_g2w_scores(arena);
            let num_features = ref_scores[0].subfeatures.len();
            PostProcessor::new_zscore(ref_scores, num_features, arena)
        }
        ScoreType::Percentile => {
            let ref_scores = load_reference_g2w_scores(arena);
            let num_features = ref_scores[0].subfeatures.len();
            PostProcessor::new_percentile(ref_scores, num_features, arena)
        }
    }
}
/// Helper for [`load_post_processor`].
/// 
/// Load a reference dataset of biophysical feature scores
/// (currently [`G2WScores`]) that corresponds to the human
/// proteome or a PDB-derived dataset.
fn load_reference_g2w_scores(arena: &Bump) -> Vec<'_, G2WScores<'_>> {
    todo!()
}

/// Load the ML model that converts residue-level feature grids
/// into sequence-level biophysical feature values in [`get_g2w_scores`].
/// `feature_names` should be used to determine the order of
/// features returned. 
/// 
/// At the moment, this replaces all of the global state logic
/// involved in setting up the global `final_models` dictionary
/// in the standalone LLPhyScore executable found in Julie
/// Forman-Kay's lab's [python package], written by Hao Cai
/// (@haocai1992).
/// 
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
fn load_llphy_model<'a>(feature_names: &[&str], arena: &'a Bump) -> Vec<'a, LLPhyFeatureOld> {
    todo!()
}

/// Given sequences and PDB statistics (currently `&[GridScorer]`),
/// compute a sequence x feature grid of biophysical scores.
///
/// This function is based on the `seqs2grids`
/// function in the standalone LLPhyScore executable
/// found in Julie Forman-Kay's lab's [python package],
/// written by Hao Cai (@haocai1992).
/// 
/// [python package]: https://github.com/julie-forman-kay-lab/LLPhyScore
fn seqs_to_grids<'a>(
    sequences: &[FastaEntry<'_>],
    pdb_statistics: &[GridScorer],
    arena: &'a Bump,
) -> Vec<'a, &'a [ResScoresOld<'a>]> {
    let mut grids = Vec::with_capacity_in(sequences.len(), arena);
    for entry in sequences {
        let mut feature_grid = Vec::with_capacity_in(pdb_statistics.len(), arena);
        for grid_scorer in pdb_statistics {
            let res_scores = grid_scorer.score_sequence_to_res_scores(entry.sequence, arena);
            feature_grid.push(res_scores);
        }
        grids.push(leak_vec(feature_grid) as &[_]);
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
    grids: &[&[ResScoresOld<'_>]],
    model: &[LLPhyFeatureOld],
    arena: &'a Bump,
) -> Vec<'a, G2WScores<'a>> {
    let mut scores = Vec::with_capacity_in(grids.len(), arena);
    for grid in grids {
        debug_assert_eq!(model.len(), grid.len());
        let mut subfeatures = Vec::with_capacity_in(model.len() * 2, arena);
        for (subfeature, grid_row) in model.iter().zip(*grid) {
            let (sr_feat, lr_feat) = subfeature.get_sr_lr_g2w_score(grid_row);
            subfeatures.push(sr_feat as f64);
            subfeatures.push(lr_feat as f64);
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
        },
        None => {
            stdout_writer = csv::Writer::from_writer(stdout().lock());
            write_row = arena.alloc(|row| stdout_writer.write_record(row))
        }
    }
    struct DropWriterWrapper<'a>(&'a mut dyn FnMut(&[&str]) -> Result<(), csv::Error>);
    impl Drop for DropWriterWrapper<'_> {
        fn drop(&mut self) {
            unsafe {
                ptr::drop_in_place(self.0);
            }
        }
    }
    let writer = DropWriterWrapper(write_row);
    let mut column_buffer = Vec::with_capacity_in(feature_names.len() + 2, arena);
    let values_buffer = leak_vec(Vec::from_iter_in(
        (0..feature_names.len() + 1).map(|_| bumpalo::format!(in arena, "{}", 2.0_f64.sqrt())),
        arena,
    ));
    let feature_sum_column_name =
        bumpalo::format!(in arena, "{}-feature sum", feature_names.len());
    column_buffer.push("tag");
    column_buffer.extend_from_slice(feature_names);
    column_buffer.push(unsafe { &*(&*feature_sum_column_name as *const str) });
    writer.0(&column_buffer)?;
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
        write!(feature_sum_buffer, "{}", output_row.feature_sum).unwrap();
        *feature_sum_slot = unsafe { &*(&**feature_sum_buffer as *const str) };
        writer.0(&column_buffer)?;
    }
    Ok(())
}
