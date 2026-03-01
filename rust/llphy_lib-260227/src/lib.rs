use std::{
    fmt::Write, fs::File, io::{StdoutLock, stdout}, mem::{self, ManuallyDrop}, path::PathBuf, ptr
};
mod io;
use crate::{
    datatypes::{
        FastaEntry, GridScoreOld, LLPhyFeatureOld, PostProcessor, ResScoresOld, ScoreType,
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
#[derive(Parser)]
#[clap(verbatim_doc_comment)]
pub struct Args {
    /// Input file, in FASTA format.
    #[arg(short, long)]
    input_file: PathBuf,
    /// Output file, a CSV.
    #[arg(short, long)]
    output_file: Option<PathBuf>,
    #[arg(short, long, default_value = "percentile")]
    score_type: ScoreType,
}
pub fn run_fasta_scorer(args: Args) -> Result<(), Error> {
    let Args {
        input_file,
        output_file,
        score_type,
    } = args;
    let arena = Bump::new();
    let grid_scorer = leak_vec(load_grid_scorer(&arena));
    let post_processor = load_post_processor(score_type, &arena);
    let model = leak_vec(load_llphy_model(&arena));
    let sequences = leak_vec(read_fasta(input_file, &arena)?);
    let feature_names = leak_vec(load_feature_names(&arena));
    let grids = leak_vec(seqs_to_grids_old(&sequences, &grid_scorer, &arena));
    let output_scores = leak_vec(get_g2w_scores(&grids, &model, &arena));
    post_process(post_processor, output_scores);
    write_output(output_file, output_scores, sequences, feature_names, &arena)?;
    Ok(())
}

fn seqs_to_grids_old<'a>(
    sequences: &[FastaEntry<'_>],
    grid_scores: &[GridScoreOld],
    arena: &'a Bump,
) -> Vec<'a, &'a [ResScoresOld<'a>]> {
    let mut grids = Vec::with_capacity_in(sequences.len(), arena);
    for entry in sequences {
        let mut feature_grid = Vec::with_capacity_in(grid_scores.len(), arena);
        for grid_score in grid_scores {
            let res_scores = grid_score.score_sequence_to_res_scores(entry.sequence, arena);
            feature_grid.push(res_scores);
        }
        grids.push(leak_vec(feature_grid) as &[_]);
    }
    grids
}

fn load_grid_scorer(arena: &Bump) -> Vec<'_, GridScoreOld> {
    todo!()
}

fn load_llphy_model(arena: &Bump) -> Vec<'_, LLPhyFeatureOld> {
    todo!()
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
pub struct G2WScoresOld<'a> {
    feature_sum: f64,
    subfeatures: &'a mut [f64],
}
fn get_g2w_scores<'a>(
    grids: &[&[ResScoresOld<'_>]],
    model: &[LLPhyFeatureOld],
    arena: &'a Bump,
) -> Vec<'a, G2WScoresOld<'a>> {
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
        scores.push(G2WScoresOld {
            feature_sum,
            subfeatures: leak_vec(subfeatures),
        })
    }
    scores
}
fn load_reference_g2w_scores(arena: &Bump) -> Vec<'_, G2WScoresOld<'_>> {
    todo!()
}
fn load_post_processor(score_type: ScoreType, arena: &Bump) -> PostProcessor<'_> {
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
fn post_process(post_processor: PostProcessor<'_>, output_scores: &mut [G2WScoresOld<'_>]) {
    for output_row in output_scores {
        post_processor.transform(output_row);
    }
}
fn load_feature_names(arena: &Bump) -> Vec<'_, &str> {
    todo!()
}
fn write_output(
    path: Option<PathBuf>,
    output_scores: &[G2WScoresOld<'_>],
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
        bumpalo::format!(in arena, "{}-feature sum", feature_names.len() / 2);
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
