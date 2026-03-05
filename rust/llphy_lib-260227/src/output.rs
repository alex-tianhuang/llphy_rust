//! Module defining [`write_output`].
use std::{io::stdout, path::PathBuf};

use anyhow::Error;
use bumpalo::{
    Bump,
    collections::{String, Vec},
};

use crate::datatypes::{FastaEntry, PostProcessedFeatureMatrix};

/// Output the final biophysical scores (currently [`G2WScores`])
/// associated with the given `sequences` to a file at the given
/// `path` or stdout.
pub fn write_output(
    path: Option<PathBuf>,
    sequences: &[FastaEntry<'_>],
    matrix: PostProcessedFeatureMatrix<'_>,
    arena: &Bump,
) -> Result<(), Error> {
    match path {
        Some(path) => {
            let mut file_writer = csv::Writer::from_path(path)?;
            write_output_virtualized(
                &mut |row| file_writer.write_record(row),
                sequences,
                matrix,
                arena,
            )
        }
        None => {
            let mut stdout_writer = csv::Writer::from_writer(stdout().lock());
            write_output_virtualized(
                &mut |row| stdout_writer.write_record(row),
                sequences,
                matrix,
                arena,
            )
        }
    }
}
/// Helper method for [`write_output`]
/// that doesn't depend on the output file handle.
fn write_output_virtualized(
    writer: &mut dyn FnMut(&[&str]) -> Result<(), csv::Error>,
    sequences: &[FastaEntry<'_>],
    matrix: PostProcessedFeatureMatrix<'_>,
    arena: &Bump,
) -> Result<(), Error> {
    let feature_names = matrix.feature_names();
    let mut column_buffer = Vec::with_capacity_in(feature_names.len() + 2, arena);
    let feature_sum_column_name = bumpalo::format!(in arena, "{}-feature sum", feature_names.len());
    let mut values_buffer = String::new_in(arena);
    column_buffer.push("tag");
    column_buffer.extend_from_slice(feature_names);
    // SAFETY: this field's lifetime does not need to
    //         outlive the current function even though rustc
    //         is convinced that it does need to
    column_buffer.push(unsafe { &*(feature_sum_column_name.as_str() as *const str) });
    writer(&column_buffer)?;
    for (entry, row) in sequences.iter().zip(matrix.rows()) {
        let [tag_slot, feature_slots @ ..] = &mut *column_buffer else {
            unreachable!()
        };
        *tag_slot = entry.header;
        row.format_comma_separated_into(&mut values_buffer);
        for (slot, value_str) in feature_slots.iter_mut().zip(values_buffer.split(",")) {
            // SAFETY: this string is not read from
            //         after the end of the loop
            *slot = unsafe { &*(value_str as *const str) };
        }
        writer(&column_buffer)?;
        values_buffer.clear();
    }
    Ok(())
}
