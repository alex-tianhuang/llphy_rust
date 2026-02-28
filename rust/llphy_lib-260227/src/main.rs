use std::{
    fs::File,
    io::{BufRead, Cursor, Read},
    path::PathBuf,
};
mod io;
use bumpalo::{Bump, collections::Vec};
use clap::Parser;
mod datatypes;
mod fasta;

#[derive(Parser)]
struct Args {
    #[arg(short, long)]
    input_file: PathBuf,
    #[arg(short, long)]
    output_file: PathBuf,
}
fn main() -> Result<(), anyhow::Error> {
    let args = Args::try_parse()?;
    Ok(())
}

fn run_fasta_scorer(args: Args) -> Result<(), anyhow::Error> {
    let Args {
        input_file,
        output_file,
    } = args;
    Ok(())
}
