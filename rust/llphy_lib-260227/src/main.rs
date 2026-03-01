use anyhow::Error;
use clap::Parser;
use llphyscore::{Args, run_fasta_scorer};

fn main() -> Result<(), Error> {
    let args = Args::try_parse()?;
    run_fasta_scorer(args)
}
