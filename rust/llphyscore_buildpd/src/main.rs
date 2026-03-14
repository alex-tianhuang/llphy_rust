//! Build script helper that sets up the package data
//! (i.e. model thresholds/signs and human reference data)
//! that these LLPhyScore packages are expecting.
use crate::{grid_decoders::build_griddecoders, grid_scorers::build_gridscorers, reference_data::build_referencedata};
use anyhow::Error;
use clap::Parser;
use std::path::PathBuf;
mod grid_decoders;
mod grid_scorers;
mod reference_data;

/// Program name stands for llphyscore-build(p)ackage(d)ata.
/// 
/// Build script helper that sets up the package data
/// (i.e. model thresholds/signs and human reference data)
/// that these LLPhyScore packages are expecting.
#[derive(Parser)]
#[clap(verbatim_doc_comment, name = "llphyscore-buildpd")]
struct Args {
    /// The input `ScoreDB` directory, from Cai's (@haocai1992)
    /// LLPhyScore repository (in the standalone package, not the repository root).
    input_db: PathBuf,
    /// The output `pkg_data` directory that my LLPhyScore
    /// implementations are expecting.
    /// 
    /// If this directory is already filled (i.e. because the script was
    /// run previously), then this program will simply do a check that
    /// loading from Cai's files yields the same as what is already there.
    output_db: PathBuf,
}
/// See [`Args`] docs.
fn main() -> Result<(), Error> {
    let Args {
        input_db: scoredb_root,
        output_db: pkg_data_root,
    } = Args::parse();
    build_gridscorers(&scoredb_root, &pkg_data_root)?;
    build_griddecoders(&scoredb_root, &pkg_data_root)?;
    build_referencedata(&scoredb_root, &pkg_data_root)?;
    Ok(())
}
