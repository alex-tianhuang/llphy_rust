//! Module defining [`build_gridscorers`].

use crate::grid_scorers::load::load_grid_scorer;
use anyhow::Error;
use borsh::BorshSerialize;
use bumpalo::Bump;
use llphyscore_core::datatypes::{GridScorer, PAIR_NAMES_AND_FEATURE_NAMES};
use std::{fs::File, io::Read, path::Path};
mod load;

/// Load all [`GridScorer`]s from Cai's (@haocai1992) old data files,
/// serialize them, verify roundtrip integrity, and then write them
/// to an appropriate output directory.
///
/// The old data files are read from a database at `scoredb_root` and the
/// new `gridscorer.bin` files are dumped to a database at `pkg_data_root`.
pub fn build_gridscorers(scoredb_root: &Path, pkg_data_root: &Path) -> Result<(), Error> {
    let mut arena = Bump::new();
    let mut bytes = std::vec::Vec::new();
    for (pair_name, _) in PAIR_NAMES_AND_FEATURE_NAMES {
        arena.reset();
        let grid_scorer = load_grid_scorer(scoredb_root, pair_name, &arena)?;
        let outdir = pkg_data_root.join("feature_pairs").join(pair_name);
        let outpath = outdir.join(format!("gridscorer.bin"));
        bytes.clear();
        if outpath.exists() {
            let mut file = File::open(&outpath)?;
            file.read_to_end(&mut bytes)?;
            let round_trip = GridScorer::deserialize(&mut &*bytes, &arena)?;
            if round_trip != grid_scorer {
                let inpath = scoredb_root.join(pair_name);
                return Err(Error::msg(format!(
                    "pre-existing file @ {} is not the same grid scorer as loaded from {}",
                    outpath.display(),
                    inpath.display()
                )));
            }
        } else {
            grid_scorer.serialize(&mut bytes)?;
            let round_trip = GridScorer::deserialize(&mut &*bytes, &arena)?;
            assert!(round_trip == grid_scorer);
            std::fs::create_dir_all(&outdir)?;
            let mut file = File::options()
                .create(true)
                .truncate(true)
                .write(true)
                .open(&outpath)?;
            grid_scorer.serialize(&mut file)?;
        }
    }
    Ok(())
}
