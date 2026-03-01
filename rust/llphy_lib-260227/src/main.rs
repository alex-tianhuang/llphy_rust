use anyhow::Error;
use clap::Parser;
use llphyscore::{Args, bin_main};

/// See [`llphyscore::bin_main`].
fn main() -> Result<(), Error> {
    let args = Args::try_parse()?;
    bin_main(args)
}
