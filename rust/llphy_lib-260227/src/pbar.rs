use std::fmt::Write;
use indicatif::{ProgressBar, ProgressState, ProgressStyle};

/// Pbar in my preferred style.
pub fn pbar(n: u64) -> ProgressBar {
    ProgressBar::new(n).with_style(
        ProgressStyle::with_template(
            "{percent}%[{wide_bar}] {pos}/{len} [{elapsed}<{eta}, {rate}]",
        )
        .unwrap()
        .with_key("rate", |state: &ProgressState, w: &mut dyn Write| {
            write!(w, "{:.2}it/s", state.per_sec()).unwrap();
        }),
    )
}
