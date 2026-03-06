//! Module defining [`load_grid_scorer`].
use crate::{
    datatypes::{AAMap, Aminoacid},
    featurizer::{
        GridScorer,
        grid_scorer::{
            AvgSdevDB, PairFreqDB, XmerSize, ZGridDB, ZGridDBEntry, ZGridSubtable, xmer_sizes,
        },
    },
    load_pkg_data::legacy::read_archive_file,
};
use anyhow::{Context, Error};
use bumpalo::Bump;
use serde::Deserialize;
use serde_pickle::DeOptions;
use std::{collections::BTreeMap, mem::MaybeUninit, path::Path, ptr::addr_of_mut};

/// A list of tag tuples that are associated to each pair.
///
/// These are the old tag names that Cai (@haocai1992)
/// used to label subfeatures in his `PCON2.FREQS.wBOOTDEV` files.
const PAIR_NAMES_AND_PAIR_FREQ_TAGS: &'static [(&'static str, [&'static str; 2])] = &[
    ("S2.SUMPI", ["srpipi", "lrpipi"]),
    ("S3.WATER.V2", ["Water", "Carbon"]),
    // "ssL" also in original but seems not intentional
    // "ssH" in database but not in here
    ("S4.SSPRED", ["ssH", "ssE"]),
    ("S5.DISO", ["disL", "disS"]),
    ("S6.CHARGE.V2", ["srELEC", "lrELEC"]),
    ("S7.ELECHB.V2", ["sr_hb", "lr_hb"]),
    ("S8.CationPi.V2", ["srCATPI", "lrCATPI"]),
    ("S9.LARKS.V2", ["larkSIM", "larkFAR"]),
];
/// Helper to add line context to errors.
macro_rules! fail_on_line {
    ($reason:expr, $line:ident) => {
        $reason.context(format!("failed to parse line: {}", $line))
    };
}
/// Helper to fail with a literal error message.
macro_rules! fail_on_line_with_msg {
    ($msg:literal, $line:ident) => {
        fail_on_line!(Error::msg($msg), $line)
    };
}
/// Helper to add line context to errors as a `map_err` closure.
macro_rules! map_err_on_line {
    ($line:ident) => {
        |e| fail_on_line!(Error::from(e), $line)
    };
}

/// Load a [`GridScorer`] from
/// Cai's (@haocai1992) old data files.
/// 
/// Allocates into the given memory arena because
/// `GridScorer`s are 1MB structs and I keep getting
/// segfaults which I suspect are stack overflows.
pub fn load_grid_scorer<'a>(pair_name: &str, arena: &'a Bump) -> Result<&'a GridScorer<'a>, Error> {
    // Since `GridScorer` always borrows from the arena and is
    // therefore able to have no drop glue, if at any point the
    // initialization of this grid scorer fails it can be forgotten
    // because the memory arena will deallocate all the relevant space.
    let grid_scorer = arena.alloc_with(<MaybeUninit<GridScorer>>::uninit);
    let target = grid_scorer.as_mut_ptr();
    let subdir = Path::new(pair_name);
    let filepath = subdir.join("PCON2.FREQS.wBOOTDEV");
    let pair_freqs =
        unsafe { &mut *addr_of_mut!((*target).pair_freqs).cast::<MaybeUninit<PairFreqDB>>() };
    load_pair_freq_db_into(pair_freqs, &filepath, pair_name)
        .with_context(|| format!("failed to load frequency pair DB @ {}", filepath.display()))?;
    let xmer_dir = subdir.join("STEP4_AVGnSDEVS");
    let avg_sdevs =
        unsafe { &mut *addr_of_mut!((*target).avg_sdevs).cast::<MaybeUninit<AvgSdevDB>>() };
    load_avg_sdev_db_into(avg_sdevs, &xmer_dir)?;
    let filepath = subdir.join("STEP6_PICKLES").join("SC_GRIDS.pickle4");
    let z_grid = unsafe { &mut *addr_of_mut!((*target).z_grid).cast::<MaybeUninit<ZGridDB>>() };
    load_z_grid_db_into(z_grid, &filepath, arena)
        .with_context(|| format!("failed to load Z grid DB @ {}", filepath.display()))?;
    // SAFETY: just initialized all fields above
    Ok(unsafe { grid_scorer.assume_init_ref() })
}
/// Given space for a [`PairFreqDB`],
/// write data from the given file into it.
///
/// It is expected that the filepath contains plaintext
/// column data.
///
/// See also [`PAIR_NAMES_AND_PAIR_FREQ_TAGS`].
fn load_pair_freq_db_into(
    this: &mut MaybeUninit<PairFreqDB>,
    filepath: &Path,
    pair_name: &str,
) -> Result<(), Error> {
    let Some((_, tags)) = PAIR_NAMES_AND_PAIR_FREQ_TAGS
        .iter()
        .find(|(known_name, _)| *known_name == pair_name)
    else {
        return Err(Error::msg(format!("unrecognized pair name: {}", pair_name)));
    };
    let bytes = read_archive_file(filepath)?;
    let s = String::from_utf8(bytes)?;
    let pair_freq_db = PairFreqDB::init_with_nans(this);
    for line in s.lines() {
        let mut parts = line.split_ascii_whitespace();
        let pair_key = parts
            .next()
            .ok_or_else(|| fail_on_line_with_msg!("expected non-empty line", line))?;
        let [aa_x, b'_', ref middle @ .., b'_', aa_y] = *pair_key.as_bytes() else {
            return Err(fail_on_line_with_msg!(
                "expected `pair_key` string of the form `{aa_x}_{N}_{aa_y}` where `N` is a number",
                line
            ));
        };
        // SAFETY: `middle` comes from a `str` with no unpaired surrogates
        //         and is flanked by two ASCII characters so therefore
        //         cannot contain unpaired surrogates.
        let separation = unsafe { str::from_utf8_unchecked(middle) };
        let aa_x = Aminoacid::try_from(aa_x).map_err(map_err_on_line!(line))?;
        let aa_y = Aminoacid::try_from(aa_y).map_err(map_err_on_line!(line))?;
        let separation = separation
            .parse::<usize>()
            .map_err(map_err_on_line!(line))?;
        while let Some(ptype) = parts.next() {
            let is_x = ptype.starts_with("X_");
            if !(is_x || ptype.starts_with("Y_")) {
                return Err(fail_on_line_with_msg!(
                    "expected `ptype` column to start with `X_` or `Y_`",
                    line
                ));
            }
            let freq = parts
                .next()
                .ok_or_else(|| fail_on_line_with_msg!("expected frequency (float)", line))?
                .parse::<f64>()
                .map_err(map_err_on_line!(line))?;
            let sdev = parts
                .next()
                .ok_or_else(|| fail_on_line_with_msg!("expected frequency (float)", line))?
                .parse::<f64>()
                .map_err(map_err_on_line!(line))?;
            let weight = freq / sdev;
            let total = 1.0 / sdev;
            debug_assert!(weight.is_finite());
            debug_assert!(total.is_finite());
            let (target, tag) = if is_x {
                (
                    &mut pair_freq_db[aa_x][separation].c_terminal_mapping[aa_y],
                    ptype.strip_prefix("X_").unwrap(),
                )
            } else {
                (
                    &mut pair_freq_db[aa_y][separation].n_terminal_mapping[aa_x],
                    ptype.strip_prefix("Y_").unwrap(),
                )
            };
            // Silently ignore tags that are not recognized.
            // If the tag names are wrong, the debug assertion
            // `pair_freq_db.is_nan_free()` will fail.
            if tags[0] == tag {
                target.set_a(weight, total);
            }
            if tags[1] == tag {
                target.set_b(weight, total);
            }
        }
    }
    debug_assert!(pair_freq_db.is_nan_free());
    Ok(())
}
/// Given space for a [`ZGridDB`],
/// write data from the given file into it.
///
/// Will use the given memory arena
/// to allocate dynamic amounts of memory.
///
/// It is expected that the given filepath contains
/// a pickle file that can be unpickled using
/// builtin python types.
fn load_z_grid_db_into<'a>(
    this: &mut MaybeUninit<ZGridDB<'a>>,
    filepath: &Path,
    arena: &'a Bump,
) -> Result<(), Error> {
    let bytes = read_archive_file(filepath)?;
    let bytes = &*bytes;
    /// Helper struct for getting sorted floats out of this pickle file.
    #[derive(Deserialize, PartialEq, PartialOrd)]
    struct SortableFloat(f64);
    impl Eq for SortableFloat {}
    impl Ord for SortableFloat {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            self.0.total_cmp(&other.0)
        }
    }
    let mut unpickled_value = serde_pickle::from_slice::<
        AAMap<Option<BTreeMap<usize, BTreeMap<SortableFloat, BTreeMap<SortableFloat, [i64; 3]>>>>>,
    >(bytes, DeOptions::default())?;
    let z_grid_db = ZGridDB::as_uninit_inner(this);
    for ((aa, target), slot) in z_grid_db.iter_mut().zip(unpickled_value.values_mut()) {
        let entryl1 = slot.take().ok_or_else(|| {
            Error::msg(format!(
                "expected dense aminoacid map, but missing entry for {}",
                aa
            ))
        })?;
        #[cfg(debug_assertions)]
        {
            use crate::datatypes::MAX_XMER;
            assert_eq!(entryl1.len(), MAX_XMER);
            let (&first, _) = entryl1.first_key_value().unwrap();
            assert_eq!(first, 1);
            let (&last, _) = entryl1.last_key_value().unwrap();
            assert_eq!(last, MAX_XMER);
        }
        for (xmer, entryl2) in entryl1 {
            let xmer = XmerSize::new(xmer).unwrap();
            let min_a = entryl2.first_key_value().unwrap().0.0;
            let max_a = entryl2.last_key_value().unwrap().0.0;
            let column_len = (max_a - min_a) * 2.0 + 1.0;
            debug_assert_eq!(column_len, column_len.round());
            let column_len = column_len as usize;
            let min_b = entryl2
                .iter()
                .map(|t| t.1.first_key_value().unwrap().0.0)
                .min_by(f64::total_cmp)
                .unwrap();
            let max_b = entryl2
                .iter()
                .map(|t| t.1.last_key_value().unwrap().0.0)
                .max_by(f64::total_cmp)
                .unwrap();
            let row_len = (max_b - min_b) * 2.0 + 1.0;
            debug_assert_eq!(row_len, row_len.round());
            let row_len = row_len as usize;
            let dbl_z_offsets = [min_a * 2.0, min_b * 2.0];
            let data =
                arena.alloc_slice_fill_copy(row_len * column_len, ZGridDBEntry::unoccupied());
            for (key_a, entryl3) in entryl2.iter() {
                let dbl_z_a = (key_a.0 - min_a) * 2.0;
                debug_assert_eq!(dbl_z_a, dbl_z_a.round());
                let idx_a = dbl_z_a as usize;
                let row = &mut data[idx_a * row_len..(idx_a + 1) * row_len];
                for (key_b, &[weight_total, weight_a, weight_b]) in entryl3.iter() {
                    let dbl_z_b = (key_b.0 - min_b) * 2.0;
                    debug_assert_eq!(dbl_z_b, dbl_z_b.round());
                    let idx_b = dbl_z_b as usize;
                    row[idx_b] = ZGridDBEntry::new_occupied(weight_total, weight_a, weight_b);
                }
            }
            target[xmer].write(ZGridSubtable::new(dbl_z_offsets, row_len, data));
        }
    }
    Ok(())
}
/// Given space for an [`AvgSdevDB`], initialize it
/// using data from the given `xmer_dir`.
///
/// It is expected that the subdirectory contains files
/// of the form `PCON.xmer{N}` for all integers `N` in
/// `1..=MAX_XMER`.
fn load_avg_sdev_db_into<'a>(
    this: &mut MaybeUninit<AvgSdevDB>,
    xmer_dir: &Path,
) -> Result<(), Error> {
    let avg_sdev_db = AvgSdevDB::init_with_nans(this);
    for xmer in xmer_sizes() {
        let filepath = xmer_dir.join(format!("PCON2.xmer{}", xmer.get()));
        load_one_xmer_avg_sdev_into(avg_sdev_db, &filepath, xmer).with_context(|| {
            format!(
                "failed to load xmer avg/stdev statistics @ {}",
                filepath.display()
            )
        })?;
    }
    debug_assert!(avg_sdev_db.is_nan_free());
    Ok(())
}
/// Load the data for one `PCON2.xmer{}`
/// file into the correct slot.
///
/// Helper function for [`load_avg_sdev_db_into`].
fn load_one_xmer_avg_sdev_into(
    this: &mut AvgSdevDB,
    filepath: &Path,
    xmer: XmerSize,
) -> Result<(), Error> {
    let bytes = read_archive_file(&filepath)?;
    let s = String::from_utf8(bytes)?;
    for line in s.lines() {
        let [aa, rest @ ..] = line.trim_ascii_start().as_bytes() else {
            return Err(fail_on_line_with_msg!("expected non-empty line", line));
        };
        let aa = Aminoacid::try_from(*aa).map_err(map_err_on_line!(line))?;
        // SAFETY: `middle` comes from a `str` with no unpaired surrogates
        //         and is flanked by two ASCII characters so therefore
        //         cannot contain unpaired surrogates.
        let rest = unsafe { str::from_utf8_unchecked(rest) };

        let mut parts = rest.trim_ascii_start().split_ascii_whitespace();
        let entry = &mut this[aa][xmer];
        /// Shorthand for parsing fields associated with feature `A` and then feature `B`.
        macro_rules! parse_ab_fields {
            ($entry:ident, $parts:ident, $($method:ident),*) => {{
                $(let _ = $parts
                    .next()
                    .ok_or_else(|| fail_on_line_with_msg!("expected a string", line))?;
                let avg = $parts
                    .next()
                    .ok_or_else(|| fail_on_line_with_msg!("expected an average (float)", line))?
                    .parse::<f64>()
                    .map_err(map_err_on_line!(line))?;
                let std = $parts
                    .next()
                    .ok_or_else(|| {
                        fail_on_line_with_msg!("expected an standard deviation (float)", line)
                    })?
                    .parse::<f64>()
                    .map_err(map_err_on_line!(line))?;
                let invstd = 1.0 / std;
                debug_assert!(invstd.is_finite());
                entry.$method(avg, invstd);)*
            }};
        }
        parse_ab_fields!(entry, parts, set_a, set_b);
    }
    Ok(())
}
