use std::{array, collections::BTreeMap, path::Path, str::FromStr};
use anyhow::{Context, Error};
use bumpalo::{Bump, collections::Vec};
use serde_pickle::DeOptions;

use crate::{
    datatypes::{AAMap, Aminoacid, GridScorer, LineKey, MAX_XMER},
    leak_vec,
    load_legacy::read_archive_file,
};

pub fn load_legacy_pdb_statistics(
    arena: &Bump,
) -> Result<
    (
        &'static [(&'static str, [&'static str; 2])],
        &[GridScorer<'_>],
    ),
    Error,
> {
    const OLD_GRID_NAMES: &'static [(&'static str, [&'static str; 2])] = &[
        ("S2.SUMPI", ["srpipi", "lrpipi"]),
        ("S3.WATER.V2", ["Water", "Carbon"]),
        // "ssL" also in original but seems not intentional
        // "ssH" in database but not in here... yikes...
        ("S4.SSPRED", ["ssH", "ssE"]),
        ("S5.DISO", ["disL", "disS"]),
        ("S6.CHARGE.V2", ["srELEC", "lrELEC"]),
        ("S7.ELECHB.V2", ["sr_hb", "lr_hb"]),
        ("S8.CationPi.V2", ["srCATPI", "lrCATPI"]),
        ("S9.LARKS.V2", ["larkSIM", "larkFAR"]),
    ];
    const NEW_GRID_NAMES: &'static [(&'static str, [&'static str; 2])] = &[
        ("S2.SUMPI", ["pi-pi (short-range)", "pi-pi (long-range)"]),
        ("S3.WATER.V2", ["protein-water", "protein-carbon"]),
        (
            "S4.SSPRED",
            ["sec. structure (helices)", "sec. structure (strands)"],
        ),
        ("S5.DISO", ["disorder (long)", "disorder (short)"]),
        (
            "S6.CHARGE.V2",
            ["electrostatic (short-range)", "electrostatic (long-range)"],
        ),
        (
            "S7.ELECHB.V2",
            ["hydrogen bond (short-range)", "hydrogen bond (long-range)"],
        ),
        (
            "S8.CationPi.V2",
            ["cation-pi (short-range)", "cation-pi (long-range)"],
        ),
        (
            "S9.LARKS.V2",
            ["K-Beta similarity", "K-Beta non-similarity"],
        ),
    ];
    
    #[cfg(debug_assertions)]
    {
        assert_eq!(OLD_GRID_NAMES.len(), NEW_GRID_NAMES.len());
        for ((old_grid_name, _), (new_grid_name, _)) in OLD_GRID_NAMES.iter().zip(NEW_GRID_NAMES) {
            assert_eq!(old_grid_name, new_grid_name)
        }
    }
    let mut pdb_statistics = Vec::with_capacity_in(OLD_GRID_NAMES.len(), arena);
    let mut warner = Warner::new(10);
    for (grid_name, [tag_sr, tag_lr]) in OLD_GRID_NAMES.iter().copied() {
        let grid_scorer = load_grid_scorer(Path::new(grid_name), tag_sr, tag_lr, arena, &mut warner)?;
        pdb_statistics.push(grid_scorer);
    }
    Ok((NEW_GRID_NAMES, leak_vec(pdb_statistics)))
}
fn load_grid_scorer<'a>(
    subdir: &Path,
    tag_sr: &'static str,
    tag_lr: &'static str,
    arena: &'a Bump,
    warner: &mut Warner,
) -> Result<GridScorer<'a>, Error> {
    let filepath = subdir.join("PCON2.FREQS.wBOOTDEV");
    let (pair_freq_db_x, pair_freq_db_y) =
        load_pair_freq_db(&filepath, tag_sr, tag_lr, arena, warner).with_context(|| {
            format!("failed to load frequency pair DB @ {}", filepath.display())
        })?;
    let filepath = subdir.join("STEP6_PICKLES").join("SC_GRIDS.pickle4");
    let z_grid_db = load_z_grid_db(&filepath, arena)
        .with_context(|| format!("failed to load Z grid DB @ {}", filepath.display()))?;
    let xmer_dir = subdir.join("STEP4_AVGnSDEVS");
    let avg_sdev_db = load_avg_sdev_db(&xmer_dir, arena)?;
    Ok(GridScorer {
        pair_freq_db_x,
        pair_freq_db_y,
        avg_sdev_db,
        z_grid_db,
    })
}

struct Warner {
    limit: usize,
    warn_count: usize,
}
impl Warner {
    pub fn new(limit: usize) -> Self {
        Self {
            limit,
            warn_count: 0,
        }
    }
    pub fn warn(&mut self, msg: std::fmt::Arguments<'_>) {
        if self.warn_count < self.limit {
            eprintln!("{}", msg);
        }
        self.warn_count += 1;
    }
}
impl Drop for Warner {
    fn drop(&mut self) {
        if self.warn_count > self.limit {
            eprintln!(
                "...\n{} more messages suppressed",
                self.warn_count - self.limit
            )
        }
    }
}
type PairFreqDBInner = AAMap<AAMap<[[f64; 2]; 2]>>;
fn load_pair_freq_db<'a>(
    filepath: &Path,
    tag_sr: &'static str,
    tag_lr: &'static str,
    arena: &'a Bump,
    warner: &mut Warner,
) -> Result<(&'a [PairFreqDBInner], &'a [PairFreqDBInner]), Error> {
    fn initialize_pair_freq_db(arena: &Bump) -> &mut [PairFreqDBInner] {
        debug_assert_eq!(align_of::<PairFreqDBInner>(), align_of::<f64>());
        let cap_in_bytes = (MAX_XMER + 1) * size_of::<PairFreqDBInner>();
        let cap_in_floats = cap_in_bytes / size_of::<f64>();
        let buf = arena.alloc_slice_fill_copy(cap_in_floats, f64::NAN);
        let pair_freq_db: &mut [PairFreqDBInner];
        unsafe {
            #[allow(unused)]
            let (start, mid, end) = buf.align_to_mut::<PairFreqDBInner>();
            debug_assert_eq!(start.len(), 0);
            debug_assert_eq!(end.len(), 0);
            pair_freq_db = mid;
        }
        #[cfg(debug_assertions)]
        check_initialized(pair_freq_db);
        pair_freq_db
    }
    #[cfg(debug_assertions)]
    fn check_initialized(db: &[PairFreqDBInner]) {
        for entry1 in db {
            for entry2 in entry1.values() {
                for entry3 in entry2.values() {
                    for slot in entry3.iter().flatten() {
                        assert!(slot.is_nan())
                    }
                }
            }
        }
    }
    #[cfg(debug_assertions)]
    fn check_filled(db: &[PairFreqDBInner]) {
        for entry1 in db {
            for entry2 in entry1.values() {
                for entry3 in entry2.values() {
                    for slot in entry3.iter().flatten() {
                        assert!(!slot.is_nan())
                    }
                }
            }
        }
    }
    let pair_freq_db_x = initialize_pair_freq_db(arena);
    let pair_freq_db_y = initialize_pair_freq_db(arena);
    let bytes = read_archive_file(filepath)?;
    let mut bytes = &*bytes;
    let mut proc_line = |line: &[u8]| {
        let mut parts = line
            .split(u8::is_ascii_whitespace)
            .filter(|b| !b.is_empty());
        let fail = || {
            Error::msg(format!(
                "failed to parse line: {}",
                String::from_utf8_lossy(line)
            ))
        };
        let pair_key = parts.next().ok_or_else(fail)?;
        let [aa_x, b'_', ref separation @ .., b'_', aa_y] = *pair_key else {
            return Err(fail());
        };
        let aa_x = Aminoacid::try_from(aa_x).map_err(|_| fail())?;
        let aa_y = Aminoacid::try_from(aa_y).map_err(|_| fail())?;
        let separation = parse_bytes::<usize>(separation).ok_or_else(fail)?;
        while let Some(ptype) = parts.next() {
            let is_x = ptype.starts_with(b"X_");
            if !(is_x || ptype.starts_with(b"Y_")) {
                return Err(fail());
            }
            let freq = parts.next().and_then(parse_bytes::<f64>).ok_or_else(fail)?;
            let sdev = parts.next().and_then(parse_bytes::<f64>).ok_or_else(fail)?;
            let (target, first_aa, second_aa, tag) = if is_x {
                (
                    &mut *pair_freq_db_x,
                    aa_x,
                    aa_y,
                    ptype.strip_prefix(b"X_").unwrap(),
                )
            } else {
                (
                    &mut *pair_freq_db_y,
                    aa_y,
                    aa_x,
                    ptype.strip_prefix(b"Y_").unwrap(),
                )
            };
            if ![tag_sr.as_bytes(), tag_lr.as_bytes()].contains(&tag) {
                warner.warn(format_args!(
                    "ignoring line with unrecognized tag = {} (available tags = {:?})",
                    String::from_utf8_lossy(tag),
                    [tag_sr, tag_lr]
                ));
                continue;
            }
            let data_block = &mut target[separation][first_aa][second_aa];
            let is_sr = tag == tag_sr.as_bytes();
            if is_sr {
                data_block[0] = [freq, sdev]
            } else {
                data_block[1] = [freq, sdev]
            }
        }
        Ok(())
    };
    while let Some(idx) = next_line_idx(bytes) {
        let (line, rest) = bytes.split_at(idx);
        bytes = &rest[1..];
        proc_line(line)?;
    }
    if !bytes.is_empty() {
        proc_line(bytes)?;
    }
    #[cfg(debug_assertions)]
    {
        check_filled(pair_freq_db_x);
        check_filled(pair_freq_db_y);
    }
    Ok((pair_freq_db_x, pair_freq_db_y))
}
type ZGridDB<'a> = AAMap<&'a [&'a [(LineKey, &'a [(LineKey, [isize; 3])])]]>;
fn load_z_grid_db<'a>(filepath: &Path, arena: &'a Bump) -> Result<ZGridDB<'a>, Error> {
    let bytes = read_archive_file(filepath)?;
    let bytes = &*bytes;
    let mut unpickled_value = serde_pickle::from_slice::<
        AAMap<Option<BTreeMap<usize, BTreeMap<LineKey, BTreeMap<LineKey, [isize; 3]>>>>>,
    >(bytes, DeOptions::default())?;
    let mut z_grid_db: ZGridDB = AAMap([const { &[] }; 20]);
    for ((aa, target), slot) in z_grid_db.iter_mut().zip(unpickled_value.values_mut()) {
        let entryl1 = slot.take().ok_or_else(|| {
            Error::msg(format!(
                "expected dense aminoacid map, but missing entry for {}",
                aa
            ))
        })?;
        #[cfg(debug_assertions)]
        {
            assert_eq!(entryl1.len(), MAX_XMER);
            let (&first, _) = entryl1.first_key_value().unwrap();
            assert_eq!(first, 1);
            let (&last, _) = entryl1.last_key_value().unwrap();
            assert_eq!(last, MAX_XMER);
        }
        *target = &*arena.alloc_slice_fill_iter(entryl1.into_values().map(|entryl2| {
            &*arena.alloc_slice_fill_iter(entryl2.into_iter().map(|(sr_gridpoint, entryl3)| {
                (sr_gridpoint, &*arena.alloc_slice_fill_iter(entryl3))
            }))
        }));
    }
    Ok(z_grid_db)
}
type AvgSdevDB<'a> = AAMap<&'a [[[f64; 2]; 2]]>;
fn load_avg_sdev_db<'a>(xmer_dir: &Path, arena: &'a Bump) -> Result<AvgSdevDB<'a>, Error> {
    let default_block = [[f64::NAN, f64::NAN], [f64::NAN, f64::NAN]];
    let mut avg_sdev_db = AAMap(array::from_fn(|_| {
        arena.alloc_slice_fill_copy(MAX_XMER, default_block)
    }));
    for idx in 0..MAX_XMER {
        let xmer = idx + 1;
        let filepath = xmer_dir.join(format!("PCON2.xmer{}", xmer));
        let bytes = read_archive_file(&filepath)?;
        let mut bytes = &*bytes;
        let mut proc_line = |line: &[u8]| {
            let fail = || {
                Error::msg(format!(
                    "failed to parse line: {}",
                    String::from_utf8_lossy(line)
                ))
                .context(format!(
                    "failed to load xmer avg/stdev statistics @ {}",
                    filepath.display()
                ))
            };
            let [aa, rest @ ..] = trim_start_whitespace(line) else {
                return Err(fail());
            };
            let aa = Aminoacid::try_from(*aa).map_err(|_| fail())?;
            let rest = trim_start_whitespace(rest);
            let mut parts = rest
                .split(u8::is_ascii_whitespace)
                .filter(|b| !b.is_empty());
            let data_block = &mut avg_sdev_db[aa][idx];
            for sub_block in data_block {
                let _ = parts.next().ok_or_else(fail)?;
                let avg = parts.next().and_then(parse_bytes::<f64>).ok_or_else(fail)?;
                let sdev = parts.next().and_then(parse_bytes::<f64>).ok_or_else(fail)?;
                *sub_block = [avg, sdev]
            }
            Ok(())
        };
        while let Some(idx) = next_line_idx(bytes) {
            let (line, rest) = bytes.split_at(idx);
            bytes = &rest[1..];
            proc_line(line)?;
        }
        if !bytes.is_empty() {
            proc_line(bytes)?;
        }
    }
    Ok(AAMap(avg_sdev_db.0.map(|v| v as &[_])))
}
fn next_line_idx(bytes: &[u8]) -> Option<usize> {
    let address = bytes.iter().find(|b| **b == b'\n')?;
    Some(address as *const u8 as usize - bytes.as_ptr() as usize)
}
fn parse_bytes<T: FromStr>(b: &[u8]) -> Option<T> {
    T::from_str(str::from_utf8(b).ok()?).ok()
}
fn trim_start_whitespace(bytes: &[u8]) -> &[u8] {
    let mut truncate_idx = bytes.len();
    for b in bytes {
        if !b.is_ascii_whitespace() {
            truncate_idx = b as *const u8 as usize - bytes.as_ptr() as usize;
            break;
        }
    }
    &bytes[truncate_idx..]
}
