//! Module defining [`build_referencedata`].
use anyhow::Error;
use borsh::BorshSerialize;
use bumpalo::Bump;
use llphyscore_core::datatypes::{ModelTrainingBase, ReferenceFeatureMatrix};
use std::{fs::File, io::Read, path::Path};

/// Load all [`ReferenceFeatureMatrix`]s from Cai's (@haocai1992) old data files,
/// serialize them, verify roundtrip integrity, and then write them
/// to an appropriate output directory.
///
/// The old data files are read from a database at `scoredb_root` and the
/// new `gridscorer.bin` files are dumped to a database at `pkg_data_root`.
pub fn build_referencedata(scoredb_root: &Path, pkg_data_root: &Path) -> Result<(), Error> {
    let mut arena = Bump::new();
    let mut bytes = Vec::new();
    for model_train_base in [
        ModelTrainingBase::Human,
        ModelTrainingBase::HumanPDB,
        ModelTrainingBase::PDB,
    ] {
        arena.reset();
        let inpath = scoredb_root.join(format!(
            "human_g2w_scores_using_{}weight.csv",
            model_train_base
        ));
        let matrix = load_reference_scores(&inpath, &arena)?;
        let outpath = pkg_data_root
            .join("human_reference_data")
            .join(format!("{}.distr.bin", model_train_base));
        bytes.clear();
        if outpath.exists() {
            let mut file = File::open(&outpath)?;
            file.read_to_end(&mut bytes)?;
            let round_trip = ReferenceFeatureMatrix::deserialize(&mut &*bytes, &arena)?;
            if round_trip != matrix {
                return Err(Error::msg(format!(
                    "pre-existing file @ {} is not the same human reference data as loaded from {}",
                    outpath.display(),
                    inpath.display()
                )));
            }
        } else {
            matrix.serialize(&mut bytes)?;
            let round_trip = ReferenceFeatureMatrix::deserialize(&mut &*bytes, &arena)?;
            assert!(round_trip == matrix);
            let mut file = File::options()
                .create(true)
                .truncate(true)
                .write(true)
                .open(&outpath)?;
            matrix.serialize(&mut file)?;
        }
    }
    Ok(())
}
/// Load `human` sequence reference scores computed from
/// the model trained on the given `model_train_base`.
///
/// Use the memory arena to allocate space for these scores.
fn load_reference_scores<'a>(
    filepath: &Path,
    arena: &'a Bump,
) -> Result<ReferenceFeatureMatrix<'a>, Error> {
    let mut bytes = Vec::new();
    File::open(filepath)?.read_to_end(&mut bytes)?;
    let mut reader = csv::Reader::from_reader(&*bytes);
    let mut record = csv::StringRecord::new();
    let num_features: usize;
    let num_ref_seqs: usize;
    let feature_names: &[&str];
    {
        let headers = reader.headers()?;
        num_features = headers
            .len()
            .checked_sub(1)
            .ok_or_else(|| Error::msg("expected a `feature_sum` column"))?;
        num_ref_seqs = bytes.iter().filter(|b| **b == b'\n').count() - 1;
        let headers = headers.into_iter();
        const UNINIT_SENTINEL: &'static str = "";
        let feature_names_uninit = arena.alloc_slice_fill_copy(num_features, UNINIT_SENTINEL);
        for (slot, header) in feature_names_uninit.iter_mut().zip(headers) {
            *slot = &*arena.alloc_str(header);
        }
        #[cfg(debug_assertions)]
        for slot in feature_names_uninit.iter_mut() {
            assert_ne!(*slot, UNINIT_SENTINEL);
        }
        feature_names = feature_names_uninit;
    }
    const UNINIT_SENTINEL: i64 = i64::MIN;
    let data = arena.alloc_slice_fill_copy((num_features + 1) * num_ref_seqs, UNINIT_SENTINEL);
    let mut cursor = 0;
    while reader.read_record(&mut record)? {
        debug_assert_eq!(record.len(), num_features + 1);
        let mut record_iter = record.iter();
        for (i, subfeature_value) in (&mut record_iter).take(num_features).enumerate() {
            let subfeature_value = subfeature_value.parse::<i64>()?;
            debug_assert_ne!(subfeature_value, UNINIT_SENTINEL);
            data[cursor + i * num_ref_seqs] = subfeature_value;
        }
        let sumfeature_value = record_iter
            .next()
            .ok_or_else(|| Error::msg("expected a `feature_sum` column"))?
            .parse::<i64>()?;
        debug_assert_ne!(sumfeature_value, UNINIT_SENTINEL);
        data[num_features * num_ref_seqs + cursor] = sumfeature_value;
        cursor += 1;
    }
    debug_assert_eq!(cursor, num_ref_seqs);
    Ok(ReferenceFeatureMatrix {
        feature_names,
        num_ref_seqs,
        data,
    })
}
