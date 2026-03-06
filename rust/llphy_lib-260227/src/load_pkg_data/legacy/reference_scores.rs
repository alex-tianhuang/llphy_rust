//! Module defining [`load_reference_scores`].

use std::path::{Path, PathBuf};

use anyhow::{Context, Error};
use bumpalo::{Bump, collections::Vec};

use crate::{datatypes::{ModelTrainingBase, ReferenceFeatureMatrix}, leak_vec, load_pkg_data::legacy::read_archive_file};

/// Load `human` sequence reference scores computed from
/// the model trained on the given `model_train_base`.
/// 
/// Use the memory arena to allocate space for these scores.
pub fn load_reference_scores(model_train_base: ModelTrainingBase, arena: &Bump) -> Result<ReferenceFeatureMatrix<'_>, Error> {
    let filepath = PathBuf::from(format!(
        "human_g2w_scores_using_{}weight.csv",
        model_train_base
    ));
    load_reference_scores_from_filepath(&filepath, arena).with_context(|| format!(
        "failed to parse reference g2w scores in archive @ {}",
        filepath.display()
    ))
}
fn load_reference_scores_from_filepath<'a>(filepath: &Path, arena: &'a Bump) -> Result<ReferenceFeatureMatrix<'a>, Error> {
    let bytes = read_archive_file(&filepath)?;
    let mut reader = csv::Reader::from_reader(&*bytes);
    let mut record = csv::StringRecord::new();
    let num_features: usize;
    let num_ref_seqs: usize;
    let mut feature_names: Vec<'_, &str>;
    {
        let headers = reader.headers()?;
        num_features = headers.len().checked_sub(1).ok_or_else(|| Error::msg("expected a `feature_sum` column"))?;
        num_ref_seqs = bytes.iter().filter(|b| **b == b'\n').count() - 1;
        let headers = headers.into_iter();
        feature_names = Vec::with_capacity_in(num_features, arena);
        feature_names.extend(headers.take(num_features).map(|s| &*arena.alloc_str(s)));
        debug_assert_eq!(feature_names.len(), num_features);
    }
    let feature_names = &*leak_vec(feature_names);
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
        let sumfeature_value = record_iter.next().ok_or_else(|| Error::msg("expected a `feature_sum` column"))?.parse::<i64>()?;
        debug_assert_ne!(sumfeature_value, UNINIT_SENTINEL);
        data[num_features * num_ref_seqs + cursor] = sumfeature_value;
        cursor += 1;
    }
    debug_assert_eq!(cursor, num_ref_seqs);
    Ok(ReferenceFeatureMatrix { feature_names, num_ref_seqs, data })
}