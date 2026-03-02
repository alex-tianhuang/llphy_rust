use std::{io::Read, path::PathBuf, str::FromStr};

use anyhow::Error;
use bumpalo::Bump;

use crate::{G2WScores, datatypes::ModelTrainingType, load_legacy::read_archive_file};

pub fn load_legacy_reference_g2w_scores<'a>(
    model_train_base: ModelTrainingType,
    arena: &'a Bump,
) -> Result<&'a [G2WScores<'a>], Error> {
    let filename = PathBuf::from(format!(
        "human_g2w_scores_using_{}weight.csv",
        model_train_base
    ));
    let bytes = read_archive_file(&filename)?;
    let mut reader = csv::Reader::from_reader(&*bytes);
    let mut record = csv::StringRecord::new();
    let fail = || {
        Error::msg(format!(
            "failed to parse reference g2w scores in archive @ {}",
            filename.display()
        ))
    };
    let headers = reader.headers().map_err(|_| fail())?;
    if headers.is_empty() {
        return Err(fail());
    }
    let num_features = headers.len() - 1;
    let cap = bytes.iter().filter(|b| **b == b'\n').count() - 1;
    let all_g2w_scores = arena.alloc_slice_fill_with(cap, |_| G2WScores {
        feature_sum: f64::NAN,
        subfeatures: arena.alloc_slice_fill_copy(num_features, f64::NAN),
    });
    let mut cursor = 0;
    while reader.read_record(&mut record).map_err(|_| fail())? {
        let entry = &mut all_g2w_scores[cursor];
        let mut cells = record.iter();
        debug_assert_eq!(record.len(), num_features + 1);
        for (subfeatures_slot, subfeature_value) in entry
            .subfeatures
            .iter_mut()
            .zip((&mut cells).take(num_features))
        {
            *subfeatures_slot = f64::from_str(subfeature_value).map_err(|_| fail())?;
        }
        entry.feature_sum = f64::from_str(cells.next().unwrap()).map_err(|_| fail())?;
        cursor += 1;
    }
    debug_assert_eq!(cap, cursor);
    Ok(all_g2w_scores)
}
