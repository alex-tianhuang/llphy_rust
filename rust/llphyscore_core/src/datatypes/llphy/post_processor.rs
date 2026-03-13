//! Module defining a struct for post-processing
//! raw `llphyscore` values into z-scores or percentiles.
use anyhow::Error;
use bumpalo::Bump;
use clap::ValueEnum;

use crate::datatypes::{FeatureMatrix, PostProcessedFeatureMatrix, ReferenceFeatureMatrix};
/// The type of score to return.
#[derive(Clone, ValueEnum)]
#[value(verbatim_doc_comment)]
pub enum ScoreType {
    /// Report raw `llphyscore` values and features.
    Raw,
    /// Report z-scores in comparison to the human reference IDRome.
    ZScore,
    /// Report values in as percentiles compared to the human reference IDRome.
    Percentile,
}
/// The struct that does the requested post-processing
/// (selected based on [`ScoreType`]).
#[derive(PartialEq)]
pub enum PostProcessor<'a> {
    /// Report raw `llphyscore` values and features.
    Raw,
    /// Report z-scores in comparison to the human reference IDRome.
    ZScore {
        feature_names: &'a [&'a str],
        sumfeature_stats: [f64; 2],
        subfeature_stats: &'a [[f64; 2]],
    },
    /// Report values in as percentiles compared to the human reference IDRome.
    Percentile(ReferenceFeatureMatrix<'a>),
}
impl<'a> PostProcessor<'a> {
    /// Make a new post-processor that turns raw feature values into z-scores.
    pub fn new_zscore(ref_scores: ReferenceFeatureMatrix<'a>, arena: &'a Bump) -> Self {
        let ReferenceFeatureMatrix {
            feature_names,
            num_ref_seqs,
            data,
        } = ref_scores;
        let num_features = feature_names.len();
        debug_assert_eq!((num_features + 1) * num_ref_seqs, data.len());
        let subfeature_stats = arena.alloc_slice_fill_with(num_features, |i| {
            mean_std(
                data[i * num_ref_seqs..(i + 1) * num_ref_seqs]
                    .iter()
                    .copied(),
            )
        });
        let sumfeature_stats = mean_std(data[num_features * num_ref_seqs..].iter().copied());
        PostProcessor::ZScore {
            feature_names,
            sumfeature_stats,
            subfeature_stats,
        }
    }
    /// Make a new post-processor that turns raw
    /// feature values into percentile values.
    pub fn new_percentile(mut ref_scores: ReferenceFeatureMatrix<'a>) -> Self {
        let ReferenceFeatureMatrix {
            feature_names,
            num_ref_seqs,
            ref mut data,
        } = ref_scores;
        let num_features = feature_names.len();
        debug_assert_eq!((num_features + 1) * num_ref_seqs, data.len());
        for feature_data in data.chunks_exact_mut(num_ref_seqs) {
            feature_data.sort();
        }
        PostProcessor::Percentile(ref_scores)
    }
    /// Apply a z-score or percentile transformation
    /// to the given matrix, if requested.
    pub fn post_process<'b, 'c>(
        &self,
        matrix: FeatureMatrix<'b, i64>,
        arena: &'c Bump,
    ) -> Result<PostProcessedFeatureMatrix<'b>, Error> {
        let get_post_processor_feat_name_index = |post_processor_feat_names: &[&str],
                                                  feat_name: &str|
         -> Result<usize, Error> {
            let (idx, _) = post_processor_feat_names.iter().enumerate().find(|(_, pp_name)| **pp_name == feat_name).ok_or_else(|| Error::msg(format!("post-processor did not recognize this feature name from feature matrix: {}", feat_name)))?;
            Ok(idx)
        };
        match *self {
            Self::Raw => Ok(PostProcessedFeatureMatrix::Raw(matrix)),
            Self::ZScore {
                feature_names,
                sumfeature_stats,
                subfeature_stats,
            } => {
                let reordered_stats =
                    arena.alloc_slice_try_fill_iter(matrix.feature_names.iter().map(
                        |&feat_name| -> Result<[f64; 2], Error> {
                            let idx = get_post_processor_feat_name_index(feature_names, feat_name)?;
                            Ok(subfeature_stats[idx])
                        },
                    ))?;
                let num_features = matrix.feature_names.len();
                let row_size = matrix.feature_names.len() + 1;
                for row in matrix.data.chunks_exact_mut(row_size) {
                    for (slot, &[mean, std]) in
                        row[..num_features].iter_mut().zip(reordered_stats.iter())
                    {
                        let raw_value = *slot;
                        let zscore = (raw_value as f64 - mean) / std;
                        write_float_to_slot(slot, zscore);
                    }
                    let raw_value = row[num_features];
                    let [mean, std] = sumfeature_stats;
                    let zscore = (raw_value as f64 - mean) / std;
                    write_float_to_slot(&mut row[num_features], zscore);
                }
                Ok(PostProcessedFeatureMatrix::Processed(
                    matrix.cast_to_float(),
                ))
            }
            Self::Percentile(ReferenceFeatureMatrix {
                feature_names,
                num_ref_seqs,
                ref data,
            }) => {
                let reordered_ranges =
                    arena.alloc_slice_try_fill_iter(matrix.feature_names.iter().map(
                        |&feat_name| -> Result<std::ops::Range<usize>, Error> {
                            let idx = get_post_processor_feat_name_index(feature_names, feat_name)?;
                            let range = (idx * num_ref_seqs)..((idx + 1) * num_ref_seqs);
                            Ok(range)
                        },
                    ))?;
                let num_features = matrix.feature_names.len();
                let row_size = matrix.feature_names.len() + 1;
                let num_sequences = matrix.data.len() / row_size;
                for (j, range) in reordered_ranges.iter_mut().enumerate() {
                    let reference_values = &data[range.clone()];
                    for i in 0..num_sequences {
                        let slot_idx = i * row_size + j;
                        let raw_value = matrix.data[slot_idx];
                        let percentile = percentile(raw_value, reference_values);
                        write_float_to_slot(&mut matrix.data[slot_idx], percentile);
                    }
                }
                let reference_values = &data[feature_names.len() * num_ref_seqs..];
                for i in 0..num_sequences {
                    let slot_idx = i * row_size + num_features;
                    let raw_value = matrix.data[slot_idx];
                    let percentile = percentile(raw_value, reference_values);
                    write_float_to_slot(&mut matrix.data[slot_idx], percentile);
                }
                Ok(PostProcessedFeatureMatrix::Processed(
                    matrix.cast_to_float(),
                ))
            }
        }
    }
}

/// Utility to help re-use the integer buffer for floats.
fn write_float_to_slot(slot: &mut i64, x: f64) {
    *slot = x.to_bits().cast_signed()
}
/// Utility for computing a mean and standard deviation.
fn mean_std(iter: impl ExactSizeIterator<Item = i64>) -> [f64; 2] {
    let n = iter.len();
    let mut sum = 0;
    let mut sum_sqr = 0;
    for value in iter {
        sum += value;
        sum_sqr += value * value;
    }
    let mean = sum as f64 / n as f64;
    let var = (sum_sqr as f64 - mean * sum as f64) / (n - 1) as f64;
    [mean, var.sqrt()]
}
/// Utility for computing a percentile from a sorted list of reference values.
fn percentile(value: i64, reference_values: &[i64]) -> f64 {
    match reference_values.binary_search(&value) {
        Ok(idx) => {
            let mut num_less = 0;
            for i in (0..idx).rev() {
                if reference_values[i] < value {
                    num_less = i + 1;
                    break;
                }
            }
            let mut num_less_or_eq = reference_values.len();
            for i in (idx + 1)..reference_values.len() {
                if reference_values[i] > value {
                    num_less_or_eq = i;
                    break;
                }
            }
            ((num_less + num_less_or_eq + 1) * 50) as f64 / reference_values.len() as f64
        }
        Err(idx) => (idx * 100) as f64 / reference_values.len() as f64,
    }
}
