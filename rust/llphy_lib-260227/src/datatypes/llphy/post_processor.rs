//! Module defining structs that turn raw biophysical feature scores
//! to z-scores or percentiles if requested.
//! 
//! See [`ScoreType`] and [`PostProcessor`].
use crate::{G2WScores, leak_vec};
use anyhow::Error;
use bumpalo::{Bump, collections::Vec};
use clap::ValueEnum;
use pyo3::FromPyObject;

/// The type of score to return.
#[derive(Clone, ValueEnum)]
#[value(verbatim_doc_comment)]
pub enum ScoreType {
    /// Report raw `llphyscore` values and features.
    Raw,
    /// Report z-scores in comparison to a reference proteome. 
    ZScore,
    /// Report values in as percentiles compared to a reference proteome.
    Percentile,
}
impl<'a, 'py> FromPyObject<'a, 'py> for ScoreType {
    type Error = Error;
    fn extract(obj: pyo3::Borrowed<'a, 'py, pyo3::PyAny>) -> Result<Self, Self::Error> {
        let s = obj.extract()?;
        Ok(ScoreType::from_str(s, false).map_err(Error::msg)?)
    }
}
/// A data struct derived from a reference dataset of
/// phase separation positive/negative proteins.
/// 
/// Depending on what score-type is asked for (see [`ScoreType`]),
/// the post-processing will require different data.
pub enum PostProcessor<'a> {
    /// No post-processing required.
    Raw,
    /// Report z-scores by normalizing using `(mean, std)` tuples.
    ZScore {
        feature_sum: (f64, f64),
        subfeatures: &'a [(f64, f64)],
    },
    /// Report percentiles by binary searching for the percentile
    /// in a sorted array of biophysical feature scores on the
    /// reference protein dataset.
    Percentile {
        feature_sums: &'a [f64],
        subfeatures: &'a [&'a [f64]],
    },
}
impl<'a> PostProcessor<'a> {
    /// Make a new [`PostProcessor`] that does not do any post-processing of
    /// biophysical feature scores.
    pub fn new_raw() -> Self {
        PostProcessor::Raw
    }
    /// Make a new [`PostProcessor`] that converts biphysical feature scores
    /// to z-scores.
    pub fn new_zscore(
        data: Vec<'_, G2WScores<'_>>,
        num_features: usize,
        arena: &'a Bump,
    ) -> Self {
        let subfeature_zscore_data =
            (0..num_features).map(|i| mean_std(data.iter().map(|entry| entry.subfeatures[i])));
        let subfeatures = Vec::from_iter_in(subfeature_zscore_data, arena);
        PostProcessor::ZScore {
            feature_sum: mean_std(data.iter().map(|entry| entry.feature_sum)),
            subfeatures: leak_vec(subfeatures),
        }
    }
    /// Make a new [`PostProcessor`] that converts biphysical feature scores
    /// to percentiles.
    pub fn new_percentile(
        data: Vec<'_, G2WScores<'_>>,
        num_features: usize,
        arena: &'a Bump,
    ) -> Self {
        let mut feature_sums = Vec::with_capacity_in(data.len(), arena);
        let subfeature_zscore_data =
            (0..num_features).map(|_| Vec::with_capacity_in(data.len(), arena));
        let mut subfeature_dests = Vec::from_iter_in(subfeature_zscore_data, arena);
        for entry in data {
            feature_sums.push(entry.feature_sum);
            for (subfeature_dest, subfeature_value) in
                subfeature_dests.iter_mut().zip(entry.subfeatures)
            {
                subfeature_dest.push(*subfeature_value)
            }
        }
        feature_sums.sort_by(f64::total_cmp);
        for subfeature_dest in subfeature_dests.iter_mut() {
            subfeature_dest.sort_by(f64::total_cmp);
        }
        let subfeatures = Vec::from_iter_in(
            subfeature_dests.into_iter().map(|v| leak_vec(v) as &[_]),
            arena,
        );
        PostProcessor::Percentile {
            feature_sums: leak_vec(feature_sums),
            subfeatures: leak_vec(subfeatures),
        }
    }
    /// Given raw scores as `data`, turn them into
    /// the requested score-type, in-place.
    pub fn transform(&self, data: &mut G2WScores<'_>) {
        match *self {
            Self::Raw => (),
            Self::ZScore {
                feature_sum,
                subfeatures,
            } => {
                data.feature_sum = (data.feature_sum - feature_sum.0) / feature_sum.1;
                for (subfeature, &(subfeature_mean, subfeature_std)) in
                    data.subfeatures.iter_mut().zip(subfeatures)
                {
                    *subfeature = (*subfeature - subfeature_mean) / subfeature_std;
                }
            }
            Self::Percentile {
                feature_sums,
                subfeatures,
            } => {
                data.feature_sum = percentile(data.feature_sum, feature_sums);
                for (subfeature, &reference_subfeatures) in
                    data.subfeatures.iter_mut().zip(subfeatures)
                {
                    *subfeature = percentile(*subfeature, reference_subfeatures);
                }
            }
        }
    }
}

/// Utility for computing a mean and standard deviation.
fn mean_std(iter: impl ExactSizeIterator<Item = f64>) -> (f64, f64) {
    let n = iter.len();
    let mut sum = 0.0;
    let mut sum_sqr = 0.0;
    for value in iter {
        sum += value;
        sum_sqr += value * value;
    }
    let mean = sum / n as f64;
    let var = (sum_sqr - mean * sum) / (n - 1) as f64;
    (mean, var.sqrt())
}
/// Utility for computing a percentile from a sorted list of reference values.
fn percentile(value: f64, reference_values: &[f64]) -> f64 {
    let (Ok(idx) | Err(idx)) = reference_values.binary_search_by(|probe| probe.total_cmp(&value));
    (idx * 100) as f64 / reference_values.len() as f64
}
