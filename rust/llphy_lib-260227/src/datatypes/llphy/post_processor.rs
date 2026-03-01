use crate::{G2WScores, leak_vec};
use bumpalo::{Bump, collections::Vec};
use clap::ValueEnum;

/// The type of score to return.
#[derive(Clone, ValueEnum)]
#[value(verbatim_doc_comment)]
pub enum ScoreType {
    /// Report raw `llphyscore` values and features (integers).
    Raw,
    /// Report z-scores in comparison to a reference proteome. 
    ZScore,
    /// Report values in as percentiles compared to a reference proteome.
    Percentile,
    
}
pub enum PostProcessor<'a> {
    ZScore {
        feature_sum: (f64, f64),
        subfeatures: &'a [(f64, f64)],
    },
    Percentile {
        feature_sums: &'a [f64],
        subfeatures: &'a [&'a [f64]],
    },
    Raw,
}
impl<'a> PostProcessor<'a> {
    pub fn new_raw() -> Self {
        PostProcessor::Raw
    }
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
fn percentile(value: f64, reference_values: &[f64]) -> f64 {
    let (Ok(idx) | Err(idx)) = reference_values.binary_search_by(|probe| probe.total_cmp(&value));
    (idx * 100) as f64 / reference_values.len() as f64
}
