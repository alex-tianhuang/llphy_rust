//! Datatypes for storing biophysical feature values
//! associated to sequence(s).

use std::slice::from_raw_parts_mut;
use std::fmt::Write;
use bumpalo::collections::String;

/// Feature values for lots of sequences,
/// along with associated feature names.
pub struct FeatureMatrix<'a, T> {
    pub feature_names: &'a [&'a str],
    /// `n_seqs x (n_features + 1)` buffer of values,
    /// grouped by sequence.
    ///
    /// The `+1` to `n_features` is a column
    /// for the sum of all features.
    pub data: &'a mut [T],
}
/// Feature values for lots of sequences,
/// that may or may not be post-processed.
///
/// Post-processing, such as converting values
/// into percentiles or zscores, turns values into
/// floats.
pub enum PostProcessedFeatureMatrix<'a> {
    Raw(FeatureMatrix<'a, i64>),
    Processed(FeatureMatrix<'a, f64>),
}
/// Feature values for one sequence,
/// that may or may not be post-processed.
///
/// Post-processing, such as converting values
/// into percentiles or zscores, turns values into
/// floats.
pub enum PostProcessedFeatureVector<'a> {
    Raw(&'a [i64]),
    Processed(&'a [f64]),
}
/// Integer feature values associated with
/// a reference set of sequences, grouped in memory
/// by feature rather than by sequence.
pub struct ReferenceFeatureMatrix<'a> {
    pub feature_names: &'a [&'a str],
    /// The number of reference sequences that
    /// this dataset is comprised of.
    pub num_ref_seqs: usize,
    /// `(n_features + 1) x num_ref_seqs` buffer of values.
    ///
    /// The `+1` to `n_features` is a column
    /// for the sum of all features.
    pub data: &'a mut [i64],
}
impl<'a> FeatureMatrix<'a, i64> {
    /// Assert that this matrix now contains floats.
    ///
    /// Not the same as `let x: i64; x as f64`.
    /// More like `let x: i64; mem::transmute(x)`.
    pub fn cast_to_float(self) -> FeatureMatrix<'a, f64> {
        let FeatureMatrix {
            feature_names,
            data,
        } = self;
        let data = unsafe { &mut *from_raw_parts_mut(data.as_mut_ptr().cast::<f64>(), data.len()) };
        FeatureMatrix {
            feature_names,
            data,
        }
    }
}
impl PostProcessedFeatureMatrix<'_> {
    /// Get feature names associated to this feature matrix.
    pub fn feature_names(&self) -> &[&str] {
        let (Self::Raw(FeatureMatrix { feature_names, .. }) | Self::Processed(FeatureMatrix { feature_names, .. })) = *self;
        feature_names
    }
    /// Equivalent to `self.data.len()`
    /// if the enum was not in the way.
    fn data_len(&self) -> usize {
        match self {
            Self::Raw(FeatureMatrix { data , .. }) => data.len(),
            Self::Processed(FeatureMatrix { data, .. }) => data.len()
        }
    }
    /// One more than the number of feature names.
    /// (because the last feature is an unnamed sum over the others).
    pub fn row_len(&self) -> usize {
        self.feature_names().len() + 1
    }
    /// Returns an iterator over feature vectors
    /// for each sequence.
    pub fn rows<'a>(&'a self) -> impl ExactSizeIterator<Item = PostProcessedFeatureVector<'a>> {
        let num_sequences = self.data_len() / self.row_len();
        debug_assert_eq!(self.data_len() % self.row_len(), 0);
        (0..num_sequences).map(move |i| match self {
            Self::Raw(FeatureMatrix { data , .. }) => PostProcessedFeatureVector::Raw(&data[i * self.row_len()..][..self.row_len()]),
            Self::Processed(FeatureMatrix { data, .. }) => PostProcessedFeatureVector::Processed(&data[i * self.row_len()..][..self.row_len()])
        })
    }
}
impl<'a> PostProcessedFeatureVector<'a> {
    /// Utility method for the [`crate::output`] module.
    /// 
    /// Write the values in this vector into a big comma separated list.
    pub fn format_comma_separated_into(&self, fmt: &mut String<'_>) {
        match *self {
            Self::Raw(data) => {
                let [subfeatures @ .., sumfeature] = data else {
                    panic!("expected non-empty data slice")
                };
                for value in subfeatures {
                    write!(fmt, "{},", value).unwrap();
                }
                write!(fmt, "{}", sumfeature).unwrap();
            }
            Self::Processed(data) => {
                let [subfeatures @ .., sumfeature] = data else {
                    panic!("expected non-empty data slice")
                };
                for value in subfeatures {
                    write!(fmt, "{},", value).unwrap();
                }
                write!(fmt, "{}", sumfeature).unwrap();
            }
        }
    }
}