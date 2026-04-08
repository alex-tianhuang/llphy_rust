use crate::PKG_DATA_ROOT;
use bumpalo::Bump;
use llphyscore_core::{
    datatypes::{GridDecoderPair, ModelTrainingBase, PostProcessedFeatureVector},
    load_pkg_data::load_grid_decoders,
};
use pyo3::{Bound, PyResult, Python, types::{PyDict, PyDictMethods, PyString}};
use std::path::Path;
/// Helper macro for defining a static struct
/// from something that is constructed using an arena.
macro_rules! self_referential_static {
    (static ref $NAME:ident: $ty:ty = {
        let $arena:ident: &$lt:lifetime Bump = &$set_arena:expr;
        $($body:stmt);*
    }) => {
        use __private::$NAME;
        mod __private {
            use super::*;
            mod struct_def {
                use super::super::*;
                type T<$lt> = $ty;
                #[allow(non_camel_case_types)]
                pub struct $NAME<$lt> {
                    #[allow(unused)]
                    pub(super) arena: bumpalo::Bump,
                    pub(super) data: std::mem::ManuallyDrop<$ty>
                }
                unsafe impl Sync for $NAME<'static> {}
                impl std::ops::Deref for $NAME<'static> {
                    type Target = T<'static>;
                    fn deref(&self) -> &Self::Target {
                        &self.data
                    }
                }
            }
            lazy_static::lazy_static! {
                pub static ref $NAME: struct_def::$NAME<'static> = {
                    let arena = $set_arena;
                    let $arena = unsafe { &*(&arena as *const bumpalo::Bump) };
                    let data = {
                        $($body);*
                    };
                    struct_def::$NAME {
                        arena,
                        data: std::mem::ManuallyDrop::new(data)
                    }
                };
            }
        }
    };
}
/// Load [`GridDecoderPair`]s and cache them in static memory.
pub fn load_grid_decoders_static(
    py: Python,
    model_training_base: ModelTrainingBase,
) -> PyResult<&'static [(&'static str, GridDecoderPair)]> {
    macro_rules! impl_branches {
        ($model_training_base:ident, $($branch:ident),+) => {
            match $model_training_base {
                $(ModelTrainingBase::$branch => {
                    self_referential_static! {
                        static ref GRID_DECODERS: PyResult<&'a [(&'a str, GridDecoderPair)]> = {
                            let arena: &'a Bump = &Bump::new();
                            match load_grid_decoders(Path::new(PKG_DATA_ROOT), ModelTrainingBase::$branch, arena) {
                                Ok(decoders) => Ok(decoders),
                                Err(e) => Err(e.into())
                            }
                        }
                    }
                    Ok(&GRID_DECODERS.as_ref().map_err(|e| e.clone_ref(py))?)
                }),+
            }
        };
    }
    impl_branches!(model_training_base, Human, HumanPDB, PDB)
}

/// Turns a row of numbers + feature names
/// into a dictionary of {feat_name, feat_value} pairs.
/// 
/// Fails if and only if `dict.__setitem__` fails, so typically does not fail.
pub fn serialize_row<'py>(
    row: PostProcessedFeatureVector,
    feature_names: &[Bound<PyString>],
    py: Python<'py>,
) -> PyResult<Bound<'py, PyDict>> {
    match row {
        PostProcessedFeatureVector::Raw(data) => {
            debug_assert_eq!(feature_names.len(), data.len());
            let output_row = PyDict::new(py);
            for (feat_name, value) in feature_names.iter().zip(data) {
                output_row.set_item(feat_name, value)?;
            }
            Ok(output_row)
        }
        PostProcessedFeatureVector::Processed(data) => {
            debug_assert_eq!(feature_names.len(), data.len());
            let output_row = PyDict::new(py);
            for (feat_name, value) in feature_names.iter().zip(data) {
                output_row.set_item(feat_name, value)?;
            }
            Ok(output_row)
        }
    }
}
