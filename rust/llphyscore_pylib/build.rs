use std::path::Path;
fn main() {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let dev_pkg_data_root = manifest_dir
        .parent()
        .expect("the `llphyscore` bin was moved out of the workspace it was created in")
        .join("llphy_lib-260227")
        .join("pkg_data");
    println!(
        "cargo:rustc-env=PKG_DATA_ROOT={}",
        dev_pkg_data_root.display()
    );
}
