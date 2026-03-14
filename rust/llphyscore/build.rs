fn main() {
    const EXPECT_MSG: &'static str = "the `llphyscore` bin was moved out of the workspace it was created in";
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let dev_pkg_data_root = manifest_dir
        .parent()
        .expect(EXPECT_MSG)
        .parent()
        .expect(EXPECT_MSG)
        .join("pkg_data");
    println!(
        "cargo:rustc-env=PKG_DATA_ROOT={}",
        dev_pkg_data_root.display()
    );
}
