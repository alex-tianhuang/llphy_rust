use std::{
    env,
    path::{Path, PathBuf},
    process::Command,
};

/// Python discovery script written by Claude Sonnet 4.6.
fn main() {
    // Re-run if any of the env vars that influence Python discovery change.
    println!("cargo:rerun-if-env-changed=PYO3_PYTHON");
    println!("cargo:rerun-if-env-changed=VIRTUAL_ENV");
    println!("cargo:rerun-if-env-changed=CONDA_PREFIX");

    // -------------------------------------------------------------------------
    // NOTE: We deliberately do NOT try to set PYO3_PYTHON here.
    //
    // PyO3's own build.rs runs in a separate process and cannot see env vars
    // set with std::env::set_var(), and cargo:rustc-env only affects the
    // compiled binary. PyO3 does its own discovery (active venv → python →
    // python3) from the environment that existed when `cargo` was invoked.
    //
    // If you need to pin a specific interpreter, set PYO3_PYTHON in your shell
    // before running cargo:
    //
    //   PYO3_PYTHON=/path/to/python cargo build
    //
    // or put it in .cargo/config.toml:
    //
    //   [env]
    //   PYO3_PYTHON = "/path/to/python"
    //
    // What this build script *does* do:
    //   1. Mirror PyO3's discovery so we know which interpreter will be used.
    //   2. Verify numpy is installed in that interpreter.
    //   3. On macOS, bake an rpath so `cargo run` / the binary can find
    //      libpython at runtime (dyld doesn't use linker search paths).
    // -------------------------------------------------------------------------

    let python = discover_python().expect(
        "Could not find a suitable Python interpreter.\n\
         Activate a venv/conda env, or set PYO3_PYTHON before invoking cargo.",
    );
    eprintln!("[build.rs] Using Python: {}", python.display());

    verify_numpy(&python);

    // On macOS, bake the directory containing libpython as an rpath so that
    // `cargo run` works without needing DYLD_LIBRARY_PATH.
    // On Linux the linker already embeds RUNPATH correctly; on Windows the DLL
    // must be on PATH at runtime (standard Windows behaviour).
    if cfg!(target_os = "macos") {
        let lib_dir = python_lib_dir(&python);
        eprintln!("[build.rs] Baking macOS rpath: {}", lib_dir.display());
        println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib_dir.display());
    }
}

// ---------------------------------------------------------------------------
// Python discovery — mirrors PyO3's own priority order so we inspect the
// same interpreter that PyO3 will link against.
// ---------------------------------------------------------------------------

fn discover_python() -> Option<PathBuf> {
    // 1. Explicit override (same var PyO3 reads).
    if let Ok(p) = env::var("PYO3_PYTHON") {
        let path = PathBuf::from(p);
        return is_usable_python(&path).then_some(path);
    }

    // 2. Active venv (PEP 405 / virtualenv).
    if let Ok(root) = env::var("VIRTUAL_ENV") {
        let p = python_binary_in(Path::new(&root));
        if is_usable_python(&p) {
            return Some(p);
        }
    }

    // 3. Active conda / mamba env (both set CONDA_PREFIX).
    if let Ok(root) = env::var("CONDA_PREFIX") {
        let p = python_binary_in(Path::new(&root));
        if is_usable_python(&p) {
            return Some(p);
        }
    }

    // 4. python3 / python on PATH (same fallback order PyO3 uses).
    for name in &["python3", "python"] {
        if let Some(p) = which(name) {
            if is_usable_python(&p) {
                return Some(p);
            }
        }
    }

    None
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Platform-appropriate Python binary path inside an env root.
fn python_binary_in(root: &Path) -> PathBuf {
    if cfg!(windows) {
        root.join("Scripts").join("python.exe")
    } else {
        root.join("bin").join("python3")
    }
}

/// Returns true if the path is an executable Python >= 3.7.
fn is_usable_python(python: &Path) -> bool {
    Command::new(python)
        .args(["-c", "import sys; assert sys.version_info >= (3, 7)"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Hard-fail the build if numpy is not importable.
fn verify_numpy(python: &Path) {
    let output = Command::new(python)
        .args(["-c", "import numpy; print(numpy.__version__)"])
        .output()
        .unwrap_or_else(|e| panic!("Failed to invoke {}: {e}", python.display()));

    if output.status.success() {
        eprintln!(
            "[build.rs] numpy {} ✓",
            String::from_utf8_lossy(&output.stdout).trim()
        );
    } else {
        panic!(
            "numpy not found in {}.\n\
             Run: pip install numpy\n\
             {}",
            python.display(),
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

/// Ask sysconfig for the directory containing libpython (used for macOS rpath).
fn python_lib_dir(python: &Path) -> PathBuf {
    let out = Command::new(python)
        .args([
            "-c",
            "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))",
        ])
        .output()
        .expect("failed to query Python LIBDIR");
    PathBuf::from(String::from_utf8_lossy(&out.stdout).trim().to_string())
}

/// Minimal `which`: search PATH for an executable by name.
fn which(name: &str) -> Option<PathBuf> {
    let path_var = env::var("PATH").ok()?;
    let exe = if cfg!(windows) {
        format!("{name}.exe")
    } else {
        name.to_string()
    };
    env::split_paths(&path_var)
        .map(|dir| dir.join(&exe))
        .find(|p| p.is_file())
}