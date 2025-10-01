use std::path::{Path, PathBuf};

pub fn gen_target_fn(name: impl AsRef<Path>) -> PathBuf {
    let target_dir = std::env::var("CARGO_TARGET_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("target"));

    target_dir.join(name)
}
