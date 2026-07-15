use std::path::Path;

fn main() {
    for directory in ["/lib64", "/usr/lib64", "/usr/lib/x86_64-linux-gnu"] {
        if Path::new(directory).exists() {
            println!("cargo:rustc-link-search=native={directory}");
        }
    }
}
