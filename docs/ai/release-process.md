# Release process

## Version policy

`Cargo.toml` is the source of truth for the fastars version. Release tags use
the same version with a `v` prefix, such as `v0.1.0`. The release workflow
rejects a tag that does not match `Cargo.toml`.

While fastars is below `1.0.0`:

- Increment the patch version for compatible bug fixes and performance work.
- Increment the minor version for new features or intentional CLI/index-format
  incompatibilities.

After `1.0.0`, follow Semantic Versioning: patch for compatible fixes, minor
for compatible features, and major for incompatible behavior.

The crate version and `.ffx` format version are separate. If an incompatible
`.ffx` change is necessary, increment the format version in
`src/ffx/format.rs` as well. Older indexes will then produce a rebuild error.

## Automated builds

`.github/workflows/release.yml` builds these archives:

- Linux x86-64
- macOS Apple Silicon
- macOS Intel

Each archive contains the `fastars` executable and README, with a corresponding
SHA-256 file. A manual workflow run uploads temporary Actions artifacts without
creating a GitHub release. Pushing a matching `v*` tag also creates a GitHub
release and attaches all archives and checksums.

The workflow uses bundled zstd code. Linux and macOS provide the zlib library
used for BGZF decompression.

## Release checklist

1. Choose the new version and update `version` in `Cargo.toml`.
2. Run `cargo check` so the root package version in `Cargo.lock` is updated.
3. Run the release checks:

   ```bash
   cargo test --locked
   cargo clippy --all-targets -- -D warnings
   cargo build --release --locked
   ```

4. Commit the version change, for example:

   ```bash
   git add Cargo.toml Cargo.lock
   git commit -m "Release fastars 0.2.0"
   ```

5. Create and push an annotated tag:

   ```bash
   git tag -a v0.2.0 -m "fastars 0.2.0"
   git push origin HEAD
   git push origin v0.2.0
   ```

6. Confirm the `Build release binaries` workflow succeeds and inspect the
   generated GitHub release assets.

If the release job cannot write releases, enable read/write workflow
permissions under the repository or organization Actions settings.
