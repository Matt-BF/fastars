# Local index benchmark

`data/METAVR_sorted.2GiB.fna.bgz` is a local 2 GiB prefix of
`METAVR_sorted.fna.bgz`. It ends at a complete BGZF block boundary and is
intended only for performance testing; its final FASTA record may be partial.
The data file is ignored by Git.

- Compressed size: 2,147,486,355 bytes
- BGZF blocks: 160,605
- Indexed records: 263,410
- SHA-256: `0a91d472ace40d651bff10f056605a44eeb6eff8d3acf453d4ec50d2e22f293c`

Build the release binary before benchmarking:

```bash
cargo build --release
```

Benchmark the default in-memory sorting path:

```bash
hyperfine --warmup 1 \
  'target/release/fastars index \
    --fasta benchmarks/data/METAVR_sorted.2GiB.fna.bgz \
    --output /tmp/fastars-metavr-sample.ffx \
    --temp-directory /tmp \
    --threads 32 \
    --sort-memory 512'
```

Use a 64 MiB budget to exercise the external-sort fallback:

```bash
hyperfine --warmup 1 \
  'target/release/fastars index \
    --fasta benchmarks/data/METAVR_sorted.2GiB.fna.bgz \
    --output /tmp/fastars-metavr-sample.ffx \
    --temp-directory /tmp \
    --threads 32 \
    --sort-memory 64'
```
