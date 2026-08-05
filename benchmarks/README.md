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

## BGZF line-scanning result

Replacing byte-at-a-time line reads with block-slice scans produced the
following local results with `--threads 32 --sort-memory 512`. Wall time varies
with the shared filesystem, while user CPU is stable.

| Reader | Wall time (s) | User CPU (s) | Peak RSS (KiB) |
| --- | ---: | ---: | ---: |
| Byte-at-a-time | 105.18 | 54.33 | 110,252 |
| Byte-at-a-time | 123.81 | 55.32 | 103,996 |
| Byte-at-a-time, post-control | 59.15 | 54.53 | 120,732 |
| Block-slice scan | 39.67 | 37.68 | 106,796 |
| Block-slice scan | 37.64 | 37.78 | 106,404 |

The block-slice implementation reduced user CPU by about 31%. Against the
immediately following pre-change control, elapsed time fell by about 35%. The
generated `.ffx` files were byte-identical.

## Block-aware indexing scanner result

The dedicated indexing scanner consumes decompressed BGZF blocks directly and
copies only lines that cross a block boundary. A four-thread run reduced user
CPU by another 10% while producing a byte-identical `.ffx` file.

| Scanner | Wall time (s) | User CPU (s) | Peak RSS (KiB) |
| --- | ---: | ---: | ---: |
| Allocated line scan | 44.06 | 36.75 | 49,924 |
| Block-aware index scan | 44.63 | 32.94 | 46,408 |

The shared filesystem dominated elapsed time during these runs, so user CPU is
the more reliable measure of the parser improvement.
