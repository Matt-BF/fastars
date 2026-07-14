# fastars

`fastars` fetches records from a large BGZF-compressed FASTA using exact full
IDs or optional short IDs. It writes FASTA records to standard output, so it
can be used directly in shell pipelines.

## Requirements

- Rust and Cargo to build the program.
- A BGZF-compressed FASTA (`.bgz`), not plain gzip or zstd compression.
- Matching `.fai` and `.gzi` files for the FASTA.
- The system `sort` command to build an `.ffx` lookup index.

Create missing FASTA indexes with samtools, for example:

```bash
samtools faidx sequences.fna.bgz
```

## Build

```bash
cargo build --release
```

The executable is `target/release/fastars`.

## Build a lookup index

An `.ffx` index is sorted by ID and points to the matching `.fai` entries.
Full FASTA IDs are always indexed. Add a regex when short-ID lookup is also
needed; capture group 1 becomes the short ID.

```bash
fastars index \
  --fai sequences.fna.bgz.fai \
  --short-id-regex '^([^|]+)'
```

For a header such as:

```text
IMGVR_UViG_2582581227_000001|2582581227|2582690522
```

the example regex creates the short key
`IMGVR_UViG_2582581227_000001`.

Use `--output PATH.ffx` to choose a different index location and
`--temp-directory DIR` to choose where the temporary sort input is written.

## Fetch records

Full-ID lookup is the default:

```bash
fastars --fasta sequences.fna.bgz \
  'IMGVR_UViG_2582581227_000001|2582581227|2582690522' > selected.fna
```

Use `--short-id` only when the `.ffx` was built with `--short-id-regex`:

```bash
fastars --fasta sequences.fna.bgz \
  --short-id IMGVR_UViG_2582581227_000001 > selected.fna
```

Supply multiple IDs as arguments or one per line with `--ids-file`:

```bash
fastars --fasta sequences.fna.bgz --ids-file requested_ids.txt > selected.fna
```

The default paths are derived from `--fasta`:

- `.fai`: `FASTA.bgz.fai`
- `.gzi`: `FASTA.bgz.gzi`
- `.ffx`: `FASTA.bgz.fai.ffx`

Override them with `--fai`, `--gzi`, and `--ffx`. Use `--sort-by-offset` for
more efficient disk access when fetching many records; output order then
follows FASTA position rather than request order.

## Notes

- `.ffx` files are generated artifacts and should be rebuilt after changing
  the short-ID regex or upgrading from an older `fastars` index format.
- The current `.ffx` gets offsets from the `.fai`, so keep the `.fai` and
  `.gzi` files alongside the BGZF FASTA.
- Plain `.gz` and `.zst` FASTA files are not supported for random retrieval.
