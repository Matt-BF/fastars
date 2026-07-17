# fastars

`fastars` fetches records from large BGZF-compressed FASTA files using a
self-contained `.ffx` index. It writes FASTA records to standard output, so it
fits directly into shell pipelines.

## Requirements

- Rust and Cargo to build the program.
- A BGZF-compressed FASTA (`.bgz`), not plain gzip or zstd compression.
- The system `sort` command to build an `.ffx` index.

Fetch does not require `.fai` or `.gzi` files. Those files are only optional
inputs for building `.ffx` faster when they already exist.

## Build

```bash
cargo build --release
```

The executable is:

```bash
target/release/fastars
```

## Build an index

The preferred path scans the BGZF FASTA directly:

```bash
fastars index --fasta sequences.fna.bgz
```

This writes:

```text
sequences.fna.bgz.ffx
```

If existing samtools indexes are available, they can be used as a build
accelerator:

```bash
fastars index \
  --fai sequences.fna.bgz.fai \
  --gzi sequences.fna.bgz.gzi \
  --output sequences.fna.bgz.ffx
```

The resulting `.ffx` is the same kind of self-contained fetch index. It stores
full IDs, BGZF virtual offsets, sequence lengths, and FASTA line layout.

After building the index, use `--id-mode prefix` to fetch IDs by literal prefix
or `--id-regexp` to select indexed IDs with a regular expression. Examples for
both modes are below.

Useful build options:

```bash
fastars index --fasta sequences.fna.bgz \
  --output custom.ffx \
  --temp-directory /scratch/tmp
```

## Fetch by exact full ID

Exact full-ID lookup is the default:

```bash
fastars --fasta sequences.fna.bgz \
  'IMGVR_UViG_2582581227_000001|2582581227|2582690522' > selected.fna
```

If `sequences.fna.bgz.ffx` is missing, `fastars` builds it automatically from
the FASTA before fetching.

## Fetch by prefix

Use prefix mode when your query is the beginning of the indexed full ID:

```bash
fastars --fasta sequences.fna.bgz \
  --id-mode prefix IMGVR_UViG_2582581227_000001 > selected.fna
```

For headers like:

```text
IMGVR_UViG_2582581227_000001|2582581227|2582690522
```

the query `IMGVR_UViG_2582581227_000001` matches because it is a literal prefix.
This is global and does not assume `|` separators.

## Fetch from an ID file

Use `-f` or `--ids-file` for one query per line:

```bash
fastars --fasta sequences.fna.bgz \
  --id-mode prefix \
  -f short_ids.txt > selected.fna
```

With `--id-mode exact`, each line must be a full exact ID. With
`--id-mode prefix`, each line is treated as a literal prefix.

## Search indexed IDs with regex

`--id-regexp` scans the indexed full IDs, not the FASTA sequence text:

```bash
fastars --fasta sequences.fna.bgz \
  --id-regexp 'GVMAG' > gvmag_records.fna
```

Invert the regex to fetch everything whose full ID does not match:

```bash
fastars --fasta sequences.fna.bgz \
  --id-regexp 'GVMAG' \
  --invert-match > non_gvmag_records.fna
```

Regex mode is useful for broad metadata-style searches. For large ID lists,
prefer exact or prefix lookup with `-f` because it uses binary search over the
sorted ID index.

## Output order

By default, exact and prefix results follow query order, and regex results
follow sorted ID order. Use `--sort-by-offset` to fetch in FASTA order, which
can reduce random disk access for many records:

```bash
fastars --fasta sequences.fna.bgz \
  --id-mode prefix \
  -f short_ids.txt \
  --sort-by-offset > selected.fna
```

## Notes

- `.ffx` is a generated artifact. Rebuild it after changing the FASTA or
  upgrading from an older `fastars` index format.
- Fetching needs only the BGZF FASTA and `.ffx`; `.fai` and `.gzi` are not
  read during fetch.
- Plain `.gz` and `.zst` FASTA files are not supported for random retrieval.
