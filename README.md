# fastars

`fastars` fetches records from large BGZF-compressed or uncompressed FASTA files using a
self-contained `.ffx` index. It writes FASTA (nucleotide and protein) records
to standard output, so it fits directly into shell pipelines.

Use `fastars index` to build an index and `fastars fetch` to retrieve records.

## Requirements if building from source

- Rust and Cargo to build the program.
- A C compiler for the bundled zstd library.
- A BGZF-compressed (`.bgz`) or uncompressed FASTA. Plain gzip and zstd
  compression are not supported.

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

Tagged releases are built for Linux x86-64 and macOS on both Apple Silicon and
Intel.

## Build an index

The preferred path scans the BGZF-compressed or uncompressed FASTA directly:

```bash
fastars index --fasta sequences.fna.bgz
```

Set `--sort-memory <MiB>` to change the budget; larger indexes automatically
fall back to a platform-independent external merge sort in `--temp-directory`.
Up to 512 MiB is used for an in-memory ID sort by default.
Temporary records are stored in compact binary runs, with each sort run bounded
by the requested memory budget. BGZF decompression and block encoding use all
available CPUs by default; use `--threads <N>` to set the worker count.

This writes:

```text
sequences.fna.bgz.ffx
```

If existing samtools indexes are available, they can be used as a build
accelerator but are not required and can be deleted after creating the `.ffx`
file:

```bash
fastars index \
  --fai sequences.fna.bgz.fai \
  --gzi sequences.fna.bgz.gzi \
  --output sequences.fna.bgz.ffx
```

The resulting `.ffx` is a compressed, self-contained fetch index. It stores
primary IDs, complete FASTA headers, BGZF virtual offsets, sequence lengths,
and FASTA line layout in independently compressed blocks for fast lookup
without loading the complete index. Building from `.fai/.gzi` cannot preserve
header descriptions because `.fai` contains only primary IDs.

After building the index, use `--id-mode prefix` to fetch IDs by literal prefix
or `--id-regexp` to select indexed IDs with a regular expression. Examples for
both modes are below.

## Fetch by exact ID or header

Exact lookup is the default. A primary ID query returns the record with its
complete original header:

```bash
fastars fetch --fasta sequences.fna.bgz \
  'IMGVR_UViG_2582581227_000001|2582581227|2582690522' > selected.fna
```

Given this header:

```text
>ABC123 hypothetical protein [Species name]
```

both `ABC123` and the quoted complete header are valid exact queries:

```bash
fastars fetch --fasta sequences.fna.bgz ABC123
fastars fetch --fasta sequences.fna.bgz \
  'ABC123 hypothetical protein [Species name]'
```

If multiple records share a primary ID, the short query returns all of them;
the complete-header query can select one specific record.

If `sequences.fna.bgz.ffx` is missing, `fastars` builds it automatically from
the FASTA before fetching.

## Fetch by prefix

Use prefix mode when your query is the beginning of a primary ID or complete
header:

```bash
fastars fetch --fasta sequences.fna.bgz \
  --id-mode prefix IMGVR_UViG_2582581227_000001 > selected.fna
```

For headers like:

```text
IMGVR_UViG_2582581227_000001|2582581227|2582690522
```

the query `IMGVR_UViG_2582581227_000001` matches because it is a literal prefix.
For the descriptive `ABC123` example above, the prefix
`ABC123 hypothetical` also matches.

## Fetch from an ID file

Use `-f` or `--ids-file` for one query per line:

```bash
fastars fetch --fasta sequences.fna.bgz \
  --id-mode prefix \
  -f short_ids.txt > selected.fna
```

With `--id-mode exact`, each line may be a primary ID or complete header. With
`--id-mode prefix`, each line is treated as a literal prefix of an ID or header.
`-m` is the short form of `--id-mode`.

For explicit query results spread across the FASTA, use `--fetch-threads <N>`
(`N` from 1 through 16) to fetch bounded batches in parallel while preserving
output order. Regex-only results remain sequential so they can stream without
buffering the full result set.

## Search complete headers with regex

`--id-regexp` scans complete indexed headers, including descriptions, but not
the FASTA sequence text:

```bash
fastars fetch --fasta sequences.fna.bgz \
  --id-regexp 'GVMAG' > gvmag_records.fna
```

Invert the regex to fetch everything whose header does not match:

```bash
fastars fetch --fasta sequences.fna.bgz \
  --id-regexp 'GVMAG' \
  -v > non_gvmag_records.fna
```

`-v` and `--invert-match` are equivalent.
`-r` is the short form of `--id-regexp`.

Regex mode is useful for broad metadata-style searches. For large ID lists,
prefer exact or prefix lookup with `-f` because it uses binary search over the
sorted ID index.

## Output order

By default, exact and prefix results follow query order, and regex results
follow sorted ID order. Use `--sort-by-offset` to fetch in FASTA order, which
can reduce random disk access for many records. `-s` is its short form:

```bash
fastars fetch --fasta sequences.fna.bgz \
  --id-mode prefix \
  -f short_ids.txt \
  --sort-by-offset > selected.fna
```

## Notes

- `.ffx` is a generated artifact. Rebuild it after changing the FASTA or
  upgrading from an older `fastars` index format. Older uncompressed indexes
  are not compatible with the current format.
- Plain `.gz` and `.zst` FASTA files are not supported for random retrieval.
