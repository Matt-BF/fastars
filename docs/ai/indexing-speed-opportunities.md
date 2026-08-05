# Index creation speed opportunities

## Summary

The next major indexing gains are likely in FASTA parsing and BGZF
decompression, not in replacing the sorted `.ffx` lookup structure with a
trie. Building the index must still read, decompress, and inspect the entire
FASTA, and that work dominates the relatively small final index.

A block-aware FASTA scanner now avoids allocating and copying every sequence
line. External-sort improvements become important when the ID records do not
fit within the configured sort-memory budget.

## Benchmark input

The local benchmark uses:

- File: `benchmarks/data/METAVR_sorted.2GiB.fna.bgz`
- Compressed size: 2,147,486,355 bytes
- BGZF blocks: 160,605
- Indexed records: 263,410
- Output `.ffx` size: 3,888,512 bytes

The sample is a complete-BGZF-block prefix of the full METAVR file. Its final
FASTA record may be partial, so it is intended only for performance testing.

## Thread-count sweep

The following single-run sweep used the release binary with a 512 MiB sort
budget. Shared-filesystem load varies substantially, so elapsed times are
directional rather than precise.

| Threads | Wall time | User CPU | Peak RSS |
| ---: | ---: | ---: | ---: |
| 1 | 121.54 s | 35.08 s | 38,060 KiB |
| 2 | 77.81 s | 36.24 s | 47,888 KiB |
| 4 | 44.06 s | 36.75 s | 49,924 KiB |
| 8 | 46.56 s | 37.20 s | 66,864 KiB |
| 16 | 58.39 s | 37.48 s | 81,996 KiB |
| 32 | 59.91 s | 38.63 s | 114,924 KiB |

An earlier 32-thread run completed in 37.64 seconds, demonstrating the size of
the storage-related variance. However, the stable aggregate CPU time and
increasing memory use show that large thread counts do not reduce the amount
of work. Four to eight decompression workers appear sufficient for this
workload, pending repeated interleaved measurements.

The previously implemented decompressed-block line scan reduced user CPU from
about 54.7 seconds to 37.7 seconds, a reduction of approximately 31%.

## Prioritized improvements

### 1. Block-aware, near-zero-copy FASTA parsing

The current reader constructs a new `Vec<u8>` for every FASTA line and copies
the line out of the decompressed BGZF block. Sequence lines are then scanned
again to validate whitespace and determine their layout. Sequence bases are
not stored in `.ffx`, so these allocations and copies are avoidable.

A dedicated indexing scanner should:

- Parse directly from decompressed block slices.
- Locate newlines and headers without copying sequence data.
- Allocate only IDs and rare lines that cross BGZF block boundaries.
- Calculate sequence lengths and line wrapping while walking bytes once.
- Preserve the current input validation and virtual-offset behavior.

This change reduced user CPU from 36.75 seconds to 32.94 seconds in a
four-thread sample run, about 10%, and produced a byte-identical `.ffx` file.
Wall time remained storage-bound during the comparison.

### 2. Reuse or replace BGZF decompression state

Every BGZF block currently allocates compressed and decompressed buffers and
initializes and destroys a zlib stream. Workers could instead retain their
zlib state and recycle buffers between jobs.

A separate experimental branch could benchmark `libdeflate` against zlib on
METAVR. It introduces a dependency and should only be adopted if it produces a
meaningful measured improvement.

### 3. Improve the external-sort path

When records exceed the in-memory budget, fastars writes textual temporary
records and invokes GNU `sort`. The current `--sort-memory` value determines
when spilling starts, but is not explicitly passed to GNU `sort` with `-S`.
The configured thread count is also not passed with `--parallel`.

The improvements can be staged:

1. Pass explicit memory and thread limits to GNU `sort`, then benchmark the
   spill path.
2. Replace textual temporary records with a compact binary representation to
   reduce formatting, parsing, and temporary I/O.
3. If sorting remains significant, radix-partition IDs into memory-sized runs
   and merge or concatenate the partitions in lexical order.

For the full METAVR file, a larger budget such as `--sort-memory 4096` is
expected to keep the projected 8-9 million index records in memory, subject to
verification on the complete input.

### 4. Analyze decompressed blocks in workers

Decompression workers could also find newline positions, invalid whitespace,
and header candidates. The ordered serial stage would then only join lines or
records that cross block boundaries and emit metadata.

This could parallelize more of the parser, but it is more complex because
FASTA lines and records may span BGZF blocks. It should follow the simpler
zero-copy scanner so its incremental value can be measured.

### 5. Reduce ID allocation and cloning

IDs are currently allocated or cloned at several stages for sortedness
tracking, record ownership, and block boundaries. Passing owned strings through
the builder and retaining references where possible would reduce allocator
work. This is more relevant to record-dense FASTA files than sequence-heavy
METAVR data.

### 6. Tune and separate thread counts

Using every available CPU is not automatically fastest. After repeated
benchmarks, consider:

- Capping automatic BGZF decompression workers near the observed useful range.
- Configuring decompression and `.ffx` encoding workers separately.
- Avoiding large encoding pools when the output contains few `.ffx` blocks.

### 7. Optimize the `.fai`/`.gzi` build route

The `.fai`/`.gzi` path currently binary-searches `.gzi` using repeated file
seeks for every FASTA record. Loading the `.gzi` entries once and searching an
in-memory array would substantially reduce small random reads. This does not
accelerate direct FASTA scanning.

## Index structure conclusion

A conventional trie would remove the requirement to sort IDs, but would not
reduce the dominant scan and decompression work. It would also introduce node
and pointer overhead and poorer locality unless implemented as a specialized
compressed radix structure.

If eliminating sorting becomes more important than ordered prefix lookup, a
hash-bucketed on-disk index or radix-partitioned build is a better candidate
than a conventional trie. The existing sorted block format remains attractive
for compact front coding, binary search, prefix queries, and predictable disk
access.

## Recommended implementation order

Keep each item in a focused commit and benchmark it independently:

1. Tune the GNU `sort` invocation for the spill path.
2. Reuse decompression state and buffers.
3. Benchmark `libdeflate` in an experimental branch.
4. Revisit parallel block analysis only if parsing remains dominant.
5. Change automatic thread defaults only after repeated measurements.
