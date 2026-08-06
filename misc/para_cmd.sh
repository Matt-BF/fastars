#!/bin/bash
# need blastcmd, gnu parallel, and bgzip
set -e

pixi shell -e dev --manifest-path=/home/neri/Documents/GitHub/rolypoly/pyproject.toml
cd /run/media/neri/ssd2/ncbi_nr

# 1. Configuration
THREADS=12
DB="/home/neri/Documents/GitHub/rolypoly/data/reference_seqs/ncbi_ribovirus/protein_taxdb/nr_blastdb/nr"
TMP_DIR="nr_chunks_tmp"
OUT_FILE=/run/media/neri/ssd2/ncbi_nr/nr.fasta.bgz

# Verify the folder actually contains volume indexes first
if ! ls "${DB_DIR}"/nr.*.pin >/dev/null 2>&1; then
    echo "Error: No volume index files (.pin) found in $DB_DIR"
    exit 1
fi

# 2. Create isolated working directory.
mkdir -p "$TMP_DIR"

# 3. Direct File Stream Loop
# We find the files and immediately strip out the extension, then use a while loop to pipe them directly into parallel via standard input (stdin).
# max I/O read rate I get from this ssd is ~100mbs.
find "$DB_DIR" -maxdepth 1 -name "nr.*.pv in" | sed 's/\.pin$//' | sort | \
parallel -j "$THREADS" \
  "echo 'Processing volume: {/}' && blastdbcmd -db {} -dbtype prot -entry all -outfmt %f 2>$TMP_DIR/{/}.err | bgzip -@ 1 > $TMP_DIR/{/}.fasta.gz"

# 4. Concatenate binary chunks (Bypasses re-compression, instant append, might not be ideal compression wise but should still be managable)
cat "$TMP_DIR"/nr.*.fasta.gz > "$OUT_FILE"

echo "list of file:\n>ls -lsh \n" > bla.log
ls -lsh >> bla.log
du ./* -sh >>

# Clean up working directory
# rm -rf "$TMP_DIR"

