mod external_sort;
mod scanner;

use crate::fasta::FastaReader;
use crate::ffx::{IndexRecord, IndexStats, IndexWriter, read_u64};
use crate::progress::{self, ProgressBar, Unit};
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::mem;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

pub const DEFAULT_SORT_MEMORY_MIB: usize = 512;
static NEXT_TEMPORARY_FILE: AtomicU64 = AtomicU64::new(0);

pub fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct RecordMetadata {
    virtual_offset: u64,
    sequence_length: u64,
    line_bases: u64,
    line_width: u64,
}

pub fn build_index_from_fasta(
    fasta_path: &str,
    output_path: &str,
    temp_directory: Option<&str>,
    sort_memory_bytes: usize,
    threads: usize,
    show_progress: bool,
) -> Result<(), Box<dyn Error>> {
    let mut reader = FastaReader::with_threads(fasta_path, threads)?;
    let mut builder = IndexBuilder::new(
        output_path,
        temp_directory,
        sort_memory_bytes,
        threads,
        show_progress,
    )?;
    let mut scan_progress =
        ProgressBar::new(show_progress, "scan", reader.total_bytes(), Unit::Bytes);
    scanner::scan_fasta(
        &mut reader,
        &mut scan_progress,
        |full_id, header, metadata| builder.add_record(&full_id, &header, metadata),
    )?;
    scan_progress.finish(reader.consumed_bytes());

    drop(reader);
    builder.finish(output_path)?;
    Ok(())
}

pub fn build_index_from_fai_gzi(
    fai_path: &str,
    gzi_path: &str,
    output_path: &str,
    temp_directory: Option<&str>,
    sort_memory_bytes: usize,
    threads: usize,
    show_progress: bool,
) -> Result<(), Box<dyn Error>> {
    let mut gzi = File::open(gzi_path)?;
    let gzi_count = read_gzi_count(&mut gzi)?;
    let fai = File::open(fai_path)?;
    let fai_len = fai.metadata()?.len();
    let mut reader = BufReader::new(fai);
    let mut builder = IndexBuilder::new(
        output_path,
        temp_directory,
        sort_memory_bytes,
        threads,
        show_progress,
    )?;
    let mut scan_progress = ProgressBar::new(show_progress, "scan", fai_len, Unit::Bytes);
    let mut consumed = 0_u64;
    let mut line = String::new();

    loop {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        consumed = consumed
            .checked_add(bytes as u64)
            .ok_or_else(|| io::Error::other("FAI input is too large"))?;
        if line.trim().is_empty() {
            scan_progress.update(consumed);
            continue;
        }
        let line = line.trim_end_matches(['\r', '\n']);
        let (full_id, mut metadata) = parse_fai_record(line)?;
        metadata.virtual_offset = gzi_virtual_offset(&mut gzi, gzi_count, metadata.virtual_offset)?;
        builder.add_record(&full_id, &full_id, metadata)?;
        scan_progress.update(consumed);
    }
    scan_progress.finish(consumed);

    builder.finish(output_path)?;
    Ok(())
}

pub(super) fn parse_header(bytes: &[u8]) -> io::Result<(String, String)> {
    let mut end = bytes.len();
    if end > 0 && bytes[end - 1] == b'\n' {
        end -= 1;
    }
    if end > 0 && bytes[end - 1] == b'\r' {
        end -= 1;
    }

    let content = &bytes[1..end];
    let id_end = content
        .iter()
        .position(|byte| byte.is_ascii_whitespace())
        .unwrap_or(content.len());
    if id_end == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTA header has no ID",
        ));
    }

    let header = String::from_utf8(content.to_vec())
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Non-UTF-8 FASTA header"))?;
    Ok((header[..id_end].to_string(), header))
}

pub(super) fn sequence_line_layout(bytes: &[u8]) -> io::Result<(u64, u64)> {
    let mut base_len = bytes.len();
    if base_len > 0 && bytes[base_len - 1] == b'\n' {
        base_len -= 1;
    }
    if base_len > 0 && bytes[base_len - 1] == b'\r' {
        base_len -= 1;
    }
    if bytes[..base_len]
        .iter()
        .any(|byte| byte.is_ascii_whitespace())
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Whitespace inside FASTA sequence lines is not supported",
        ));
    }
    Ok((base_len as u64, bytes.len() as u64))
}

fn parse_fai_record(line: &str) -> io::Result<(String, RecordMetadata)> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 5 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Malformed .fai line",
        ));
    }
    Ok((
        fields[0].to_string(),
        RecordMetadata {
            virtual_offset: fields[2]
                .parse()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid .fai offset"))?,
            sequence_length: fields[1].parse().map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "Invalid .fai sequence length")
            })?,
            line_bases: fields[3].parse().map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "Invalid .fai line_bases")
            })?,
            line_width: fields[4].parse().map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "Invalid .fai line_width")
            })?,
        },
    ))
}

struct IndexBuilder {
    records: Vec<IndexRecord>,
    buffered_bytes: usize,
    sort_memory_bytes: usize,
    threads: usize,
    sort_input: Option<BufWriter<File>>,
    record_count: u64,
    last_id: Option<String>,
    input_is_sorted: bool,
    sort_input_path: PathBuf,
    show_progress: bool,
}

impl IndexBuilder {
    fn new(
        output_path: &str,
        temp_directory: Option<&str>,
        sort_memory_bytes: usize,
        threads: usize,
        show_progress: bool,
    ) -> io::Result<Self> {
        if sort_memory_bytes == 0 || threads == 0 {
            return Err(io::Error::other(
                "Sort memory and thread count must be greater than zero",
            ));
        }
        let output = Path::new(output_path);
        let temp_root = temp_directory
            .map(PathBuf::from)
            .or_else(|| output.parent().map(PathBuf::from))
            .filter(|path| !path.as_os_str().is_empty())
            .unwrap_or_else(|| PathBuf::from("."));
        fs::create_dir_all(&temp_root)?;

        let pid = std::process::id();
        let nonce = NEXT_TEMPORARY_FILE.fetch_add(1, Ordering::Relaxed);
        let sort_input_path = temp_root.join(format!("fastars-{pid}-{nonce}.sort-input.tmp"));

        Ok(Self {
            records: Vec::new(),
            buffered_bytes: 0,
            sort_memory_bytes,
            threads,
            sort_input: None,
            record_count: 0,
            last_id: None,
            input_is_sorted: true,
            sort_input_path,
            show_progress,
        })
    }

    fn add_record(
        &mut self,
        full_id: &str,
        header: &str,
        metadata: RecordMetadata,
    ) -> io::Result<()> {
        if full_id.bytes().any(|byte| byte.is_ascii_whitespace()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Whitespace inside the primary FASTA ID is not supported",
            ));
        }
        let header_suffix = header.strip_prefix(full_id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTA header does not start with its primary ID",
            )
        })?;
        if header_suffix
            .bytes()
            .next()
            .is_some_and(|byte| !byte.is_ascii_whitespace())
            || header_suffix.contains(['\n', '\r'])
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid complete FASTA header",
            ));
        }

        if self
            .last_id
            .as_deref()
            .is_some_and(|previous| previous > full_id)
        {
            self.input_is_sorted = false;
        }
        self.last_id = Some(full_id.to_string());

        let record = IndexRecord {
            full_id: full_id.to_string(),
            header_suffix: (!header_suffix.is_empty()).then(|| header_suffix.into()),
            virtual_offset: metadata.virtual_offset,
            sequence_length: metadata.sequence_length,
            line_bases: metadata.line_bases,
            line_width: metadata.line_width,
        };
        if let Some(sort_input) = &mut self.sort_input {
            external_sort::write_record(sort_input, &record)?;
        } else {
            let record_bytes = sort_memory_cost(&record);
            if self
                .buffered_bytes
                .checked_add(record_bytes)
                .is_some_and(|total| total <= self.sort_memory_bytes)
            {
                self.buffered_bytes += record_bytes;
                self.records.push(record);
            } else {
                self.spill_records(record)?;
            }
        }

        self.record_count = self
            .record_count
            .checked_add(1)
            .ok_or_else(|| io::Error::other("Too many FASTA records"))?;
        Ok(())
    }

    fn finish(mut self, output_path: &str) -> Result<(), Box<dyn Error>> {
        let stats = if let Some(mut sort_input) = self.sort_input.take() {
            sort_input.flush()?;
            drop(sort_input);
            self.write_spilled_index(output_path)?
        } else {
            if !self.input_is_sorted {
                let started = Instant::now();
                progress::status(
                    self.show_progress,
                    "sort",
                    &format!("sorting {} records in memory", self.record_count),
                );
                self.records
                    .sort_by(|left, right| left.full_id.cmp(&right.full_id));
                progress::status(
                    self.show_progress,
                    "sort",
                    &format!("done in {}", progress::format_duration(started.elapsed())),
                );
            } else {
                eprintln!("[INFO] input IDs are already sorted; skipping sort");
            }
            let records = mem::take(&mut self.records);
            write_record_index(
                records,
                output_path,
                self.record_count,
                self.threads,
                self.show_progress,
            )?
        };
        report_stats(&stats);
        Ok(())
    }

    fn spill_records(&mut self, record: IndexRecord) -> io::Result<()> {
        let mut sort_input = BufWriter::new(File::create(&self.sort_input_path)?);
        for buffered in self.records.drain(..) {
            external_sort::write_record(&mut sort_input, &buffered)?;
        }
        external_sort::write_record(&mut sort_input, &record)?;
        self.buffered_bytes = 0;
        self.sort_input = Some(sort_input);
        Ok(())
    }

    fn write_spilled_index(&self, output_path: &str) -> Result<IndexStats, Box<dyn Error>> {
        Ok(external_sort::write_index(
            &self.sort_input_path,
            output_path,
            self.record_count,
            self.sort_memory_bytes,
            self.threads,
            self.input_is_sorted,
            self.show_progress,
        )?)
    }
}

fn sort_memory_cost(record: &IndexRecord) -> usize {
    record
        .full_id
        .len()
        .saturating_add(
            record
                .header_suffix
                .as_deref()
                .map(str::len)
                .unwrap_or_default(),
        )
        .saturating_add(3 * mem::size_of::<IndexRecord>())
}

fn report_stats(stats: &IndexStats) {
    let ratio = if stats.raw_block_bytes == 0 {
        0.0
    } else {
        stats.stored_block_bytes as f64 / stats.raw_block_bytes as f64
    };
    eprintln!(
        "[INFO] indexed {} records in {} blocks; {} bytes on disk (blocks compressed to {:.1}%)",
        stats.record_count,
        stats.block_count,
        stats.file_size,
        ratio * 100.0
    );
}

impl Drop for IndexBuilder {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.sort_input_path);
    }
}

fn write_record_index(
    records: Vec<IndexRecord>,
    output_path: &str,
    expected_records: u64,
    threads: usize,
    show_progress: bool,
) -> io::Result<IndexStats> {
    if records.len() as u64 != expected_records {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Buffered index record count does not match scanned FASTA records",
        ));
    }
    let mut writer = IndexWriter::with_threads(output_path, threads)?;
    let mut progress = ProgressBar::new(show_progress, "write", expected_records, Unit::Records);
    let mut count = 0_u64;
    for record in records {
        writer.add(record)?;
        count += 1;
        progress.update(count);
    }
    let stats = writer.finish()?;
    progress.finish(count);
    Ok(stats)
}

fn read_gzi_count(file: &mut File) -> io::Result<u64> {
    file.seek(SeekFrom::Start(0))?;
    read_u64(file)
}

fn read_gzi_pair(file: &mut File, index: u64) -> io::Result<(u64, u64)> {
    file.seek(SeekFrom::Start(8 + index * 16))?;
    Ok((read_u64(file)?, read_u64(file)?))
}

fn gzi_virtual_offset(file: &mut File, entry_count: u64, sequence_offset: u64) -> io::Result<u64> {
    let mut low = 0;
    let mut high = entry_count;
    while low < high {
        let middle = low + (high - low) / 2;
        let (_, uncompressed) = read_gzi_pair(file, middle)?;
        if uncompressed <= sequence_offset {
            low = middle + 1;
        } else {
            high = middle;
        }
    }

    let (compressed, uncompressed) = if low == 0 {
        (0, 0)
    } else {
        read_gzi_pair(file, low - 1)?
    };
    let within_block = sequence_offset
        .checked_sub(uncompressed)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid .gzi offset"))?;
    if within_block >= 65_536 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Sequence offset exceeds BGZF block",
        ));
    }
    Ok((compressed << 16) | within_block)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ffx::FfxIndex;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_PATH: AtomicU64 = AtomicU64::new(0);

    fn test_path(suffix: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "fastars-indexer-{}-{}-{suffix}",
            std::process::id(),
            NEXT_PATH.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn fai_builder_sorts_and_indexes_records() {
        let fai = test_path("input.fai");
        let gzi = test_path("input.gzi");
        let output = test_path("output.ffx");
        fs::write(
            &fai,
            "zeta\t10\t30\t10\t11\nalpha\t20\t0\t20\t21\nduplicate\t5\t20\t5\t6\nduplicate\t5\t10\t5\t6\n",
        )
        .unwrap();
        fs::write(&gzi, 0_u64.to_le_bytes()).unwrap();

        build_index_from_fai_gzi(
            fai.to_str().unwrap(),
            gzi.to_str().unwrap(),
            output.to_str().unwrap(),
            None,
            1,
            2,
            false,
        )
        .unwrap();

        let mut index = FfxIndex::open(output.to_str().unwrap()).unwrap();
        assert_eq!(index.find_exact("duplicate").unwrap().len(), 2);
        assert_eq!(index.find_prefix("a").unwrap()[0].full_id, "alpha");
        assert_eq!(index.find_prefix("z").unwrap()[0].virtual_offset, 30);

        for path in [fai, gzi, output] {
            fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn plain_fasta_builder_indexes_byte_offsets() {
        let fasta = test_path("input.fna");
        let output = test_path("plain-output.ffx");
        fs::write(&fasta, b">zeta description\nACGT\nAC\n>alpha\nTTT\n").unwrap();

        build_index_from_fasta(
            fasta.to_str().unwrap(),
            output.to_str().unwrap(),
            None,
            usize::MAX,
            4,
            false,
        )
        .unwrap();

        let mut index = FfxIndex::open(output.to_str().unwrap()).unwrap();
        let alpha = index.find_exact("alpha").unwrap().pop().unwrap();
        let zeta = index.find_exact("zeta").unwrap().pop().unwrap();
        let zeta_by_header = index.find_exact("zeta description").unwrap().pop().unwrap();
        assert_eq!(alpha.virtual_offset, 33);
        assert_eq!(alpha.sequence_length, 3);
        assert_eq!(zeta.virtual_offset, 18);
        assert_eq!(zeta.header, "zeta description");
        assert_eq!(zeta_by_header.virtual_offset, zeta.virtual_offset);
        assert_eq!(zeta.sequence_length, 6);
        assert_eq!((zeta.line_bases, zeta.line_width), (4, 5));

        for path in [fasta, output] {
            fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn builder_detects_sorted_and_unsorted_ids() {
        let directory = test_path("ordering-directory");
        fs::create_dir(&directory).unwrap();
        let output = directory.join("output.ffx");
        let metadata = RecordMetadata {
            virtual_offset: 0,
            sequence_length: 1,
            line_bases: 1,
            line_width: 2,
        };
        let mut builder = IndexBuilder::new(
            output.to_str().unwrap(),
            directory.to_str(),
            usize::MAX,
            1,
            false,
        )
        .unwrap();

        builder.add_record("alpha", "alpha", metadata).unwrap();
        builder.add_record("alpha", "alpha", metadata).unwrap();
        builder.add_record("zeta", "zeta", metadata).unwrap();
        assert!(builder.input_is_sorted);

        builder.add_record("beta", "beta", metadata).unwrap();
        assert!(!builder.input_is_sorted);

        drop(builder);
        fs::remove_dir(directory).unwrap();
    }

    #[test]
    fn in_memory_sort_preserves_duplicate_input_order() {
        let output = test_path("memory-sort.ffx");
        let metadata = |virtual_offset| RecordMetadata {
            virtual_offset,
            sequence_length: 1,
            line_bases: 1,
            line_width: 2,
        };
        let mut builder =
            IndexBuilder::new(output.to_str().unwrap(), None, usize::MAX, 2, false).unwrap();
        builder.add_record("zeta", "zeta", metadata(30)).unwrap();
        builder
            .add_record("duplicate", "duplicate", metadata(20))
            .unwrap();
        builder
            .add_record("duplicate", "duplicate", metadata(10))
            .unwrap();
        builder.add_record("alpha", "alpha", metadata(0)).unwrap();
        builder.finish(output.to_str().unwrap()).unwrap();

        let mut index = FfxIndex::open(output.to_str().unwrap()).unwrap();
        let offsets = index
            .find_exact("duplicate")
            .unwrap()
            .into_iter()
            .map(|record| record.virtual_offset)
            .collect::<Vec<_>>();
        assert_eq!(offsets, vec![20, 10]);

        fs::remove_file(output).unwrap();
    }

    #[test]
    fn dropping_builder_removes_temporary_files() {
        let directory = test_path("temporary-directory");
        fs::create_dir(&directory).unwrap();
        let output = directory.join("output.ffx");
        let builder = IndexBuilder::new(
            output.to_str().unwrap(),
            directory.to_str(),
            usize::MAX,
            1,
            false,
        )
        .unwrap();
        let input = builder.sort_input_path.clone();
        drop(builder);
        assert!(!input.exists());
        fs::remove_dir(directory).unwrap();
    }

    #[test]
    fn external_sort_preserves_duplicate_input_order() {
        let output = test_path("external-sort.ffx");
        let metadata = |virtual_offset| RecordMetadata {
            virtual_offset,
            sequence_length: 1,
            line_bases: 1,
            line_width: 2,
        };
        let mut builder = IndexBuilder::new(output.to_str().unwrap(), None, 1, 2, false).unwrap();
        builder.add_record("zeta", "zeta", metadata(30)).unwrap();
        builder
            .add_record("duplicate", "duplicate", metadata(20))
            .unwrap();
        builder.add_record("alpha", "alpha", metadata(0)).unwrap();
        builder
            .add_record("duplicate", "duplicate", metadata(10))
            .unwrap();
        builder.finish(output.to_str().unwrap()).unwrap();

        let mut index = FfxIndex::open(output.to_str().unwrap()).unwrap();
        let offsets = index
            .find_exact("duplicate")
            .unwrap()
            .into_iter()
            .map(|record| record.virtual_offset)
            .collect::<Vec<_>>();
        assert_eq!(offsets, vec![20, 10]);
        fs::remove_file(output).unwrap();
    }
}
