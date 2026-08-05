use super::sort_memory_cost;
use crate::ffx::{IndexRecord, IndexStats, IndexWriter};
use crate::progress::{ProgressBar, Unit};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const MAX_MERGE_RUNS: usize = 64;
const HAS_HEADER_SUFFIX: u64 = 1 << 63;

pub(super) fn write_record<W: Write>(writer: &mut W, record: &IndexRecord) -> io::Result<()> {
    let mut id_length = u64::try_from(record.full_id.len())
        .map_err(|_| io::Error::other("FASTA ID is too long"))?;
    let header_suffix = record.header_suffix.as_deref().unwrap_or("");
    if id_length >= HAS_HEADER_SUFFIX {
        return Err(io::Error::other("FASTA ID is too long"));
    }
    let suffix_length = u32::try_from(header_suffix.len())
        .map_err(|_| io::Error::other("FASTA header is too long"))?;
    if record.header_suffix.is_some() {
        id_length |= HAS_HEADER_SUFFIX;
    }
    writer.write_all(&id_length.to_le_bytes())?;
    writer.write_all(record.full_id.as_bytes())?;
    if record.header_suffix.is_some() {
        writer.write_all(&suffix_length.to_le_bytes())?;
        writer.write_all(header_suffix.as_bytes())?;
    }
    for value in [
        record.virtual_offset,
        record.sequence_length,
        record.line_bases,
        record.line_width,
    ] {
        writer.write_all(&value.to_le_bytes())?;
    }
    Ok(())
}

pub(super) fn write_index(
    input_path: &Path,
    output_path: &str,
    expected_records: u64,
    memory_bytes: usize,
    threads: usize,
    input_is_sorted: bool,
    show_progress: bool,
) -> io::Result<IndexStats> {
    if input_is_sorted {
        eprintln!("[INFO] input IDs are already sorted; skipping external sort");
        return write_binary_index(
            input_path,
            output_path,
            expected_records,
            threads,
            show_progress,
        );
    }

    let mut temporary = TemporaryRuns::new(input_path);
    let mut runs = create_runs(
        input_path,
        memory_bytes,
        expected_records,
        &mut temporary,
        show_progress,
    )?;
    let mut pass = 1;
    while runs.len() > MAX_MERGE_RUNS {
        let phase = format!("merge pass {pass}");
        let mut progress = ProgressBar::new(show_progress, phase, runs.len() as u64, Unit::Runs);
        let mut merged = Vec::new();
        let mut consumed_runs = 0_u64;
        for group in runs.chunks(MAX_MERGE_RUNS) {
            let path = temporary.next_path();
            let mut writer = BufWriter::new(File::create(&path)?);
            merge_runs(group, |record| write_record(&mut writer, &record))?;
            writer.flush()?;
            merged.push(path);
            consumed_runs += group.len() as u64;
            progress.update(consumed_runs);
        }
        progress.finish(consumed_runs);
        remove_runs(&runs);
        runs = merged;
        pass += 1;
    }

    let mut writer = IndexWriter::with_threads(output_path, threads)?;
    let mut progress = ProgressBar::new(show_progress, "write", expected_records, Unit::Records);
    let mut count = 0_u64;
    merge_runs(&runs, |record| {
        writer.add(record)?;
        count = count
            .checked_add(1)
            .ok_or_else(|| io::Error::other("Too many sorted records"))?;
        progress.update(count);
        Ok(())
    })?;
    validate_count(count, expected_records)?;
    let stats = writer.finish()?;
    progress.finish(count);
    Ok(stats)
}

fn create_runs(
    input_path: &Path,
    memory_bytes: usize,
    expected_records: u64,
    temporary: &mut TemporaryRuns,
    show_progress: bool,
) -> io::Result<Vec<PathBuf>> {
    let mut input = BufReader::new(File::open(input_path)?);
    let mut pending = None;
    let mut runs = Vec::new();
    let mut count = 0_u64;
    let mut progress = ProgressBar::new(show_progress, "sort", expected_records, Unit::Records);

    loop {
        let mut records = Vec::new();
        let mut buffered_bytes = 0_usize;
        if let Some(record) = pending.take() {
            buffered_bytes = sort_memory_cost(&record);
            records.push(record);
        }

        while let Some(record) = read_record(&mut input)? {
            count = count
                .checked_add(1)
                .ok_or_else(|| io::Error::other("Too many records in sort input"))?;
            progress.update(count);
            let record_bytes = sort_memory_cost(&record);
            if !records.is_empty()
                && buffered_bytes
                    .checked_add(record_bytes)
                    .is_none_or(|total| total > memory_bytes)
            {
                pending = Some(record);
                break;
            }
            buffered_bytes = buffered_bytes.saturating_add(record_bytes);
            records.push(record);
        }

        if records.is_empty() {
            break;
        }
        records.sort_by(|left, right| left.full_id.cmp(&right.full_id));
        let path = temporary.next_path();
        let mut writer = BufWriter::new(File::create(&path)?);
        for record in records {
            write_record(&mut writer, &record)?;
        }
        writer.flush()?;
        runs.push(path);
    }

    validate_count(count, expected_records)?;
    progress.finish(count);
    Ok(runs)
}

fn write_binary_index(
    input_path: &Path,
    output_path: &str,
    expected_records: u64,
    threads: usize,
    show_progress: bool,
) -> io::Result<IndexStats> {
    let mut input = BufReader::new(File::open(input_path)?);
    let mut writer = IndexWriter::with_threads(output_path, threads)?;
    let mut progress = ProgressBar::new(show_progress, "write", expected_records, Unit::Records);
    let mut count = 0_u64;
    while let Some(record) = read_record(&mut input)? {
        writer.add(record)?;
        count = count
            .checked_add(1)
            .ok_or_else(|| io::Error::other("Too many sorted records"))?;
        progress.update(count);
    }
    validate_count(count, expected_records)?;
    let stats = writer.finish()?;
    progress.finish(count);
    Ok(stats)
}

fn merge_runs<F>(paths: &[PathBuf], mut visit: F) -> io::Result<()>
where
    F: FnMut(IndexRecord) -> io::Result<()>,
{
    let mut readers = paths
        .iter()
        .map(|path| File::open(path).map(BufReader::new))
        .collect::<io::Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::new();
    for (run, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = read_record(reader)? {
            heap.push(RunHead { run, record });
        }
    }

    while let Some(head) = heap.pop() {
        let run = head.run;
        visit(head.record)?;
        if let Some(record) = read_record(&mut readers[run])? {
            heap.push(RunHead { run, record });
        }
    }
    Ok(())
}

fn read_record<R: Read>(reader: &mut R) -> io::Result<Option<IndexRecord>> {
    let Some(encoded_id_length) = read_optional_u64(reader)? else {
        return Ok(None);
    };
    let has_header_suffix = encoded_id_length & HAS_HEADER_SUFFIX != 0;
    let id_length = encoded_id_length & !HAS_HEADER_SUFFIX;
    let id_length = usize::try_from(id_length).map_err(|_| {
        io::Error::new(io::ErrorKind::InvalidData, "Temporary FASTA ID is too long")
    })?;
    if id_length == 0 {
        return invalid("Temporary sort record has an empty FASTA ID");
    }
    let mut id = vec![0_u8; id_length];
    reader.read_exact(&mut id)?;
    let full_id = String::from_utf8(id).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "Temporary sort record has a non-UTF-8 FASTA ID",
        )
    })?;
    let header_suffix = if has_header_suffix {
        let suffix_length = read_u32(reader)? as usize;
        let mut suffix = vec![0_u8; suffix_length];
        reader.read_exact(&mut suffix)?;
        let suffix = String::from_utf8(suffix).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Temporary sort record has a non-UTF-8 FASTA header",
            )
        })?;
        Some(suffix.into_boxed_str())
    } else {
        None
    };
    Ok(Some(IndexRecord {
        full_id,
        header_suffix,
        virtual_offset: read_u64(reader)?,
        sequence_length: read_u64(reader)?,
        line_bases: read_u64(reader)?,
        line_width: read_u64(reader)?,
    }))
}

fn read_optional_u64<R: Read>(reader: &mut R) -> io::Result<Option<u64>> {
    let mut bytes = [0_u8; 8];
    let read = reader.read(&mut bytes)?;
    if read == 0 {
        return Ok(None);
    }
    reader.read_exact(&mut bytes[read..])?;
    Ok(Some(u64::from_le_bytes(bytes)))
}

fn read_u64<R: Read>(reader: &mut R) -> io::Result<u64> {
    let mut bytes = [0_u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn read_u32<R: Read>(reader: &mut R) -> io::Result<u32> {
    let mut bytes = [0_u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn validate_count(actual: u64, expected: u64) -> io::Result<()> {
    if actual != expected {
        return invalid("Sorted index record count does not match scanned FASTA records");
    }
    Ok(())
}

fn remove_runs(paths: &[PathBuf]) {
    for path in paths {
        let _ = fs::remove_file(path);
    }
}

struct RunHead {
    run: usize,
    record: IndexRecord,
}

impl Eq for RunHead {}

impl PartialEq for RunHead {
    fn eq(&self, other: &Self) -> bool {
        self.run == other.run && self.record.full_id == other.record.full_id
    }
}

impl Ord for RunHead {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .record
            .full_id
            .cmp(&self.record.full_id)
            .then_with(|| other.run.cmp(&self.run))
    }
}

impl PartialOrd for RunHead {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

struct TemporaryRuns {
    input_path: PathBuf,
    paths: Vec<PathBuf>,
}

impl TemporaryRuns {
    fn new(input_path: &Path) -> Self {
        Self {
            input_path: input_path.to_path_buf(),
            paths: Vec::new(),
        }
    }

    fn next_path(&mut self) -> PathBuf {
        let mut name = self
            .input_path
            .file_name()
            .unwrap_or_default()
            .to_os_string();
        name.push(format!(".run-{}", self.paths.len()));
        let path = self.input_path.with_file_name(name);
        self.paths.push(path.clone());
        path
    }
}

impl Drop for TemporaryRuns {
    fn drop(&mut self) {
        remove_runs(&self.paths);
    }
}

fn invalid<T>(message: &str) -> io::Result<T> {
    Err(io::Error::new(io::ErrorKind::InvalidData, message))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ffx::FfxIndex;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_PATH: AtomicU64 = AtomicU64::new(0);

    fn test_path(suffix: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "fastars-external-sort-{}-{}-{suffix}",
            std::process::id(),
            NEXT_PATH.fetch_add(1, Ordering::Relaxed)
        ))
    }

    fn record(id: &str, offset: u64) -> IndexRecord {
        IndexRecord {
            full_id: id.to_string(),
            header_suffix: None,
            virtual_offset: offset,
            sequence_length: 10,
            line_bases: 10,
            line_width: 11,
        }
    }

    #[test]
    fn binary_records_round_trip() {
        let mut expected = record("alpha", 42);
        expected.header_suffix = Some(" description with spaces".into());
        let mut bytes = Vec::new();
        write_record(&mut bytes, &expected).unwrap();
        let mut input = bytes.as_slice();
        assert_eq!(read_record(&mut input).unwrap(), Some(expected));
        assert_eq!(read_record(&mut input).unwrap(), None);
    }

    #[test]
    fn truncated_binary_record_is_rejected() {
        let mut bytes = Vec::new();
        write_record(&mut bytes, &record("alpha", 42)).unwrap();
        bytes.pop();
        assert!(read_record(&mut bytes.as_slice()).is_err());
    }

    #[test]
    fn merge_is_stable_across_runs() {
        let paths = [test_path("merge-0"), test_path("merge-1")];
        for (path, records) in [
            (&paths[0], vec![record("duplicate", 10), record("zeta", 30)]),
            (&paths[1], vec![record("alpha", 0), record("duplicate", 20)]),
        ] {
            let mut writer = BufWriter::new(File::create(path).unwrap());
            for record in records {
                write_record(&mut writer, &record).unwrap();
            }
            writer.flush().unwrap();
        }

        let mut merged = Vec::new();
        merge_runs(&paths, |record| {
            merged.push((record.full_id, record.virtual_offset));
            Ok(())
        })
        .unwrap();
        assert_eq!(
            merged,
            vec![
                ("alpha".to_string(), 0),
                ("duplicate".to_string(), 10),
                ("duplicate".to_string(), 20),
                ("zeta".to_string(), 30),
            ]
        );

        remove_runs(&paths);
    }

    #[test]
    fn external_sort_merges_more_than_one_batch() {
        let input = test_path("input.tmp");
        let output = test_path("output.ffx");
        let mut writer = BufWriter::new(File::create(&input).unwrap());
        for value in (0..70).rev() {
            write_record(&mut writer, &record(&format!("id-{value:03}"), value)).unwrap();
        }
        writer.flush().unwrap();

        write_index(&input, output.to_str().unwrap(), 70, 1, 2, false, false).unwrap();
        let mut index = FfxIndex::open(output.to_str().unwrap()).unwrap();
        assert_eq!(index.find_exact("id-000").unwrap()[0].virtual_offset, 0);
        assert_eq!(index.find_exact("id-069").unwrap()[0].virtual_offset, 69);
        let run_prefix = format!(
            "{}.run-",
            input.file_name().unwrap_or_default().to_string_lossy()
        );
        assert!(
            !fs::read_dir(input.parent().unwrap())
                .unwrap()
                .filter_map(Result::ok)
                .any(|entry| entry.file_name().to_string_lossy().starts_with(&run_prefix))
        );

        fs::remove_file(input).unwrap();
        fs::remove_file(output).unwrap();
    }
}
