mod block;
mod format;
mod writer;

pub(crate) use block::IndexRecord;
use block::decode_block;
pub(crate) use format::read_u64;
use format::{DirectoryEntry, read_directory, read_header};
use regex::Regex;
use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
pub(crate) use writer::{IndexStats, IndexWriter};

const CACHE_BLOCKS: usize = 4;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FfxRecord {
    pub record_id: u64,
    /// The first whitespace-delimited token in the FASTA header.
    pub full_id: String,
    /// The complete FASTA header without the leading `>` or line ending.
    pub header: String,
    pub virtual_offset: u64,
    pub sequence_length: u64,
    pub line_bases: u64,
    pub line_width: u64,
}

pub struct FfxIndex {
    file: File,
    directory: Vec<DirectoryEntry>,
    cache: HashMap<usize, Vec<FfxRecord>>,
    cache_order: VecDeque<usize>,
}

impl FfxIndex {
    pub fn open(path: &str) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let file_len = file.metadata()?.len();
        let header = read_header(&mut file)?;
        let directory = read_directory(&mut file, header, file_len)?;
        Ok(Self {
            file,
            directory,
            cache: HashMap::new(),
            cache_order: VecDeque::new(),
        })
    }

    pub fn find_exact(&mut self, query: &str) -> io::Result<Vec<FfxRecord>> {
        let id = primary_id(query);
        let mut matches = self.find_exact_id(id)?;
        if query.len() != id.len() {
            matches.retain(|record| record.header == query);
        }
        Ok(matches)
    }

    fn find_exact_id(&mut self, id: &str) -> io::Result<Vec<FfxRecord>> {
        let mut matches = Vec::new();
        let mut block_index = self.lower_bound_block(id);
        while block_index < self.directory.len() {
            let records = self.read_block(block_index)?;
            let start = records.partition_point(|record| record.full_id.as_str() < id);
            for record in records.into_iter().skip(start) {
                match record.full_id.as_str().cmp(id) {
                    Ordering::Less => continue,
                    Ordering::Equal => matches.push(record),
                    Ordering::Greater => return Ok(matches),
                }
            }
            block_index += 1;
        }
        Ok(matches)
    }

    pub fn find_prefix(&mut self, prefix: &str) -> io::Result<Vec<FfxRecord>> {
        let id_prefix = primary_id(prefix);
        if prefix.len() != id_prefix.len() {
            let mut matches = self.find_exact_id(id_prefix)?;
            matches.retain(|record| record.header.starts_with(prefix));
            return Ok(matches);
        }

        let mut matches = Vec::new();
        let mut block_index = self.lower_bound_block(id_prefix);
        while block_index < self.directory.len() {
            let records = self.read_block(block_index)?;
            let start = records.partition_point(|record| record.full_id.as_str() < id_prefix);
            for record in records.into_iter().skip(start) {
                if record.full_id.starts_with(id_prefix) {
                    matches.push(record);
                } else {
                    return Ok(matches);
                }
            }
            block_index += 1;
        }
        Ok(matches)
    }

    pub fn find_regex(&mut self, regex: &Regex, invert: bool) -> io::Result<Vec<FfxRecord>> {
        let mut records = Vec::new();
        self.for_each_regex(regex, invert, |record| {
            records.push(record);
            Ok(())
        })?;
        Ok(records)
    }

    pub fn for_each_regex<F>(&mut self, regex: &Regex, invert: bool, mut visit: F) -> io::Result<()>
    where
        F: FnMut(FfxRecord) -> io::Result<()>,
    {
        for block_index in 0..self.directory.len() {
            for record in self.read_block(block_index)? {
                if regex.is_match(&record.header) != invert {
                    visit(record)?;
                }
            }
        }
        Ok(())
    }

    fn lower_bound_block(&self, target: &str) -> usize {
        self.directory
            .partition_point(|entry| entry.last_id.as_str() < target)
    }

    fn read_block(&mut self, block_index: usize) -> io::Result<Vec<FfxRecord>> {
        if let Some(records) = self.cache.get(&block_index) {
            return Ok(records.clone());
        }
        let entry = self.directory.get(block_index).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "Block index is out of bounds")
        })?;
        let stored_len = usize::try_from(entry.stored_len)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Block is too large"))?;
        self.file.seek(SeekFrom::Start(entry.block_offset))?;
        let mut stored = vec![0_u8; stored_len];
        self.file.read_exact(&mut stored)?;
        let records = decode_block(&stored, entry)?;

        if self.cache.len() == CACHE_BLOCKS
            && let Some(expired) = self.cache_order.pop_front()
        {
            self.cache.remove(&expired);
        }
        self.cache.insert(block_index, records.clone());
        self.cache_order.push_back(block_index);
        Ok(records)
    }
}

fn primary_id(header: &str) -> &str {
    let end = header
        .bytes()
        .position(|byte| byte.is_ascii_whitespace())
        .unwrap_or(header.len());
    &header[..end]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

    static NEXT_PATH: AtomicU64 = AtomicU64::new(0);

    fn path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "fastars-{name}-{}-{}.ffx",
            std::process::id(),
            NEXT_PATH.fetch_add(1, AtomicOrdering::Relaxed)
        ))
    }

    fn record(id: &str, offset: u64) -> IndexRecord {
        IndexRecord {
            full_id: id.to_string(),
            header_suffix: None,
            virtual_offset: offset,
            sequence_length: u32::MAX as u64 + 10,
            line_bases: 1_000_000,
            line_width: 1_000_002,
        }
    }

    fn described_record(id: &str, header: &str, offset: u64) -> IndexRecord {
        let mut record = record(id, offset);
        record.header_suffix = Some(header.strip_prefix(id).unwrap().into());
        record
    }

    #[test]
    fn queries_work_across_small_block_boundaries() {
        let output = path("queries");
        let mut writer = IndexWriter::with_limits(output.to_str().unwrap(), 2, 1024).unwrap();
        for value in [
            record("a", 10),
            described_record("duplicate", "duplicate first description", 20),
            described_record("duplicate", "duplicate second description", 30),
            record("duplicate", 40),
            record("prefix-1", 50),
            record("prefix-2", 60),
            record("z", 70),
        ] {
            writer.add(value).unwrap();
        }
        writer.finish().unwrap();

        let mut index = FfxIndex::open(output.to_str().unwrap()).unwrap();
        assert_eq!(index.find_exact("duplicate").unwrap().len(), 3);
        assert_eq!(
            index.find_exact("duplicate second description").unwrap()[0].virtual_offset,
            30
        );
        assert_eq!(index.find_prefix("prefix-").unwrap().len(), 2);
        assert_eq!(
            index.find_prefix("duplicate first").unwrap()[0].header,
            "duplicate first description"
        );
        assert!(index.find_exact("missing").unwrap().is_empty());
        let regex = Regex::new("^(a|z)$").unwrap();
        assert_eq!(index.find_regex(&regex, false).unwrap().len(), 2);
        assert_eq!(index.find_regex(&regex, true).unwrap().len(), 5);
        let description = Regex::new("second description$").unwrap();
        assert_eq!(index.find_regex(&description, false).unwrap().len(), 1);
        fs::remove_file(output).unwrap();
    }

    #[test]
    fn empty_index_round_trips() {
        let output = path("empty");
        IndexWriter::with_limits(output.to_str().unwrap(), 2, 1024)
            .unwrap()
            .finish()
            .unwrap();
        let mut index = FfxIndex::open(output.to_str().unwrap()).unwrap();
        assert!(index.find_exact("anything").unwrap().is_empty());
        fs::remove_file(output).unwrap();
    }

    #[test]
    fn writer_rejects_unsorted_records() {
        let output = path("unsorted");
        let mut writer = IndexWriter::with_limits(output.to_str().unwrap(), 2, 1024).unwrap();
        writer.add(record("z", 1)).unwrap();
        assert!(writer.add(record("a", 2)).is_err());
        drop(writer);
        assert!(!output.exists());
    }

    #[test]
    fn parallel_writer_matches_serial_output() {
        let serial_output = path("serial");
        let parallel_output = path("parallel");
        let records = (0..20)
            .map(|value| record(&format!("id-{value:03}"), value))
            .collect::<Vec<_>>();

        let mut serial =
            IndexWriter::with_limits(serial_output.to_str().unwrap(), 2, 1024).unwrap();
        let mut parallel =
            IndexWriter::with_limits_and_threads(parallel_output.to_str().unwrap(), 2, 1024, 4)
                .unwrap();
        for value in records {
            serial.add(value.clone()).unwrap();
            parallel.add(value).unwrap();
        }
        serial.finish().unwrap();
        parallel.finish().unwrap();

        assert_eq!(
            fs::read(&serial_output).unwrap(),
            fs::read(&parallel_output).unwrap()
        );
        fs::remove_file(serial_output).unwrap();
        fs::remove_file(parallel_output).unwrap();
    }

    #[test]
    fn very_long_id_can_exceed_the_normal_block_target() {
        let output = path("long-id");
        let id = "x".repeat(8_192);
        let mut writer = IndexWriter::with_limits(output.to_str().unwrap(), 2, 128).unwrap();
        writer.add(record(&id, 1)).unwrap();
        writer.finish().unwrap();
        let mut index = FfxIndex::open(output.to_str().unwrap()).unwrap();
        assert_eq!(index.find_exact(&id).unwrap()[0].full_id, id);
        fs::remove_file(output).unwrap();
    }

    #[test]
    fn truncated_directory_is_rejected() {
        let output = path("truncated");
        let mut writer = IndexWriter::with_limits(output.to_str().unwrap(), 2, 1024).unwrap();
        writer.add(record("id", 1)).unwrap();
        let stats = writer.finish().unwrap();
        File::options()
            .write(true)
            .open(&output)
            .unwrap()
            .set_len(stats.file_size - 1)
            .unwrap();
        assert!(FfxIndex::open(output.to_str().unwrap()).is_err());
        fs::remove_file(output).unwrap();
    }
}
