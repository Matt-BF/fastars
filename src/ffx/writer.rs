use super::block::{EncodedBlock, IndexRecord, encode_block};
use super::format::{DirectoryEntry, HEADER_SIZE, Header, write_directory, write_header};
use std::fs::{self, File};
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

const TARGET_RECORDS: u32 = 4_096;
const MAX_RAW_BLOCK_SIZE: usize = 4 * 1024 * 1024;

pub(crate) struct IndexStats {
    pub record_count: u64,
    pub block_count: u64,
    pub file_size: u64,
    pub raw_block_bytes: u64,
    pub stored_block_bytes: u64,
}

pub(crate) struct IndexWriter {
    output_path: PathBuf,
    temporary_path: PathBuf,
    output: BufWriter<File>,
    directory: Vec<DirectoryEntry>,
    pending: Vec<IndexRecord>,
    pending_size: usize,
    last_id: Option<String>,
    record_count: u64,
    blocks_len: u64,
    raw_block_bytes: u64,
    target_records: u32,
    max_raw_block_size: usize,
    committed: bool,
}

impl IndexWriter {
    pub(crate) fn new(output_path: &str) -> io::Result<Self> {
        Self::with_limits(output_path, TARGET_RECORDS, MAX_RAW_BLOCK_SIZE)
    }

    pub(super) fn with_limits(
        output_path: &str,
        target_records: u32,
        max_raw_block_size: usize,
    ) -> io::Result<Self> {
        if target_records == 0 || max_raw_block_size == 0 {
            return Err(io::Error::other("Invalid .ffx block configuration"));
        }
        let output_path = PathBuf::from(output_path);
        let parent = output_path
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        fs::create_dir_all(parent)?;
        let name = output_path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("fastars.ffx");
        let temporary_path = parent.join(format!(".{name}.{}.tmp", std::process::id()));
        let mut output = BufWriter::new(File::create(&temporary_path)?);
        output.write_all(&vec![0_u8; HEADER_SIZE as usize])?;

        Ok(Self {
            output_path,
            temporary_path,
            output,
            directory: Vec::new(),
            pending: Vec::new(),
            pending_size: 0,
            last_id: None,
            record_count: 0,
            blocks_len: 0,
            raw_block_bytes: 0,
            target_records,
            max_raw_block_size,
            committed: false,
        })
    }

    pub(crate) fn add(&mut self, record: IndexRecord) -> io::Result<()> {
        if record.full_id.is_empty() || record.full_id.contains(['\t', '\n', '\r']) {
            return invalid("Invalid FASTA ID in sorted index input");
        }
        if record.line_bases == 0 || record.line_width < record.line_bases {
            return invalid("Invalid FASTA line layout in sorted index input");
        }
        if let Some(last_id) = &self.last_id
            && last_id > &record.full_id
        {
            return invalid("Sorted index input is not lexically ordered");
        }

        // A plain-ID record uses at most ten bytes for each of its five
        // varints, so this is a conservative bound with room for framing.
        let estimated_size = record.full_id.len().saturating_add(64);
        if !self.pending.is_empty()
            && (self.pending.len() >= self.target_records as usize
                || self.pending_size.saturating_add(estimated_size) > self.max_raw_block_size)
        {
            self.flush_block()?;
        }
        self.pending_size = self.pending_size.saturating_add(estimated_size);
        self.last_id = Some(record.full_id.clone());
        self.pending.push(record);
        Ok(())
    }

    pub(crate) fn finish(mut self) -> io::Result<IndexStats> {
        self.flush_block()?;
        let directory_offset = HEADER_SIZE
            .checked_add(self.blocks_len)
            .ok_or_else(|| io::Error::other(".ffx file is too large"))?;
        let directory_len = write_directory(&mut self.output, &self.directory)?;
        let header = Header {
            record_count: self.record_count,
            block_count: self.directory.len() as u64,
            blocks_offset: HEADER_SIZE,
            blocks_len: self.blocks_len,
            directory_offset,
            directory_len,
            target_records: self.target_records,
        };
        self.output.flush()?;
        self.output.seek(SeekFrom::Start(0))?;
        write_header(&mut self.output, header)?;
        self.output.flush()?;
        self.output.get_ref().sync_all()?;

        let file_size = directory_offset
            .checked_add(directory_len)
            .ok_or_else(|| io::Error::other(".ffx file is too large"))?;
        fs::rename(&self.temporary_path, &self.output_path)?;
        self.committed = true;
        Ok(IndexStats {
            record_count: self.record_count,
            block_count: self.directory.len() as u64,
            file_size,
            raw_block_bytes: self.raw_block_bytes,
            stored_block_bytes: self.blocks_len,
        })
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let encoded = encode_block(&self.pending)?;
        let stored_len = encoded.stored.len() as u64;
        let block_offset = HEADER_SIZE
            .checked_add(self.blocks_len)
            .ok_or_else(|| io::Error::other(".ffx file is too large"))?;
        let record_count = self.pending.len() as u64;
        let last_id = self.pending.last().unwrap().full_id.clone();
        self.output.write_all(&encoded.stored)?;
        self.directory.push(directory_entry(
            block_offset,
            stored_len,
            self.record_count,
            record_count,
            last_id,
            &encoded,
        ));
        self.blocks_len = self
            .blocks_len
            .checked_add(stored_len)
            .ok_or_else(|| io::Error::other(".ffx file is too large"))?;
        self.raw_block_bytes = self
            .raw_block_bytes
            .checked_add(encoded.raw_len)
            .ok_or_else(|| io::Error::other(".ffx file is too large"))?;
        self.record_count = self
            .record_count
            .checked_add(record_count)
            .ok_or_else(|| io::Error::other("Too many FASTA records"))?;
        self.pending.clear();
        self.pending_size = 0;
        Ok(())
    }
}

impl Drop for IndexWriter {
    fn drop(&mut self) {
        if !self.committed {
            let _ = fs::remove_file(&self.temporary_path);
        }
    }
}

fn directory_entry(
    block_offset: u64,
    stored_len: u64,
    first_record_id: u64,
    record_count: u64,
    last_id: String,
    encoded: &EncodedBlock,
) -> DirectoryEntry {
    DirectoryEntry {
        block_offset,
        stored_len,
        raw_len: encoded.raw_len,
        first_record_id,
        record_count,
        flags: encoded.flags,
        last_id,
    }
}

fn invalid<T>(message: &str) -> io::Result<T> {
    Err(io::Error::new(io::ErrorKind::InvalidData, message))
}
