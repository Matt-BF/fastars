use super::block::{EncodedBlock, IndexRecord, encode_block};
use super::format::{DirectoryEntry, HEADER_SIZE, Header, write_directory, write_header};
use std::fs::{self, File};
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::mem;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread::{self, JoinHandle};

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
    encoder: Option<ParallelEncoder>,
    committed: bool,
}

struct EncodeJob {
    sequence: usize,
    records: Vec<IndexRecord>,
}

struct EncodeResult {
    sequence: usize,
    encoded: io::Result<EncodedBlock>,
    record_count: u64,
    last_id: String,
}

struct ParallelEncoder {
    senders: Vec<Sender<EncodeJob>>,
    results: Receiver<EncodeResult>,
    workers: Vec<JoinHandle<()>>,
    next_worker: usize,
    submitted: usize,
    written: usize,
    ready: Vec<EncodeResult>,
    max_in_flight: usize,
}

impl IndexWriter {
    pub(crate) fn with_threads(output_path: &str, threads: usize) -> io::Result<Self> {
        Self::with_limits_and_threads(output_path, TARGET_RECORDS, MAX_RAW_BLOCK_SIZE, threads)
    }

    #[cfg(test)]
    pub(super) fn with_limits(
        output_path: &str,
        target_records: u32,
        max_raw_block_size: usize,
    ) -> io::Result<Self> {
        Self::with_limits_and_threads(output_path, target_records, max_raw_block_size, 1)
    }

    pub(super) fn with_limits_and_threads(
        output_path: &str,
        target_records: u32,
        max_raw_block_size: usize,
        threads: usize,
    ) -> io::Result<Self> {
        if target_records == 0 || max_raw_block_size == 0 || threads == 0 {
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
            encoder: (threads > 1).then(|| ParallelEncoder::new(threads)),
            committed: false,
        })
    }

    pub(crate) fn add(&mut self, record: IndexRecord) -> io::Result<()> {
        if record.full_id.is_empty()
            || record
                .full_id
                .bytes()
                .any(|byte| byte.is_ascii_whitespace())
        {
            return invalid("Invalid FASTA ID in sorted index input");
        }
        let header_suffix = record.header_suffix.as_deref().unwrap_or("");
        if header_suffix
            .bytes()
            .next()
            .is_some_and(|byte| !byte.is_ascii_whitespace())
            || header_suffix.contains(['\n', '\r'])
        {
            return invalid("Invalid complete FASTA header in sorted index input");
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
        let estimated_size = record
            .full_id
            .len()
            .saturating_add(header_suffix.len())
            .saturating_add(64);
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
        while self
            .encoder
            .as_ref()
            .is_some_and(ParallelEncoder::has_pending)
        {
            self.write_next_encoded()?;
        }
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
        let records = mem::take(&mut self.pending);
        self.pending_size = 0;
        if let Some(encoder) = &mut self.encoder {
            encoder.submit(records)?;
            if encoder.in_flight() >= encoder.max_in_flight {
                self.write_next_encoded()?;
            }
            return Ok(());
        }

        let record_count = records.len() as u64;
        let last_id = records.last().unwrap().full_id.clone();
        let encoded = encode_block(&records)?;
        self.write_encoded(encoded, record_count, last_id)
    }

    fn write_next_encoded(&mut self) -> io::Result<()> {
        let result = self
            .encoder
            .as_mut()
            .ok_or_else(|| io::Error::other("Parallel encoder is not configured"))?
            .receive_next()?;
        self.write_encoded(result.encoded?, result.record_count, result.last_id)
    }

    fn write_encoded(
        &mut self,
        encoded: EncodedBlock,
        record_count: u64,
        last_id: String,
    ) -> io::Result<()> {
        let stored_len = encoded.stored.len() as u64;
        let block_offset = HEADER_SIZE
            .checked_add(self.blocks_len)
            .ok_or_else(|| io::Error::other(".ffx file is too large"))?;
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
        Ok(())
    }
}

impl ParallelEncoder {
    fn new(threads: usize) -> Self {
        let (result_sender, results) = mpsc::channel();
        let mut senders = Vec::with_capacity(threads);
        let mut workers = Vec::with_capacity(threads);
        for _ in 0..threads {
            let (sender, receiver) = mpsc::channel::<EncodeJob>();
            let result_sender = result_sender.clone();
            workers.push(thread::spawn(move || {
                while let Ok(job) = receiver.recv() {
                    let record_count = job.records.len() as u64;
                    let last_id = job.records.last().unwrap().full_id.clone();
                    let encoded = encode_block(&job.records);
                    if result_sender
                        .send(EncodeResult {
                            sequence: job.sequence,
                            encoded,
                            record_count,
                            last_id,
                        })
                        .is_err()
                    {
                        break;
                    }
                }
            }));
            senders.push(sender);
        }
        drop(result_sender);
        Self {
            senders,
            results,
            workers,
            next_worker: 0,
            submitted: 0,
            written: 0,
            ready: Vec::new(),
            max_in_flight: threads.saturating_mul(2),
        }
    }

    fn submit(&mut self, records: Vec<IndexRecord>) -> io::Result<()> {
        let sequence = self.submitted;
        self.senders[self.next_worker]
            .send(EncodeJob { sequence, records })
            .map_err(|_| io::Error::other("Parallel block encoder stopped unexpectedly"))?;
        self.next_worker = (self.next_worker + 1) % self.senders.len();
        self.submitted += 1;
        Ok(())
    }

    fn receive_next(&mut self) -> io::Result<EncodeResult> {
        loop {
            if let Some(position) = self
                .ready
                .iter()
                .position(|result| result.sequence == self.written)
            {
                self.written += 1;
                return Ok(self.ready.swap_remove(position));
            }
            let result = self
                .results
                .recv()
                .map_err(|_| io::Error::other("Parallel block encoder stopped unexpectedly"))?;
            if result.sequence == self.written {
                self.written += 1;
                return Ok(result);
            }
            self.ready.push(result);
        }
    }

    fn in_flight(&self) -> usize {
        self.submitted - self.written
    }

    fn has_pending(&self) -> bool {
        self.written < self.submitted
    }
}

impl Drop for ParallelEncoder {
    fn drop(&mut self) {
        self.senders.clear();
        for worker in self.workers.drain(..) {
            let _ = worker.join();
        }
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
