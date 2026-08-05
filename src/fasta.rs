use crate::bgzf::BgzfReader;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

const PLAIN_CHUNK_SIZE: usize = 1024 * 1024;

pub(crate) struct FastaChunk {
    pub offset: u64,
    pub bytes: Vec<u8>,
}

pub(crate) enum FastaReader {
    Bgzf(BgzfReader),
    Plain(PlainFastaReader),
}

impl FastaReader {
    pub(crate) fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::with_threads(path, 1)
    }

    pub(crate) fn with_threads<P: AsRef<Path>>(path: P, threads: usize) -> io::Result<Self> {
        if threads == 0 {
            return Err(io::Error::other(
                "FASTA reader thread count must be greater than zero",
            ));
        }
        if is_gzip(&path)? {
            return Ok(Self::Bgzf(BgzfReader::with_threads(path, threads)?));
        }
        Ok(Self::Plain(PlainFastaReader::new(path)?))
    }

    pub(crate) fn read_chunk(&mut self) -> io::Result<Option<FastaChunk>> {
        match self {
            Self::Bgzf(reader) => reader.read_block().map(|block| {
                block.map(|block| FastaChunk {
                    offset: block.address << 16,
                    bytes: block.bytes,
                })
            }),
            Self::Plain(reader) => reader.read_chunk(),
        }
    }

    pub(crate) fn recycle_chunk(&mut self, chunk: FastaChunk) {
        match self {
            Self::Bgzf(reader) => reader.recycle_block(chunk.bytes),
            Self::Plain(reader) => reader.recycle_chunk(chunk.bytes),
        }
    }

    pub(crate) fn seek(&mut self, offset: u64) -> io::Result<()> {
        match self {
            Self::Bgzf(reader) => reader.seek_virtual(offset),
            Self::Plain(reader) => reader.seek(offset),
        }
    }

    pub(crate) fn read_sequence(
        &mut self,
        sequence_length: u64,
        line_bases: u64,
        line_width: u64,
    ) -> io::Result<Vec<u8>> {
        match self {
            Self::Bgzf(reader) => reader.read_sequence(sequence_length, line_bases, line_width),
            Self::Plain(reader) => reader.read_sequence(sequence_length, line_bases, line_width),
        }
    }
}

pub(crate) struct PlainFastaReader {
    reader: BufReader<File>,
    position: u64,
    spare_chunk: Vec<u8>,
}

impl PlainFastaReader {
    fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            reader: BufReader::new(File::open(path)?),
            position: 0,
            spare_chunk: Vec::new(),
        })
    }

    fn read_chunk(&mut self) -> io::Result<Option<FastaChunk>> {
        let offset = self.position;
        let mut bytes = std::mem::take(&mut self.spare_chunk);
        bytes.resize(PLAIN_CHUNK_SIZE, 0);
        let read = self.reader.read(&mut bytes)?;
        if read == 0 {
            return Ok(None);
        }
        bytes.truncate(read);
        self.position = self
            .position
            .checked_add(read as u64)
            .ok_or_else(|| io::Error::other("FASTA offset overflow"))?;
        Ok(Some(FastaChunk { offset, bytes }))
    }

    fn recycle_chunk(&mut self, bytes: Vec<u8>) {
        self.spare_chunk = bytes;
    }

    fn seek(&mut self, offset: u64) -> io::Result<()> {
        self.reader.seek(SeekFrom::Start(offset))?;
        self.position = offset;
        Ok(())
    }

    fn read_sequence(
        &mut self,
        sequence_length: u64,
        line_bases: u64,
        line_width: u64,
    ) -> io::Result<Vec<u8>> {
        if line_bases == 0 || line_width < line_bases {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid FASTA line layout in .ffx",
            ));
        }

        let capacity = usize::try_from(sequence_length)
            .map_err(|_| io::Error::other("FASTA sequence is too large"))?;
        let mut sequence = Vec::with_capacity(capacity);
        let mut remaining = sequence_length;
        while remaining > 0 {
            let bases = remaining.min(line_bases) as usize;
            let start = sequence.len();
            sequence.resize(start + bases, 0);
            self.reader.read_exact(&mut sequence[start..])?;
            remaining -= bases as u64;
            if remaining > 0 {
                let separator = usize::try_from(line_width - line_bases)
                    .map_err(|_| io::Error::other("FASTA line separator is too large"))?;
                let mut discarded = vec![0_u8; separator];
                self.reader.read_exact(&mut discarded)?;
            }
        }
        Ok(sequence)
    }
}

fn is_gzip<P: AsRef<Path>>(path: P) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut magic = [0_u8; 2];
    let read = file.read(&mut magic)?;
    Ok(read == magic.len() && magic == [31, 139])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_PATH: AtomicU64 = AtomicU64::new(0);

    fn test_path() -> PathBuf {
        std::env::temp_dir().join(format!(
            "fastars-plain-fasta-{}-{}.fna",
            std::process::id(),
            NEXT_PATH.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn plain_reader_tracks_offsets_and_fetches_sequence() {
        let path = test_path();
        fs::write(&path, b">alpha description\nACGT\nAC\n>beta\nTT\n").unwrap();

        let mut reader = FastaReader::with_threads(&path, 4).unwrap();
        let chunk = reader.read_chunk().unwrap().unwrap();
        assert_eq!(chunk.offset, 0);
        assert_eq!(chunk.bytes, b">alpha description\nACGT\nAC\n>beta\nTT\n");

        reader.seek(19).unwrap();
        assert_eq!(reader.read_sequence(6, 4, 5).unwrap(), b"ACGTAC");

        fs::remove_file(path).unwrap();
    }
}
