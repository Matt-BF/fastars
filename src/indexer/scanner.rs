use super::{RecordMetadata, parse_header_id, sequence_line_layout};
use crate::fasta::FastaReader;
use std::io;

pub(super) fn scan_fasta<F>(reader: &mut FastaReader, add_record: F) -> io::Result<()>
where
    F: FnMut(String, RecordMetadata) -> io::Result<()>,
{
    let mut scanner = FastaScanner::new(add_record);
    while let Some(chunk) = reader.read_chunk()? {
        scanner.push_chunk(chunk.offset, &chunk.bytes)?;
    }
    scanner.finish()
}

struct FastaScanner<F> {
    add_record: F,
    current: Option<PendingRecord>,
    partial_line: Vec<u8>,
    partial_offset: u64,
}

impl<F> FastaScanner<F>
where
    F: FnMut(String, RecordMetadata) -> io::Result<()>,
{
    fn new(add_record: F) -> Self {
        Self {
            add_record,
            current: None,
            partial_line: Vec::new(),
            partial_offset: 0,
        }
    }

    fn push_chunk(&mut self, chunk_offset: u64, bytes: &[u8]) -> io::Result<()> {
        let mut position = 0;
        if !self.partial_line.is_empty() {
            match memchr::memchr(b'\n', bytes) {
                Some(newline) => {
                    self.partial_line.extend_from_slice(&bytes[..=newline]);
                    let line = std::mem::take(&mut self.partial_line);
                    self.process_line(self.partial_offset, &line)?;
                    position = newline + 1;
                }
                None => {
                    self.partial_line.extend_from_slice(bytes);
                    return Ok(());
                }
            }
        }

        while position < bytes.len() {
            let remaining = &bytes[position..];
            if let Some(newline) = memchr::memchr(b'\n', remaining) {
                let end = position + newline + 1;
                self.process_line(offset_at(chunk_offset, position)?, &bytes[position..end])?;
                position = end;
            } else {
                self.partial_offset = offset_at(chunk_offset, position)?;
                self.partial_line.extend_from_slice(remaining);
                break;
            }
        }
        Ok(())
    }

    fn finish(mut self) -> io::Result<()> {
        if !self.partial_line.is_empty() {
            let line = std::mem::take(&mut self.partial_line);
            self.process_line(self.partial_offset, &line)?;
        }
        self.finish_record()
    }

    fn process_line(&mut self, offset: u64, bytes: &[u8]) -> io::Result<()> {
        if bytes.starts_with(b">") {
            self.finish_record()?;
            self.current = Some(PendingRecord::new(parse_header_id(bytes)?));
            return Ok(());
        }

        let record = self.current.as_mut().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Expected FASTA header line starting with '>'",
            )
        })?;
        record.add_sequence_line(offset, bytes)
    }

    fn finish_record(&mut self) -> io::Result<()> {
        let Some(record) = self.current.take() else {
            return Ok(());
        };
        let (id, metadata) = record.finish()?;
        (self.add_record)(id, metadata)
    }
}

struct PendingRecord {
    id: String,
    sequence_offset: Option<u64>,
    sequence_length: u64,
    line_bases: u64,
    line_width: u64,
    previous_layout: Option<(u64, u64)>,
}

impl PendingRecord {
    fn new(id: String) -> Self {
        Self {
            id,
            sequence_offset: None,
            sequence_length: 0,
            line_bases: 0,
            line_width: 0,
            previous_layout: None,
        }
    }

    fn add_sequence_line(&mut self, offset: u64, bytes: &[u8]) -> io::Result<()> {
        let (bases, width) = sequence_line_layout(bytes)?;
        if bases == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Empty FASTA sequence lines are not supported",
            ));
        }
        if self
            .previous_layout
            .is_some_and(|layout| layout != (self.line_bases, self.line_width))
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Variable-width FASTA wrapping is not supported",
            ));
        }
        if self.sequence_offset.is_none() {
            self.sequence_offset = Some(offset);
            self.line_bases = bases;
            self.line_width = width;
        }
        self.sequence_length = self
            .sequence_length
            .checked_add(bases)
            .ok_or_else(|| io::Error::other("FASTA sequence is too long"))?;
        self.previous_layout = Some((bases, width));
        Ok(())
    }

    fn finish(self) -> io::Result<(String, RecordMetadata)> {
        Ok((
            self.id,
            RecordMetadata {
                virtual_offset: self.sequence_offset.ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "FASTA record has no sequence")
                })?,
                sequence_length: self.sequence_length,
                line_bases: self.line_bases,
                line_width: self.line_width,
            },
        ))
    }
}

fn offset_at(chunk_offset: u64, position: usize) -> io::Result<u64> {
    chunk_offset
        .checked_add(position as u64)
        .ok_or_else(|| io::Error::other("FASTA offset overflow"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scanner_joins_lines_across_chunks_without_changing_offsets() {
        let mut records = Vec::new();
        let mut scanner = FastaScanner::new(|id, metadata| {
            records.push((id, metadata));
            Ok(())
        });
        scanner.push_chunk(0, b">alpha descr").unwrap();
        scanner.push_chunk(12, b"iption\nACG").unwrap();
        scanner.push_chunk(23, b"T\nAC\n>beta\nTT").unwrap();
        scanner.finish().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].0, "alpha");
        assert_eq!(records[0].1.virtual_offset, 19);
        assert_eq!(records[0].1.sequence_length, 6);
        assert_eq!((records[0].1.line_bases, records[0].1.line_width), (4, 5));
        assert_eq!(records[1].0, "beta");
        assert_eq!(records[1].1.virtual_offset, 34);
        assert_eq!(records[1].1.sequence_length, 2);
    }
}
