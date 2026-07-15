use regex::Regex;
use std::cmp::Ordering;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, Write};

const MAGIC: &[u8; 8] = b"FASTARS1";
const VERSION: u64 = 1;
pub(crate) const HEADER_SIZE: u64 = 64;
pub(crate) const RECORD_SIZE: u64 = 64;

#[derive(Clone)]
pub struct FfxRecord {
    pub record_id: u64,
    pub full_id: String,
    pub virtual_offset: u64,
    pub sequence_length: u64,
    pub line_bases: u64,
    pub line_width: u64,
}

#[derive(Clone, Copy)]
pub(crate) struct Header {
    pub record_count: u64,
    pub records_offset: u64,
    pub strings_offset: u64,
    pub strings_len: u64,
    pub order_offset: u64,
    pub order_len: u64,
}

#[derive(Clone, Copy)]
pub(crate) struct RecordEntry {
    pub full_id_offset: u64,
    pub full_id_len: u64,
    pub virtual_offset: u64,
    pub sequence_length: u64,
    pub line_bases: u64,
    pub line_width: u64,
}

pub struct FfxIndex {
    file: File,
    header: Header,
}

impl FfxIndex {
    pub fn open(path: &str) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let header = read_header(&mut file)?;
        Ok(Self { file, header })
    }

    pub fn find_exact(&mut self, id: &str) -> io::Result<Vec<FfxRecord>> {
        let mut position = self.lower_bound(id)?;
        let mut records = Vec::new();
        while position < self.header.record_count {
            let record_id = self.read_ordered_record_id(position)?;
            let record = self.read_record(record_id)?;
            if record.full_id != id {
                break;
            }
            records.push(record);
            position += 1;
        }
        Ok(records)
    }

    pub fn find_prefix(&mut self, prefix: &str) -> io::Result<Vec<FfxRecord>> {
        let mut position = self.lower_bound(prefix)?;
        let mut records = Vec::new();
        while position < self.header.record_count {
            let record_id = self.read_ordered_record_id(position)?;
            let record = self.read_record(record_id)?;
            if !record.full_id.starts_with(prefix) {
                break;
            }
            records.push(record);
            position += 1;
        }
        Ok(records)
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
        for position in 0..self.header.record_count {
            let record_id = self.read_ordered_record_id(position)?;
            let record = self.read_record(record_id)?;
            if regex.is_match(&record.full_id) != invert {
                visit(record)?;
            }
        }
        Ok(())
    }

    fn lower_bound(&mut self, target: &str) -> io::Result<u64> {
        let mut low = 0;
        let mut high = self.header.record_count;
        while low < high {
            let middle = low + (high - low) / 2;
            let record_id = self.read_ordered_record_id(middle)?;
            let full_id = self.read_full_id(record_id)?;
            match full_id.as_str().cmp(target) {
                Ordering::Less => low = middle + 1,
                Ordering::Equal | Ordering::Greater => high = middle,
            }
        }
        Ok(low)
    }

    fn read_ordered_record_id(&mut self, position: u64) -> io::Result<u64> {
        if position >= self.header.record_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Ordered index position is out of bounds",
            ));
        }
        self.file
            .seek(SeekFrom::Start(self.header.order_offset + position * 8))?;
        read_u64(&mut self.file)
    }

    fn read_record(&mut self, record_id: u64) -> io::Result<FfxRecord> {
        let entry = self.read_record_entry(record_id)?;
        let full_id = self.read_full_id_from_entry(entry)?;
        Ok(FfxRecord {
            record_id,
            full_id,
            virtual_offset: entry.virtual_offset,
            sequence_length: entry.sequence_length,
            line_bases: entry.line_bases,
            line_width: entry.line_width,
        })
    }

    fn read_full_id(&mut self, record_id: u64) -> io::Result<String> {
        let entry = self.read_record_entry(record_id)?;
        self.read_full_id_from_entry(entry)
    }

    fn read_record_entry(&mut self, record_id: u64) -> io::Result<RecordEntry> {
        if record_id >= self.header.record_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Record ID is out of bounds",
            ));
        }
        self.file.seek(SeekFrom::Start(
            self.header.records_offset + record_id * RECORD_SIZE,
        ))?;
        read_record_entry(&mut self.file)
    }

    fn read_full_id_from_entry(&mut self, entry: RecordEntry) -> io::Result<String> {
        let end = entry
            .full_id_offset
            .checked_add(entry.full_id_len)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid string offset"))?;
        if end > self.header.strings_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "String offset is out of bounds",
            ));
        }

        self.file.seek(SeekFrom::Start(
            self.header.strings_offset + entry.full_id_offset,
        ))?;
        let mut bytes = vec![0_u8; entry.full_id_len as usize];
        self.file.read_exact(&mut bytes)?;
        String::from_utf8(bytes)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Non-UTF-8 FASTA ID in .ffx"))
    }
}

pub(crate) fn write_header<W: Write>(writer: &mut W, header: Header) -> io::Result<()> {
    writer.write_all(MAGIC)?;
    write_u64(writer, VERSION)?;
    write_u64(writer, header.record_count)?;
    write_u64(writer, header.records_offset)?;
    write_u64(writer, header.strings_offset)?;
    write_u64(writer, header.strings_len)?;
    write_u64(writer, header.order_offset)?;
    write_u64(writer, header.order_len)?;
    Ok(())
}

pub(crate) fn write_record_entry<W: Write>(writer: &mut W, entry: RecordEntry) -> io::Result<()> {
    write_u64(writer, entry.full_id_offset)?;
    write_u64(writer, entry.full_id_len)?;
    write_u64(writer, entry.virtual_offset)?;
    write_u64(writer, entry.sequence_length)?;
    write_u64(writer, entry.line_bases)?;
    write_u64(writer, entry.line_width)?;
    write_u64(writer, 0)?;
    write_u64(writer, 0)?;
    Ok(())
}

pub(crate) fn read_u64<R: Read>(reader: &mut R) -> io::Result<u64> {
    let mut bytes = [0_u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

pub(crate) fn write_u64<W: Write>(writer: &mut W, value: u64) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

fn read_header(file: &mut File) -> io::Result<Header> {
    let mut bytes = [0_u8; HEADER_SIZE as usize];
    file.read_exact(&mut bytes)?;
    if &bytes[..8] != MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid or legacy .ffx index. Rebuild it with: fastars index --fasta FASTA.bgz",
        ));
    }
    let version = u64_from(&bytes, 8);
    if version != VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Unsupported .ffx version. Rebuild this index with the current fastars",
        ));
    }

    let header = Header {
        record_count: u64_from(&bytes, 16),
        records_offset: u64_from(&bytes, 24),
        strings_offset: u64_from(&bytes, 32),
        strings_len: u64_from(&bytes, 40),
        order_offset: u64_from(&bytes, 48),
        order_len: u64_from(&bytes, 56),
    };

    if header.records_offset != HEADER_SIZE || header.order_len != header.record_count * 8 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Corrupt .ffx header",
        ));
    }
    Ok(header)
}

fn read_record_entry<R: Read>(reader: &mut R) -> io::Result<RecordEntry> {
    let full_id_offset = read_u64(reader)?;
    let full_id_len = read_u64(reader)?;
    let virtual_offset = read_u64(reader)?;
    let sequence_length = read_u64(reader)?;
    let line_bases = read_u64(reader)?;
    let line_width = read_u64(reader)?;
    read_u64(reader)?;
    read_u64(reader)?;
    Ok(RecordEntry {
        full_id_offset,
        full_id_len,
        virtual_offset,
        sequence_length,
        line_bases,
        line_width,
    })
}

fn u64_from(bytes: &[u8], start: usize) -> u64 {
    u64::from_le_bytes(bytes[start..start + 8].try_into().unwrap())
}
