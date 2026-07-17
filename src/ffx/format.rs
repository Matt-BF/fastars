use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, Write};

const MAGIC: &[u8; 8] = b"FASTARS2";
const VERSION: u32 = 2;
pub(crate) const HEADER_SIZE: u64 = 80;
const DIRECTORY_FIXED_SIZE: u64 = 56;
pub(crate) const FLAG_COMPRESSED: u64 = 1;
pub(crate) const FLAG_FRONT_CODED: u64 = 1 << 1;
pub(crate) const FLAG_DELTA_OFFSETS: u64 = 1 << 2;
const KNOWN_FLAGS: u64 = FLAG_COMPRESSED | FLAG_FRONT_CODED | FLAG_DELTA_OFFSETS;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Header {
    pub record_count: u64,
    pub block_count: u64,
    pub blocks_offset: u64,
    pub blocks_len: u64,
    pub directory_offset: u64,
    pub directory_len: u64,
    pub target_records: u32,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct DirectoryEntry {
    pub block_offset: u64,
    pub stored_len: u64,
    pub raw_len: u64,
    pub first_record_id: u64,
    pub record_count: u64,
    pub flags: u64,
    pub last_id: String,
}

pub(crate) fn write_header<W: Write>(writer: &mut W, header: Header) -> io::Result<()> {
    writer.write_all(MAGIC)?;
    write_u32(writer, VERSION)?;
    write_u32(writer, HEADER_SIZE as u32)?;
    write_u64(writer, header.record_count)?;
    write_u64(writer, header.block_count)?;
    write_u64(writer, header.blocks_offset)?;
    write_u64(writer, header.blocks_len)?;
    write_u64(writer, header.directory_offset)?;
    write_u64(writer, header.directory_len)?;
    write_u32(writer, header.target_records)?;
    write_u32(writer, 0)?;
    write_u64(writer, 0)
}

pub(crate) fn read_header(file: &mut File) -> io::Result<Header> {
    file.seek(SeekFrom::Start(0))?;
    let mut bytes = [0_u8; HEADER_SIZE as usize];
    file.read_exact(&mut bytes)?;
    if &bytes[..8] != MAGIC {
        return invalid("Invalid or obsolete .ffx index; rebuild it with fastars index");
    }
    if u32_from(&bytes, 8) != VERSION || u32_from(&bytes, 12) != HEADER_SIZE as u32 {
        return invalid("Unsupported .ffx version; rebuild it with the current fastars");
    }
    if u32_from(&bytes, 68) != 0 || u64_from(&bytes, 72) != 0 {
        return invalid("Unsupported .ffx header flags");
    }

    Ok(Header {
        record_count: u64_from(&bytes, 16),
        block_count: u64_from(&bytes, 24),
        blocks_offset: u64_from(&bytes, 32),
        blocks_len: u64_from(&bytes, 40),
        directory_offset: u64_from(&bytes, 48),
        directory_len: u64_from(&bytes, 56),
        target_records: u32_from(&bytes, 64),
    })
}

pub(crate) fn write_directory<W: Write>(
    writer: &mut W,
    entries: &[DirectoryEntry],
) -> io::Result<u64> {
    let mut written = 0_u64;
    for entry in entries {
        let id_len = u64::try_from(entry.last_id.len())
            .map_err(|_| io::Error::other("FASTA ID is too long"))?;
        write_u64(writer, entry.block_offset)?;
        write_u64(writer, entry.stored_len)?;
        write_u64(writer, entry.raw_len)?;
        write_u64(writer, entry.first_record_id)?;
        write_u64(writer, entry.record_count)?;
        write_u64(writer, id_len)?;
        write_u64(writer, entry.flags)?;
        writer.write_all(entry.last_id.as_bytes())?;
        written = written
            .checked_add(DIRECTORY_FIXED_SIZE)
            .and_then(|value| value.checked_add(id_len))
            .ok_or_else(|| io::Error::other("Directory is too large"))?;
    }
    Ok(written)
}

pub(crate) fn read_directory(
    file: &mut File,
    header: Header,
    file_len: u64,
) -> io::Result<Vec<DirectoryEntry>> {
    validate_header(header, file_len)?;
    let capacity = usize::try_from(header.block_count)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Too many .ffx blocks"))?;
    let mut entries = Vec::with_capacity(capacity);
    let mut consumed = 0_u64;
    let mut expected_offset = header.blocks_offset;
    let mut expected_record_id = 0_u64;
    let mut previous_last_id: Option<String> = None;
    file.seek(SeekFrom::Start(header.directory_offset))?;

    for _ in 0..header.block_count {
        let block_offset = read_u64(file)?;
        let stored_len = read_u64(file)?;
        let raw_len = read_u64(file)?;
        let first_record_id = read_u64(file)?;
        let record_count = read_u64(file)?;
        let id_len = read_u64(file)?;
        let flags = read_u64(file)?;
        consumed = consumed
            .checked_add(DIRECTORY_FIXED_SIZE)
            .and_then(|value| value.checked_add(id_len))
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid directory size"))?;
        if consumed > header.directory_len {
            return invalid("Directory extends beyond its declared size");
        }
        let id_len = usize::try_from(id_len)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "FASTA ID is too long"))?;
        let mut id = vec![0_u8; id_len];
        file.read_exact(&mut id)?;
        let last_id = String::from_utf8(id)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Non-UTF-8 ID in .ffx"))?;

        if stored_len == 0 || raw_len == 0 || record_count == 0 {
            return invalid("Invalid empty .ffx block");
        }
        if flags & !KNOWN_FLAGS != 0 {
            return invalid("Unsupported .ffx block flags");
        }
        if flags & FLAG_COMPRESSED == 0 && stored_len != raw_len {
            return invalid("Raw .ffx block length does not match its payload");
        }
        if block_offset != expected_offset || first_record_id != expected_record_id {
            return invalid("Non-contiguous .ffx directory entry");
        }
        expected_offset = expected_offset
            .checked_add(stored_len)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid block offset"))?;
        expected_record_id = expected_record_id
            .checked_add(record_count)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid record count"))?;
        if let Some(previous) = &previous_last_id
            && previous > &last_id
        {
            return invalid("Block IDs are not lexically sorted");
        }
        previous_last_id = Some(last_id.clone());
        entries.push(DirectoryEntry {
            block_offset,
            stored_len,
            raw_len,
            first_record_id,
            record_count,
            flags,
            last_id,
        });
    }

    if consumed != header.directory_len
        || expected_offset != header.directory_offset
        || expected_record_id != header.record_count
    {
        return invalid(".ffx directory totals do not match the header");
    }
    Ok(entries)
}

fn validate_header(header: Header, file_len: u64) -> io::Result<()> {
    let blocks_end = header
        .blocks_offset
        .checked_add(header.blocks_len)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid block region"))?;
    let directory_end = header
        .directory_offset
        .checked_add(header.directory_len)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid directory region"))?;
    if header.blocks_offset != HEADER_SIZE
        || blocks_end != header.directory_offset
        || directory_end != file_len
        || (header.block_count == 0) != (header.record_count == 0)
        || header.target_records == 0
    {
        return invalid("Corrupt .ffx header");
    }
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

fn write_u32<W: Write>(writer: &mut W, value: u32) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

fn u32_from(bytes: &[u8], start: usize) -> u32 {
    u32::from_le_bytes(bytes[start..start + 4].try_into().unwrap())
}

fn u64_from(bytes: &[u8], start: usize) -> u64 {
    u64::from_le_bytes(bytes[start..start + 8].try_into().unwrap())
}

fn invalid<T>(message: &str) -> io::Result<T> {
    Err(io::Error::new(io::ErrorKind::InvalidData, message))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_is_exactly_the_declared_size() {
        let mut bytes = Vec::new();
        write_header(
            &mut bytes,
            Header {
                record_count: 0,
                block_count: 0,
                blocks_offset: HEADER_SIZE,
                blocks_len: 0,
                directory_offset: HEADER_SIZE,
                directory_len: 0,
                target_records: 4096,
            },
        )
        .unwrap();
        assert_eq!(bytes.len(), HEADER_SIZE as usize);
    }
}
