use super::FfxRecord;
use super::format::{DirectoryEntry, FLAG_COMPRESSED, FLAG_DELTA_OFFSETS, FLAG_FRONT_CODED};
use std::io::{self, Cursor, Read, Write};

const ZSTD_LEVEL: i32 = 3;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct IndexRecord {
    pub full_id: String,
    pub virtual_offset: u64,
    pub sequence_length: u64,
    pub line_bases: u64,
    pub line_width: u64,
}

pub(crate) struct EncodedBlock {
    pub stored: Vec<u8>,
    pub raw_len: u64,
    pub flags: u64,
}

pub(crate) fn encode_block(records: &[IndexRecord]) -> io::Result<EncodedBlock> {
    if records.is_empty() {
        return invalid("Cannot encode an empty .ffx block");
    }

    let mut candidates = Vec::with_capacity(4);
    for front_coded in [false, true] {
        for delta_offsets in [false, true] {
            if let Some(raw) = encode_records(records, front_coded, delta_offsets)? {
                let mut flags = 0;
                if front_coded {
                    flags |= FLAG_FRONT_CODED;
                }
                if delta_offsets {
                    flags |= FLAG_DELTA_OFFSETS;
                }
                candidates.push((raw, flags));
            }
        }
    }
    let (raw, mut flags) = candidates
        .into_iter()
        .min_by_key(|(bytes, _)| bytes.len())
        .ok_or_else(|| io::Error::other("Could not encode .ffx block"))?;
    let raw_len = u64::try_from(raw.len()).map_err(|_| io::Error::other("Block is too large"))?;
    let mut encoder = zstd::stream::Encoder::new(Vec::new(), ZSTD_LEVEL)?;
    encoder.include_checksum(true)?;
    encoder.write_all(&raw)?;
    let compressed = encoder.finish()?;
    let stored = if compressed.len() < raw.len() {
        flags |= FLAG_COMPRESSED;
        compressed
    } else {
        raw
    };
    Ok(EncodedBlock {
        stored,
        raw_len,
        flags,
    })
}

pub(crate) fn decode_block(stored: &[u8], entry: &DirectoryEntry) -> io::Result<Vec<FfxRecord>> {
    let raw = if entry.flags & FLAG_COMPRESSED != 0 {
        decompress(stored, entry.raw_len)?
    } else {
        if u64::try_from(stored.len()).ok() != Some(entry.raw_len) {
            return invalid("Raw .ffx block has an invalid length");
        }
        stored.to_vec()
    };
    decode_records(&raw, entry)
}

fn encode_records(
    records: &[IndexRecord],
    front_coded: bool,
    delta_offsets: bool,
) -> io::Result<Option<Vec<u8>>> {
    let mut output = Vec::new();
    let mut previous_id: &[u8] = &[];
    let mut previous_offset = 0_u64;

    for (position, record) in records.iter().enumerate() {
        let id = record.full_id.as_bytes();
        if front_coded && position > 0 {
            let prefix_len = common_prefix(previous_id, id);
            write_varint(&mut output, prefix_len as u64);
            write_varint(&mut output, (id.len() - prefix_len) as u64);
            output.extend_from_slice(&id[prefix_len..]);
        } else {
            write_varint(&mut output, id.len() as u64);
            output.extend_from_slice(id);
        }

        if delta_offsets && position > 0 {
            let difference = record.virtual_offset as i128 - previous_offset as i128;
            let Ok(difference) = i64::try_from(difference) else {
                return Ok(None);
            };
            write_varint(&mut output, zigzag_encode(difference));
        } else {
            write_varint(&mut output, record.virtual_offset);
        }
        write_varint(&mut output, record.sequence_length);
        write_varint(&mut output, record.line_bases);
        write_varint(&mut output, record.line_width);
        previous_id = id;
        previous_offset = record.virtual_offset;
    }
    Ok(Some(output))
}

fn decode_records(raw: &[u8], entry: &DirectoryEntry) -> io::Result<Vec<FfxRecord>> {
    let capacity = usize::try_from(entry.record_count)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Too many records in block"))?;
    let mut records = Vec::with_capacity(capacity);
    let mut cursor = 0_usize;
    let mut previous_id = Vec::new();
    let mut previous_offset = 0_u64;

    for position in 0..entry.record_count {
        let id = if entry.flags & FLAG_FRONT_CODED != 0 && position > 0 {
            let prefix_len = read_length(raw, &mut cursor)?;
            let suffix_len = read_length(raw, &mut cursor)?;
            if prefix_len > previous_id.len() {
                return invalid("Invalid front-coded ID prefix");
            }
            let suffix = take(raw, &mut cursor, suffix_len)?;
            let mut id = Vec::with_capacity(prefix_len + suffix_len);
            id.extend_from_slice(&previous_id[..prefix_len]);
            id.extend_from_slice(suffix);
            id
        } else {
            let id_len = read_length(raw, &mut cursor)?;
            take(raw, &mut cursor, id_len)?.to_vec()
        };
        let full_id = String::from_utf8(id.clone())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Non-UTF-8 ID in .ffx"))?;
        if let Some(previous) = records.last().map(|record: &FfxRecord| &record.full_id)
            && previous > &full_id
        {
            return invalid("IDs within .ffx block are not sorted");
        }

        let encoded_offset = read_varint(raw, &mut cursor)?;
        let virtual_offset = if entry.flags & FLAG_DELTA_OFFSETS != 0 && position > 0 {
            let difference = zigzag_decode(encoded_offset) as i128;
            let offset = previous_offset as i128 + difference;
            if !(0..=u64::MAX as i128).contains(&offset) {
                return invalid("Virtual-offset delta overflows u64");
            }
            offset as u64
        } else {
            encoded_offset
        };
        let sequence_length = read_varint(raw, &mut cursor)?;
        let line_bases = read_varint(raw, &mut cursor)?;
        let line_width = read_varint(raw, &mut cursor)?;
        if line_bases == 0 || line_width < line_bases {
            return invalid("Invalid FASTA line layout in .ffx block");
        }
        let record_id = entry
            .first_record_id
            .checked_add(position)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Record ID overflow"))?;
        records.push(FfxRecord {
            record_id,
            full_id,
            virtual_offset,
            sequence_length,
            line_bases,
            line_width,
        });
        previous_id = id;
        previous_offset = virtual_offset;
    }

    if cursor != raw.len() || records.last().map(|record| &record.full_id) != Some(&entry.last_id) {
        return invalid(".ffx block contents do not match its directory entry");
    }
    Ok(records)
}

fn decompress(stored: &[u8], raw_len: u64) -> io::Result<Vec<u8>> {
    let limit = raw_len
        .checked_add(1)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid raw block length"))?;
    let decoder = zstd::stream::read::Decoder::new(Cursor::new(stored))?;
    let mut raw = Vec::new();
    decoder.take(limit).read_to_end(&mut raw)?;
    if u64::try_from(raw.len()).ok() != Some(raw_len) {
        return invalid("Decompressed .ffx block length does not match its directory entry");
    }
    Ok(raw)
}

fn common_prefix(left: &[u8], right: &[u8]) -> usize {
    left.iter()
        .zip(right)
        .take_while(|(left, right)| left == right)
        .count()
}

fn write_varint(output: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        output.push((value as u8 & 0x7f) | 0x80);
        value >>= 7;
    }
    output.push(value as u8);
}

fn read_varint(input: &[u8], cursor: &mut usize) -> io::Result<u64> {
    let mut value = 0_u64;
    for shift in (0..=63).step_by(7) {
        let byte = *input
            .get(*cursor)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Truncated varint"))?;
        *cursor += 1;
        if shift == 63 && byte > 1 {
            return invalid("Varint overflows u64");
        }
        value |= u64::from(byte & 0x7f) << shift;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
    }
    invalid("Varint is too long")
}

fn read_length(input: &[u8], cursor: &mut usize) -> io::Result<usize> {
    usize::try_from(read_varint(input, cursor)?)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Length exceeds address space"))
}

fn take<'a>(input: &'a [u8], cursor: &mut usize, len: usize) -> io::Result<&'a [u8]> {
    let end = cursor
        .checked_add(len)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid block length"))?;
    let bytes = input
        .get(*cursor..end)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Truncated block"))?;
    *cursor = end;
    Ok(bytes)
}

fn zigzag_encode(value: i64) -> u64 {
    ((value << 1) ^ (value >> 63)) as u64
}

fn zigzag_decode(value: u64) -> i64 {
    ((value >> 1) as i64) ^ -((value & 1) as i64)
}

fn invalid<T>(message: &str) -> io::Result<T> {
    Err(io::Error::new(io::ErrorKind::InvalidData, message))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(id: &str, offset: u64) -> IndexRecord {
        IndexRecord {
            full_id: id.to_string(),
            virtual_offset: offset,
            sequence_length: 100_000,
            line_bases: 100_000,
            line_width: 100_001,
        }
    }

    fn entry(encoded: &EncodedBlock, records: &[IndexRecord]) -> DirectoryEntry {
        DirectoryEntry {
            block_offset: 80,
            stored_len: encoded.stored.len() as u64,
            raw_len: encoded.raw_len,
            first_record_id: 7,
            record_count: records.len() as u64,
            flags: encoded.flags,
            last_id: records.last().unwrap().full_id.clone(),
        }
    }

    #[test]
    fn varints_round_trip_boundaries() {
        for value in [0, 1, 127, 128, 16_383, 16_384, u32::MAX as u64, u64::MAX] {
            let mut bytes = Vec::new();
            write_varint(&mut bytes, value);
            let mut cursor = 0;
            assert_eq!(read_varint(&bytes, &mut cursor).unwrap(), value);
            assert_eq!(cursor, bytes.len());
        }
    }

    #[test]
    fn zigzag_round_trips_boundaries() {
        for value in [i64::MIN, -1, 0, 1, i64::MAX] {
            assert_eq!(zigzag_decode(zigzag_encode(value)), value);
        }
    }

    #[test]
    fn block_round_trip_preserves_records_and_large_layouts() {
        let records = vec![
            record("α-long-prefix-1", 1 << 40),
            record("α-long-prefix-2", (1 << 40) + 10),
        ];
        let encoded = encode_block(&records).unwrap();
        let decoded = decode_block(&encoded.stored, &entry(&encoded, &records)).unwrap();
        assert_eq!(decoded.len(), 2);
        assert_eq!(decoded[0].record_id, 7);
        assert_eq!(decoded[1].full_id, records[1].full_id);
        assert_eq!(decoded[1].line_bases, 100_000);
    }

    #[test]
    fn repetitive_ids_select_front_coding() {
        let records = (0..100)
            .map(|value| record(&format!("shared-prefix-{value:04}"), value))
            .collect::<Vec<_>>();
        let encoded = encode_block(&records).unwrap();
        assert_ne!(encoded.flags & FLAG_FRONT_CODED, 0);
    }

    #[test]
    fn unrelated_ids_and_extreme_offsets_use_bounded_fallbacks() {
        let records = vec![
            record("a-random-value", 0),
            record("z-other-value", u64::MAX),
        ];
        let encoded = encode_block(&records).unwrap();
        assert_eq!(encoded.flags & FLAG_FRONT_CODED, 0);
        assert_eq!(encoded.flags & FLAG_DELTA_OFFSETS, 0);
        let decoded = decode_block(&encoded.stored, &entry(&encoded, &records)).unwrap();
        assert_eq!(decoded[1].virtual_offset, u64::MAX);
    }

    #[test]
    fn small_incompressible_block_uses_raw_fallback() {
        let records = vec![record("x", u64::MAX)];
        let encoded = encode_block(&records).unwrap();
        assert_eq!(encoded.flags & FLAG_COMPRESSED, 0);
    }

    #[test]
    fn compressed_block_checksum_detects_corruption() {
        let records = (0..100)
            .map(|value| record(&format!("shared-prefix-{value:04}"), value))
            .collect::<Vec<_>>();
        let mut encoded = encode_block(&records).unwrap();
        assert_ne!(encoded.flags & FLAG_COMPRESSED, 0);
        *encoded.stored.last_mut().unwrap() ^= 1;
        assert!(decode_block(&encoded.stored, &entry(&encoded, &records)).is_err());
    }
}
