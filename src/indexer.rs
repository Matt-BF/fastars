use crate::bgzf::{BgzfLine, BgzfReader};
use crate::ffx::{IndexRecord, IndexWriter, read_u64};
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Clone, Copy)]
struct RecordMetadata {
    virtual_offset: u64,
    sequence_length: u64,
    line_bases: u64,
    line_width: u64,
}

pub fn build_index_from_fasta(
    fasta_path: &str,
    output_path: &str,
    temp_directory: Option<&str>,
) -> Result<(), Box<dyn Error>> {
    let mut reader = BgzfReader::new(fasta_path)?;
    let mut builder = IndexBuilder::new(output_path, temp_directory)?;
    let mut pending_header = None;

    loop {
        let header = match pending_header.take() {
            Some(line) => line,
            None => match reader.read_line()? {
                Some(line) => line,
                None => break,
            },
        };

        if !header.bytes.starts_with(b">") {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Expected FASTA header line starting with '>'",
            )
            .into());
        }

        let full_id = parse_header_id(&header.bytes)?;
        let (metadata, next_header) = read_record_metadata(&mut reader)?;
        builder.add_record(&full_id, metadata)?;
        pending_header = next_header;
    }

    builder.finish(output_path)?;
    Ok(())
}

pub fn build_index_from_fai_gzi(
    fai_path: &str,
    gzi_path: &str,
    output_path: &str,
    temp_directory: Option<&str>,
) -> Result<(), Box<dyn Error>> {
    let mut gzi = File::open(gzi_path)?;
    let gzi_count = read_gzi_count(&mut gzi)?;
    let reader = BufReader::new(File::open(fai_path)?);
    let mut builder = IndexBuilder::new(output_path, temp_directory)?;

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let (full_id, mut metadata) = parse_fai_record(&line)?;
        metadata.virtual_offset = gzi_virtual_offset(&mut gzi, gzi_count, metadata.virtual_offset)?;
        builder.add_record(&full_id, metadata)?;
    }

    builder.finish(output_path)?;
    Ok(())
}

fn read_record_metadata(reader: &mut BgzfReader) -> io::Result<(RecordMetadata, Option<BgzfLine>)> {
    let mut sequence_offset = None;
    let mut sequence_length = 0;
    let mut line_bases = 0;
    let mut line_width = 0;
    let mut previous_bases = None;
    let mut previous_width = None;

    while let Some(line) = reader.read_line()? {
        if line.bytes.starts_with(b">") {
            let metadata =
                finish_metadata(sequence_offset, sequence_length, line_bases, line_width)?;
            return Ok((metadata, Some(line)));
        }

        let (bases, width) = sequence_line_layout(&line.bytes)?;
        if bases == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Empty FASTA sequence lines are not supported",
            ));
        }

        if let (Some(previous_bases), Some(previous_width)) = (previous_bases, previous_width)
            && (previous_bases != line_bases || previous_width != line_width)
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Variable-width FASTA wrapping is not supported",
            ));
        }

        if sequence_offset.is_none() {
            sequence_offset = Some(line.start_virtual);
            line_bases = bases;
            line_width = width;
        }
        sequence_length += bases;
        previous_bases = Some(bases);
        previous_width = Some(width);
    }

    let metadata = finish_metadata(sequence_offset, sequence_length, line_bases, line_width)?;
    Ok((metadata, None))
}

fn finish_metadata(
    sequence_offset: Option<u64>,
    sequence_length: u64,
    line_bases: u64,
    line_width: u64,
) -> io::Result<RecordMetadata> {
    Ok(RecordMetadata {
        virtual_offset: sequence_offset.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "FASTA record has no sequence")
        })?,
        sequence_length,
        line_bases,
        line_width,
    })
}

fn parse_header_id(bytes: &[u8]) -> io::Result<String> {
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

    String::from_utf8(content[..id_end].to_vec())
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Non-UTF-8 FASTA ID"))
}

fn sequence_line_layout(bytes: &[u8]) -> io::Result<(u64, u64)> {
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
    sort_input: BufWriter<File>,
    record_count: u64,
    sort_input_path: PathBuf,
    sort_output_path: PathBuf,
    temp_directory: Option<PathBuf>,
}

impl IndexBuilder {
    fn new(output_path: &str, temp_directory: Option<&str>) -> io::Result<Self> {
        let output = Path::new(output_path);
        let temp_root = temp_directory
            .map(PathBuf::from)
            .or_else(|| output.parent().map(PathBuf::from))
            .filter(|path| !path.as_os_str().is_empty())
            .unwrap_or_else(|| PathBuf::from("."));
        fs::create_dir_all(&temp_root)?;

        let pid = std::process::id();
        let sort_input_path = temp_root.join(format!("fastars-{pid}.sort-input.tmp"));
        let sort_output_path = temp_root.join(format!("fastars-{pid}.sort-output.tmp"));

        Ok(Self {
            sort_input: BufWriter::new(File::create(&sort_input_path)?),
            record_count: 0,
            sort_input_path,
            sort_output_path,
            temp_directory: temp_directory.map(PathBuf::from),
        })
    }

    fn add_record(&mut self, full_id: &str, metadata: RecordMetadata) -> io::Result<()> {
        if full_id.contains(['\t', '\n', '\r']) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTA IDs containing tabs or newlines are not supported",
            ));
        }

        writeln!(
            self.sort_input,
            "{full_id}\t{}\t{}\t{}\t{}",
            metadata.virtual_offset,
            metadata.sequence_length,
            metadata.line_bases,
            metadata.line_width
        )?;

        self.record_count = self
            .record_count
            .checked_add(1)
            .ok_or_else(|| io::Error::other("Too many FASTA records"))?;
        Ok(())
    }

    fn finish(mut self, output_path: &str) -> Result<(), Box<dyn Error>> {
        self.sort_input.flush()?;

        let mut sort = Command::new("sort");
        sort.env("LC_ALL", "C")
            .arg("-t")
            .arg("\t")
            .arg("-k1,1")
            .arg("-s");
        if let Some(directory) = &self.temp_directory {
            sort.arg("-T").arg(directory);
        }
        let status = sort
            .arg(&self.sort_input_path)
            .arg("-o")
            .arg(&self.sort_output_path)
            .status()?;
        if !status.success() {
            return Err(io::Error::other("sort failed while building .ffx").into());
        }

        let sorted = BufReader::new(File::open(&self.sort_output_path)?);
        let mut writer = IndexWriter::new(output_path)?;
        let mut ordered_count = 0;
        for line in sorted.lines() {
            writer.add(parse_sorted_record(&line?)?)?;
            ordered_count += 1;
        }

        if ordered_count != self.record_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Sorted index record count does not match scanned FASTA records",
            )
            .into());
        }
        let stats = writer.finish()?;
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
        Ok(())
    }
}

impl Drop for IndexBuilder {
    fn drop(&mut self) {
        for path in [&self.sort_input_path, &self.sort_output_path] {
            let _ = fs::remove_file(path);
        }
    }
}

fn parse_sorted_record(line: &str) -> io::Result<IndexRecord> {
    let mut fields = line.split('\t');
    let full_id = fields.next().unwrap_or_default();
    let virtual_offset = parse_sorted_u64(fields.next(), "virtual offset")?;
    let sequence_length = parse_sorted_u64(fields.next(), "sequence length")?;
    let line_bases = parse_sorted_u64(fields.next(), "line bases")?;
    let line_width = parse_sorted_u64(fields.next(), "line width")?;
    if full_id.is_empty() || fields.next().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Malformed sorted index row",
        ));
    }
    Ok(IndexRecord {
        full_id: full_id.to_string(),
        virtual_offset,
        sequence_length,
        line_bases,
        line_width,
    })
}

fn parse_sorted_u64(value: Option<&str>, field: &str) -> io::Result<u64> {
    value
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, format!("Invalid {field}")))
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
    fn sorted_row_parser_rejects_extra_fields() {
        assert!(parse_sorted_record("id\t1\t2\t3\t4\textra").is_err());
        assert!(parse_sorted_record("id\t1\t2\t3").is_err());
    }

    #[test]
    fn dropping_builder_removes_temporary_files() {
        let directory = test_path("temporary-directory");
        fs::create_dir(&directory).unwrap();
        let output = directory.join("output.ffx");
        let builder = IndexBuilder::new(output.to_str().unwrap(), directory.to_str()).unwrap();
        let input = builder.sort_input_path.clone();
        let sorted = builder.sort_output_path.clone();
        drop(builder);
        assert!(!input.exists());
        assert!(!sorted.exists());
        fs::remove_dir(directory).unwrap();
    }
}
