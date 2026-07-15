use crate::bgzf::{BgzfLine, BgzfReader};
use crate::ffx::{
    HEADER_SIZE, Header, RECORD_SIZE, RecordEntry, read_u64, write_header, write_record_entry,
    write_u64,
};
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

        if let (Some(previous_bases), Some(previous_width)) = (previous_bases, previous_width) {
            if previous_bases != line_bases || previous_width != line_width {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Variable-width FASTA wrapping is not supported",
                ));
            }
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
    records: BufWriter<File>,
    strings: BufWriter<File>,
    sort_input: BufWriter<File>,
    record_count: u64,
    strings_len: u64,
    records_path: PathBuf,
    strings_path: PathBuf,
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
        let records_path = temp_root.join(format!("fastars-{pid}.records.tmp"));
        let strings_path = temp_root.join(format!("fastars-{pid}.strings.tmp"));
        let sort_input_path = temp_root.join(format!("fastars-{pid}.sort-input.tmp"));
        let sort_output_path = temp_root.join(format!("fastars-{pid}.sort-output.tmp"));

        Ok(Self {
            records: BufWriter::new(File::create(&records_path)?),
            strings: BufWriter::new(File::create(&strings_path)?),
            sort_input: BufWriter::new(File::create(&sort_input_path)?),
            record_count: 0,
            strings_len: 0,
            records_path,
            strings_path,
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

        let record_id = self.record_count;
        let id_bytes = full_id.as_bytes();
        let entry = RecordEntry {
            full_id_offset: self.strings_len,
            full_id_len: id_bytes.len() as u64,
            virtual_offset: metadata.virtual_offset,
            sequence_length: metadata.sequence_length,
            line_bases: metadata.line_bases,
            line_width: metadata.line_width,
        };

        write_record_entry(&mut self.records, entry)?;
        self.strings.write_all(id_bytes)?;
        writeln!(self.sort_input, "{full_id}\t{record_id}")?;

        self.record_count += 1;
        self.strings_len += id_bytes.len() as u64;
        Ok(())
    }

    fn finish(self, output_path: &str) -> Result<(), Box<dyn Error>> {
        let records_path = self.records_path.clone();
        let strings_path = self.strings_path.clone();
        let sort_input_path = self.sort_input_path.clone();
        let sort_output_path = self.sort_output_path.clone();
        let result = self.finish_inner(output_path);

        for path in [
            records_path.as_path(),
            strings_path.as_path(),
            sort_input_path.as_path(),
            sort_output_path.as_path(),
        ] {
            let _ = fs::remove_file(path);
        }

        result
    }

    fn finish_inner(self, output_path: &str) -> Result<(), Box<dyn Error>> {
        let IndexBuilder {
            mut records,
            mut strings,
            mut sort_input,
            record_count,
            strings_len,
            records_path,
            strings_path,
            sort_input_path,
            sort_output_path,
            temp_directory,
        } = self;

        records.flush()?;
        strings.flush()?;
        sort_input.flush()?;
        drop(records);
        drop(strings);
        drop(sort_input);

        let mut sort = Command::new("sort");
        sort.env("LC_ALL", "C")
            .arg("-t")
            .arg("\t")
            .arg("-k1,1")
            .arg("-s");
        if let Some(directory) = temp_directory {
            sort.arg("-T").arg(directory);
        }
        let status = sort
            .arg(&sort_input_path)
            .arg("-o")
            .arg(&sort_output_path)
            .status()?;
        if !status.success() {
            return Err(io::Error::other("sort failed while building .ffx").into());
        }

        let records_offset = HEADER_SIZE;
        let strings_offset = records_offset + record_count * RECORD_SIZE;
        let order_offset = strings_offset + strings_len;
        let order_len = record_count * 8;
        let header = Header {
            record_count,
            records_offset,
            strings_offset,
            strings_len,
            order_offset,
            order_len,
        };

        let mut output = BufWriter::new(File::create(output_path)?);
        write_header(&mut output, header)?;

        let mut records_file = File::open(&records_path)?;
        io::copy(&mut records_file, &mut output)?;

        let mut strings_file = File::open(&strings_path)?;
        io::copy(&mut strings_file, &mut output)?;

        let sorted = BufReader::new(File::open(&sort_output_path)?);
        let mut ordered_count = 0;
        for line in sorted.lines() {
            let line = line?;
            let record_id = line
                .rsplit_once('\t')
                .and_then(|(_, value)| value.parse::<u64>().ok())
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "Malformed sorted index row")
                })?;
            write_u64(&mut output, record_id)?;
            ordered_count += 1;
        }

        if ordered_count != record_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Sorted index record count does not match scanned FASTA records",
            )
            .into());
        }

        output.flush()?;
        Ok(())
    }
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
