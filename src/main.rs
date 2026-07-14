use regex::Regex;
use std::env;
use std::error::Error;
use std::ffi::c_void;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::mem::{size_of, zeroed};
use std::os::raw::{c_char, c_int, c_uint, c_ulong};
use std::path::Path;
use std::process::Command;

const Z_STREAM_END: c_int = 1;
const Z_FINISH: c_int = 4;

#[repr(C)]
struct ZStream {
    next_in: *mut u8,
    avail_in: c_uint,
    total_in: c_ulong,
    next_out: *mut u8,
    avail_out: c_uint,
    total_out: c_ulong,
    msg: *mut c_char,
    state: *mut c_void,
    zalloc: *mut c_void,
    zfree: *mut c_void,
    opaque: *mut c_void,
    data_type: c_int,
    adler: c_ulong,
    reserved: c_ulong,
}

#[link(name = "z")]
unsafe extern "C" {
    fn zlibVersion() -> *const c_char;
    fn inflateInit2_(
        stream: *mut ZStream,
        window_bits: c_int,
        version: *const c_char,
        size: c_int,
    ) -> c_int;
    fn inflate(stream: *mut ZStream, flush: c_int) -> c_int;
    fn inflateEnd(stream: *mut ZStream) -> c_int;
}

#[derive(Clone)]
struct Record {
    full_id: String,
    sequence_length: u64,
    sequence_offset: u64,
    line_bases: u64,
    line_width: u64,
    virtual_offset: u64,
}

enum Query {
    Short(String),
    Full(String),
}

#[derive(Clone, Copy, PartialEq)]
enum IndexKeyKind {
    Full,
    Short,
}

impl IndexKeyKind {
    fn marker(self) -> &'static str {
        match self {
            Self::Full => "F",
            Self::Short => "S",
        }
    }
}

struct Arguments {
    fasta: String,
    fai: String,
    gzi: String,
    ffx: String,
    ids_file: Option<String>,
    sort_by_offset: bool,
    verbose_missing: bool,
    short_id: bool,
    ids: Vec<String>,
}

#[allow(dead_code)]
fn usage_legacy() -> &'static str {
    "Usage: fai-fetch --fasta FASTA.bgz [--fai FASTA.bgz.fai] [--gzi FASTA.bgz.gzi] \\\n+        [--ids-file IDs.txt] [--sort-by-offset] ID [ID ...]\n\n\
The .fai must be lexically sorted by full sequence ID. FASTA records are written to stdout."
}

fn usage() -> &'static str {
    r#"fastars — fast lookup and retrieval for BGZF-compressed FASTA files

USAGE:
    fastars --fasta <FASTA.bgz> [OPTIONS] <ID>...
    fastars index --fai <FASTA.bgz.fai> [OPTIONS]

COMMANDS:
    index    Build a sorted .ffx lookup index from a FASTA .fai file.

FETCH OPTIONS:
    --fasta <FILE>       BGZF-compressed FASTA to search. Required.
    --fai <FILE>         FASTA index. Default: <FASTA>.fai
    --gzi <FILE>         BGZF index. Default: <FASTA>.gzi
    --ffx <FILE>         Lookup index. Default: <FASTA>.fai.ffx
    --ids-file <FILE>    Read one requested ID per line, in addition to IDs given directly.
    --short-id           Treat requested IDs as short IDs indexed with --short-id-regex.
    --full-id            Treat requested IDs as exact full FASTA IDs. Default.
    --sort-by-offset     Fetch records in FASTA order to reduce random disk access.
    --verbose-missing    Print every requested ID that was not found.
    -h, --help           Print this help message.

INDEX OPTIONS:
    --fai <FILE>                 Plain FASTA .fai file to index. Required.
    --output <FILE>              Output lookup index. Default: <FAI>.ffx
    --short-id-regex <REGEX>     Add short-ID keys using capture group 1.
    --temp-directory <DIR>       Directory for sort temporary files. Default: current directory
    -h, --help                   Print index-specific help.

Full FASTA IDs are always indexed. A short-ID lookup requires an index built
with --short-id-regex. FASTA records are written to standard output.
"#
}

fn index_usage() -> &'static str {
    r#"fastars index — build a sorted lookup index

USAGE:
    fastars index --fai <FASTA.bgz.fai> [OPTIONS]

OPTIONS:
    --fai <FILE>                 Plain FASTA .fai file to index. Required.
    --output <FILE>              Output lookup index. Default: <FAI>.ffx
    --short-id-regex <REGEX>     Add short-ID keys using capture group 1.
    --temp-directory <DIR>       Directory for sort temporary files. Default: current directory
    -h, --help                   Print this help message.

The index always contains exact full FASTA IDs. If --short-id-regex is set,
capture group 1 is added as an additional short-ID key.
"#
}

fn default_ffx_path(fai: &str) -> String {
    let adjacent = format!("{fai}.ffx");
    if Path::new(&adjacent).exists() {
        return adjacent;
    }
    if let Some(name) = Path::new(&adjacent)
        .file_name()
        .and_then(|name| name.to_str())
    {
        if Path::new(name).exists() {
            return name.to_string();
        }
        return name.to_string();
    }
    adjacent
}

fn parse_arguments() -> Result<Arguments, String> {
    let mut fasta = None;
    let mut fai = None;
    let mut gzi = None;
    let mut ffx = None;
    let mut ids_file = None;
    let mut sort_by_offset = false;
    let mut verbose_missing = false;
    let mut short_id = false;
    let mut ids = Vec::new();
    let mut arguments = env::args().skip(1);

    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--fasta" => fasta = arguments.next(),
            "--fai" => fai = arguments.next(),
            "--gzi" => gzi = arguments.next(),
            "--ffx" => ffx = arguments.next(),
            "--ids-file" => ids_file = arguments.next(),
            "--sort-by-offset" => sort_by_offset = true,
            "--verbose-missing" => verbose_missing = true,
            "--short-id" => short_id = true,
            "--full-id" => short_id = false,
            "-h" | "--help" => return Err(usage().to_string()),
            _ if argument.starts_with('-') => {
                return Err(format!("Unknown option: {argument}\n\n{}", usage()));
            }
            _ => ids.push(argument),
        }
    }

    let fasta = fasta.ok_or_else(|| format!("--fasta is required\n\n{}", usage()))?;
    let fai = fai.unwrap_or_else(|| format!("{fasta}.fai"));
    Ok(Arguments {
        ffx: ffx.unwrap_or_else(|| default_ffx_path(&fai)),
        fai,
        gzi: gzi.unwrap_or_else(|| format!("{fasta}.gzi")),
        fasta,
        ids_file,
        sort_by_offset,
        verbose_missing,
        short_id,
        ids,
    })
}

fn read_queries(arguments: &Arguments) -> io::Result<Vec<Query>> {
    let mut ids = arguments.ids.clone();
    if let Some(path) = &arguments.ids_file {
        let reader = BufReader::new(File::open(path)?);
        for line in reader.lines() {
            let value = line?;
            if !value.trim().is_empty() {
                ids.push(value);
            }
        }
    }

    Ok(ids
        .into_iter()
        .map(|value| {
            let id = value
                .trim()
                .trim_start_matches('>')
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            if arguments.short_id {
                Query::Short(id)
            } else {
                Query::Full(id)
            }
        })
        .collect())
}

fn read_line_at_or_after(file: &mut File, offset: u64) -> io::Result<Option<(u64, u64, String)>> {
    let length = file.metadata()?.len();
    if offset >= length {
        return Ok(None);
    }

    let mut start = offset;
    if offset > 0 {
        file.seek(SeekFrom::Start(offset - 1))?;
        let mut previous = [0_u8; 1];
        file.read_exact(&mut previous)?;
        if previous[0] != b'\n' {
            file.seek(SeekFrom::Start(offset))?;
            let mut byte = [0_u8; 1];
            while file.read(&mut byte)? == 1 {
                if byte[0] == b'\n' {
                    break;
                }
            }
            start = file.stream_position()?;
            if start >= length {
                return Ok(None);
            }
        }
    }

    file.seek(SeekFrom::Start(start))?;
    let mut bytes = Vec::new();
    let mut byte = [0_u8; 1];
    while file.read(&mut byte)? == 1 {
        if byte[0] == b'\n' {
            break;
        }
        bytes.push(byte[0]);
    }
    let next = file.stream_position()?;
    let line = String::from_utf8(bytes)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Non-UTF-8 .fai entry"))?;
    Ok(Some((start, next, line)))
}

fn validate_lookup_index(index: &mut File) -> io::Result<()> {
    let line = read_line_at_or_after(index, 0)?
        .map(|(_, _, line)| line)
        .unwrap_or_default();
    let mut fields = line.split('\t');
    let is_valid = fields.next().is_some()
        && matches!(fields.next(), Some("F" | "S"))
        && fields
            .next()
            .and_then(|value| value.parse::<u64>().ok())
            .is_some();
    if !is_valid {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Legacy or invalid .ffx index. Rebuild it with: fastars index --fai FASTA.bgz.fai [--short-id-regex REGEX]",
        ));
    }
    Ok(())
}

fn fai_key(line: &str) -> &str {
    line.split('\t').next().unwrap_or("")
}

fn lower_bound(file: &mut File, target: &str) -> io::Result<Option<(u64, u64, String)>> {
    let mut low = 0;
    let mut high = file.metadata()?.len();
    while low < high {
        let middle = low + (high - low) / 2;
        let Some((start, next, line)) = read_line_at_or_after(file, middle)? else {
            high = middle;
            continue;
        };
        if fai_key(&line) < target {
            low = next;
        } else if start < high {
            high = start;
        } else {
            high = middle;
        }
    }
    read_line_at_or_after(file, low)
}

fn parse_record(line: &str) -> io::Result<Record> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 5 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Malformed .fai line",
        ));
    }
    Ok(Record {
        full_id: fields[0].to_string(),
        sequence_length: fields[1]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid sequence length"))?,
        sequence_offset: fields[2]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid sequence offset"))?,
        line_bases: fields[3]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid line base count"))?,
        line_width: fields[4]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid line width"))?,
        virtual_offset: 0,
    })
}

fn find_index_offsets(index: &mut File, key: &str, kind: IndexKeyKind) -> io::Result<Vec<u64>> {
    let Some((_, mut next, line)) = lower_bound(index, key)? else {
        return Ok(Vec::new());
    };
    let mut offsets = Vec::new();
    let mut current = Some(line);
    while let Some(line) = current {
        let mut fields = line.split('\t');
        if fields.next() != Some(key) {
            break;
        }
        let entry_kind = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Malformed .ffx line"))?
            .to_string();
        let offset = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Malformed .ffx line"))?
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid .ffx offset"))?;
        if entry_kind == kind.marker() {
            offsets.push(offset);
        }
        current = read_line_at_or_after(index, next)?.map(|(_, position, line)| {
            next = position;
            line
        });
    }
    Ok(offsets)
}

fn extract_short_id(regex: &Regex, full_id: &str) -> io::Result<String> {
    regex
        .captures(full_id)
        .and_then(|captures| captures.get(1))
        .map(|capture| capture.as_str().to_string())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Short-ID regex did not match capture group 1",
            )
        })
}

fn find_records(
    fai: &mut File,
    index: &mut File,
    key: &str,
    kind: IndexKeyKind,
) -> io::Result<Vec<Record>> {
    let mut records = Vec::new();
    for offset in find_index_offsets(index, key, kind)? {
        let (_, _, line) = read_line_at_or_after(fai, offset)?
            .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "Invalid .ffx offset"))?;
        records.push(parse_record(&line)?);
    }
    Ok(records)
}

fn find_full(fai: &mut File, index: &mut File, full_id: &str) -> io::Result<Vec<Record>> {
    Ok(find_records(fai, index, full_id, IndexKeyKind::Full)?
        .into_iter()
        .filter(|record| record.full_id == full_id)
        .collect())
}

fn build_lookup_index(
    fai_path: &str,
    output_path: &str,
    temp_directory: Option<&str>,
    short_id_regex: Option<&Regex>,
) -> Result<(), Box<dyn Error>> {
    let temp_root = temp_directory.unwrap_or(".");
    let temp_path = format!("{temp_root}/fastars-{}.unsorted", std::process::id());
    let mut input = BufReader::new(File::open(fai_path)?);
    let mut output = BufWriter::new(File::create(&temp_path)?);
    let mut offset = 0_u64;
    let mut line = Vec::new();
    loop {
        line.clear();
        let bytes_read = input.read_until(b'\n', &mut line)?;
        if bytes_read == 0 {
            break;
        }
        let text = std::str::from_utf8(&line)?.trim_end_matches(['\r', '\n']);
        let full_id = fai_key(text);
        writeln!(
            output,
            "{full_id}\t{}\t{offset}",
            IndexKeyKind::Full.marker()
        )?;
        if let Some(regex) = short_id_regex {
            let short_id = extract_short_id(regex, full_id)?;
            if short_id != full_id {
                writeln!(
                    output,
                    "{short_id}\t{}\t{offset}",
                    IndexKeyKind::Short.marker()
                )?;
            }
        }
        offset += bytes_read as u64;
    }
    output.flush()?;

    let mut sort = Command::new("sort");
    sort.env("LC_ALL", "C").arg("-t").arg("\t").arg("-k1,1");
    if let Some(directory) = temp_directory {
        sort.arg("-T").arg(directory);
    }
    let status = sort.arg(&temp_path).arg("-o").arg(output_path).status()?;
    fs::remove_file(&temp_path)?;
    if !status.success() {
        return Err(io::Error::other("sort failed while building .ffx index").into());
    }
    Ok(())
}

fn run_index_command() -> Result<(), Box<dyn Error>> {
    let mut fai = None;
    let mut output = None;
    let mut temp_directory = None;
    let mut short_id_regex = None;
    let mut arguments = env::args().skip(2);
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--fai" => fai = arguments.next(),
            "--output" => output = arguments.next(),
            "--temp-directory" => temp_directory = arguments.next(),
            "--short-id-regex" => {
                short_id_regex = Some(arguments.next().ok_or("--short-id-regex needs a value")?);
            }
            "-h" | "--help" => {
                println!("{}", index_usage());
                return Ok(());
            }
            _ => return Err(io::Error::other(format!("Unknown index option: {argument}")).into()),
        }
    }
    let fai = fai.ok_or_else(|| io::Error::other("index requires --fai"))?;
    let output = output.unwrap_or_else(|| format!("{fai}.ffx"));
    let regex = short_id_regex.as_deref().map(Regex::new).transpose()?;
    build_lookup_index(&fai, &output, temp_directory.as_deref(), regex.as_ref())?;
    eprintln!("Wrote {output}");
    Ok(())
}

fn read_gzi_pair(file: &mut File, index: u64) -> io::Result<(u64, u64)> {
    file.seek(SeekFrom::Start(8 + index * 16))?;
    let mut bytes = [0_u8; 16];
    file.read_exact(&mut bytes)?;
    Ok((
        u64::from_le_bytes(bytes[..8].try_into().unwrap()),
        u64::from_le_bytes(bytes[8..].try_into().unwrap()),
    ))
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

fn read_gzi_count(file: &mut File) -> io::Result<u64> {
    file.seek(SeekFrom::Start(0))?;
    let mut bytes = [0_u8; 8];
    file.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn gzip_decompress(block: &[u8]) -> io::Result<Vec<u8>> {
    let mut output = vec![0_u8; 65_536];
    let mut stream: ZStream = unsafe { zeroed() };
    stream.next_in = block.as_ptr() as *mut u8;
    stream.avail_in = block.len() as c_uint;
    stream.next_out = output.as_mut_ptr();
    stream.avail_out = output.len() as c_uint;
    let initialized = unsafe {
        inflateInit2_(
            &mut stream,
            31,
            zlibVersion(),
            size_of::<ZStream>() as c_int,
        )
    };
    if initialized != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Could not initialize zlib",
        ));
    }
    let result = unsafe { inflate(&mut stream, Z_FINISH) };
    unsafe { inflateEnd(&mut stream) };
    if result != Z_STREAM_END {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Could not decompress BGZF block",
        ));
    }
    output.truncate(stream.total_out as usize);
    Ok(output)
}

struct BgzfReader {
    file: File,
    block_address: u64,
    block_size: u64,
    block: Vec<u8>,
    position: usize,
}

impl BgzfReader {
    fn new(path: &str) -> io::Result<Self> {
        Ok(Self {
            file: File::open(path)?,
            block_address: 0,
            block_size: 0,
            block: Vec::new(),
            position: 0,
        })
    }

    fn load_block(&mut self, address: u64) -> io::Result<()> {
        self.file.seek(SeekFrom::Start(address))?;
        let mut header = [0_u8; 12];
        self.file.read_exact(&mut header)?;
        if header[..4] != [31, 139, 8, 4] {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid BGZF block",
            ));
        }
        let extra_length = u16::from_le_bytes(header[10..12].try_into().unwrap()) as usize;
        let mut extra = vec![0_u8; extra_length];
        self.file.read_exact(&mut extra)?;
        let mut position = 0;
        let mut block_size = None;
        while position + 4 <= extra.len() {
            let field_length =
                u16::from_le_bytes(extra[position + 2..position + 4].try_into().unwrap()) as usize;
            if position + 4 + field_length > extra.len() {
                break;
            }
            if &extra[position..position + 2] == b"BC" && field_length == 2 {
                block_size = Some(
                    u16::from_le_bytes(extra[position + 4..position + 6].try_into().unwrap())
                        as u64
                        + 1,
                );
                break;
            }
            position += 4 + field_length;
        }
        let block_size = block_size
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing BGZF BC field"))?;
        self.file.seek(SeekFrom::Start(address))?;
        let mut compressed = vec![0_u8; block_size as usize];
        self.file.read_exact(&mut compressed)?;
        self.block = gzip_decompress(&compressed)?;
        self.block_address = address;
        self.block_size = block_size;
        self.position = 0;
        Ok(())
    }

    fn seek_virtual(&mut self, offset: u64) -> io::Result<()> {
        self.load_block(offset >> 16)?;
        self.position = (offset & 0xffff) as usize;
        if self.position > self.block.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid virtual offset",
            ));
        }
        Ok(())
    }

    fn read_raw(&mut self, count: usize) -> io::Result<Vec<u8>> {
        let mut result = Vec::with_capacity(count);
        while result.len() < count {
            if self.position == self.block.len() {
                self.load_block(self.block_address + self.block_size)?;
            }
            let available = (count - result.len()).min(self.block.len() - self.position);
            if available == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Empty BGZF block",
                ));
            }
            result.extend_from_slice(&self.block[self.position..self.position + available]);
            self.position += available;
        }
        Ok(result)
    }

    fn read_sequence(&mut self, record: &Record) -> io::Result<Vec<u8>> {
        if record.line_bases == 0 || record.line_width < record.line_bases {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid FASTA line width",
            ));
        }
        let mut sequence = Vec::with_capacity(record.sequence_length as usize);
        let mut remaining = record.sequence_length;
        while remaining > 0 {
            let bases = remaining.min(record.line_bases) as usize;
            sequence.extend_from_slice(&self.read_raw(bases)?);
            remaining -= bases as u64;
            if remaining > 0 {
                self.read_raw((record.line_width - record.line_bases) as usize)?;
            }
        }
        Ok(sequence)
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let arguments = parse_arguments().map_err(io::Error::other)?;
    let queries = read_queries(&arguments)?;
    if queries.is_empty() {
        return Err(io::Error::other("Provide at least one ID or --ids-file").into());
    }
    if !Path::new(&arguments.ffx).exists() {
        if arguments.short_id {
            return Err(io::Error::other(format!(
                "Short-ID lookup requires a prebuilt .ffx with short-ID keys. Run: fastars index --fai {} --short-id-regex 'YOUR_REGEX'",
                arguments.fai
            ))
            .into());
        }
        eprintln!("[INFO] Building lookup index {}", arguments.ffx);
        build_lookup_index(&arguments.fai, &arguments.ffx, Some("/tmp"), None)?;
    }
    let mut fai = File::open(&arguments.fai)?;
    let mut ffx = File::open(&arguments.ffx)?;
    validate_lookup_index(&mut ffx)?;
    let mut gzi = File::open(&arguments.gzi)?;
    let gzi_count = read_gzi_count(&mut gzi)?;
    let mut records = Vec::new();
    let mut missing = Vec::new();
    for query in queries {
        let (id, matches) = match query {
            Query::Short(id) => {
                let matches = find_records(&mut fai, &mut ffx, &id, IndexKeyKind::Short)?;
                (id, matches)
            }
            Query::Full(id) => {
                let matches = find_full(&mut fai, &mut ffx, &id)?;
                (id, matches)
            }
        };
        if matches.is_empty() {
            missing.push(id);
        }
        for mut record in matches {
            record.virtual_offset =
                gzi_virtual_offset(&mut gzi, gzi_count, record.sequence_offset)?;
            records.push(record);
        }
    }
    if !missing.is_empty() {
        eprintln!("[WARN] {} IDs were not found", missing.len());
        if arguments.verbose_missing {
            for id in missing {
                eprintln!("[WARN] ID not found: {id}");
            }
        }
    }
    if arguments.sort_by_offset {
        records.sort_by_key(|record| record.virtual_offset);
    }

    let mut reader = BgzfReader::new(&arguments.fasta)?;
    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());
    for record in records {
        reader.seek_virtual(record.virtual_offset)?;
        let sequence = reader.read_sequence(&record)?;
        writeln!(output, ">{}", record.full_id)?;
        for line in sequence.chunks(60) {
            output.write_all(line)?;
            output.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn main() {
    let first_argument = env::args().nth(1);
    if matches!(first_argument.as_deref(), Some("-h" | "--help")) {
        println!("{}", usage());
        return;
    }

    let result = if first_argument.as_deref() == Some("index") {
        run_index_command()
    } else {
        run()
    };
    if let Err(error) = result {
        eprintln!("[ERROR] {error}");
        std::process::exit(1);
    }
}
