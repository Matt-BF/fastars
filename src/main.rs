mod bgzf;
mod fasta;
mod ffx;
mod indexer;
mod progress;

use crate::bgzf::BgzfWriter;
use crate::fasta::FastaReader;
use crate::ffx::{FfxIndex, FfxRecord};
use crate::indexer::{
    DEFAULT_SORT_MEMORY_MIB, build_index_from_fai_gzi, build_index_from_fasta, default_threads,
};
use regex::RegexBuilder;
use std::collections::HashSet;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

const MAX_FETCH_THREADS: usize = 16;

#[derive(Clone, Copy)]
enum IdMode {
    Exact,
    Prefix,
    Partial,
}

struct FetchArgs {
    fasta: String,
    ffx: String,
    ids_file: Option<String>,
    sort_by_offset: bool,
    verbose_missing: bool,
    id_mode: IdMode,
    ignore_case: bool,
    id_regexp: Option<String>,
    invert_match: bool,
    temp_directory: Option<String>,
    sort_memory_bytes: usize,
    threads: usize,
    show_progress: bool,
    fetch_threads: usize,
    ids: Vec<String>,
}

struct IndexArgs {
    fasta: Option<String>,
    fai: Option<String>,
    gzi: Option<String>,
    output: Option<String>,
    temp_directory: Option<String>,
    sort_memory_bytes: usize,
    threads: usize,
    show_progress: bool,
}

struct GenerateArgs {
    output: String,
    number_sequences: u64,
    number_bases: u64,
}

fn usage() -> &'static str {
    r#"fastars — fast random retrieval from BGZF-compressed or plain FASTA files

USAGE:
    fastars fetch --fasta <FASTA.bgz> [OPTIONS] [ID ...]
    fastars index --fasta <FASTA.bgz> [OPTIONS]
    fastars index --fai <FASTA.bgz.fai> --gzi <FASTA.bgz.gzi> [OPTIONS]
    fastars generate --output <FASTA.bgz> --number-seqs <N> --number-bases <N>

COMMANDS:
    fetch    Fetch FASTA records by exact ID, prefix, partial text, or regular expression.
    index    Build a compressed, self-contained .ffx index.
    generate Generate a deterministic BGZF FASTA benchmark fixture.

FETCH OPTIONS:
    --fasta <FILE>             BGZF-compressed or plain FASTA. Required.
    --ffx <FILE>               Self-contained index. Default: <FASTA>.ffx
    -f, --ids-file <FILE>      Read one query ID per line.
    -m, --id-mode <MODE>       Query mode: exact, prefix, or partial. Default: exact
    -i, --ignore-case          ASCII case-insensitive ID and regex matching.
    -r, --id-regexp <REGEX>    Select indexed full IDs matching this regex.
    -v, --invert-match         With --id-regexp, select IDs that do not match.
    -s, --sort-by-offset       Fetch in FASTA order instead of request/index order.
    --verbose-missing          Print every ID query with no matches.
    --temp-directory <DIR>     Directory for temporary files if auto-indexing.
    --sort-memory <MiB>        Auto-index sort memory budget. Default: 512
    --threads <N>              Auto-index worker threads. Default: available CPUs
    --no-progress              Disable indexing progress on stderr.
    --fetch-threads <N>        Fetch explicit query results with up to 16 workers.
    -h, --help                 Print this help message.

INDEX OPTIONS:
    --fasta <FILE>             Build .ffx by scanning a BGZF or plain FASTA.
    --fai <FILE>               Optional build source: existing FASTA .fai.
    --gzi <FILE>               Required with --fai to convert offsets once.
    --output <FILE>            Output index. Default: <FASTA>.ffx or <FAI>.ffx
    --temp-directory <DIR>     Directory for external sort temporary files.
    --sort-memory <MiB>        Memory budget before external sort. Default: 512
    --threads <N>              Index worker threads. Default: available CPUs
    --no-progress              Disable progress output on stderr.
    -h, --help                 Print index-specific help.

The .ffx stores full IDs, FASTA offsets, sequence lengths, and line layout.
Fetching needs only the FASTA and .ffx; .fai/.gzi files are not
read during fetch. Use --id-mode prefix for prefix queries, --id-mode partial
for literal substring queries, or --id-regexp for regular-expression queries.
"#
}

fn generate_usage() -> &'static str {
    r#"fastars generate — create a deterministic BGZF FASTA benchmark fixture

USAGE:
    fastars generate --output <FASTA.bgz> --number-seqs <N> --number-bases <N>

OPTIONS:
    -o, --output <FILE>        BGZF FASTA output path. Required.
    --ns, --number-seqs <N>    Number of sequences. Required.
    --nbp, --number-bases <N>  Bases in every sequence. Required.
    -h, --help                 Print this help message.

Records are named seqid_1 through seqid_N. Sequence bases are deterministic.
"#
}

fn index_usage() -> &'static str {
    r#"fastars index — build a compressed, self-contained lookup index

USAGE:
    fastars index --fasta <FASTA.bgz> [OPTIONS]
    fastars index --fai <FASTA.bgz.fai> --gzi <FASTA.bgz.gzi> [OPTIONS]

OPTIONS:
    --fasta <FILE>             Build by scanning a BGZF or plain FASTA.
    --fai <FILE>               Build from an existing .fai file.
    --gzi <FILE>               BGZF .gzi file. Required with --fai.
    --output <FILE>            Output index. Default: <FASTA>.ffx or <FAI>.ffx
    --temp-directory <DIR>     Directory for external sort temporary files.
    --sort-memory <MiB>        Memory budget before external sort. Default: 512
    --threads <N>              Index worker threads. Default: available CPUs
    --no-progress              Disable progress output on stderr.
    -h, --help                 Print this help message.

The --fasta path is pure Rust and does not need samtools indexes. The --fai
and --gzi path is an accelerator when those files already exist; the resulting
.ffx is still self-contained for fetch. Use --id-mode prefix to fetch IDs by
literal prefix, or --id-regexp to select indexed IDs with a regular expression.
"#
}

fn default_ffx_path(fasta: &str) -> String {
    format!("{fasta}.ffx")
}

fn normalize_id(value: &str) -> String {
    value
        .trim()
        .strip_prefix('>')
        .unwrap_or(value.trim())
        .trim()
        .to_string()
}

fn parse_id_mode(value: &str) -> Result<IdMode, String> {
    match value {
        "exact" => Ok(IdMode::Exact),
        "prefix" => Ok(IdMode::Prefix),
        "partial" => Ok(IdMode::Partial),
        _ => Err(format!(
            "Invalid --id-mode value: {value}. Expected exact, prefix, or partial"
        )),
    }
}

fn parse_fetch_args() -> Result<FetchArgs, String> {
    let mut fasta = None;
    let mut ffx = None;
    let mut ids_file = None;
    let mut sort_by_offset = false;
    let mut verbose_missing = false;
    let mut id_mode = IdMode::Exact;
    let mut ignore_case = false;
    let mut id_regexp = None;
    let mut invert_match = false;
    let mut temp_directory = None;
    let mut sort_memory_mib = DEFAULT_SORT_MEMORY_MIB;
    let mut threads = default_threads();
    let mut show_progress = true;
    let mut fetch_threads = 1;
    let mut ids = Vec::new();
    let mut arguments = env::args().skip(2);

    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--fasta" => fasta = arguments.next(),
            "--ffx" => ffx = arguments.next(),
            "-f" | "--ids-file" => ids_file = arguments.next(),
            "-m" | "--id-mode" => {
                let value = arguments
                    .next()
                    .ok_or_else(|| "--id-mode needs a value".to_string())?;
                id_mode = parse_id_mode(&value)?;
            }
            "-i" | "--ignore-case" => ignore_case = true,
            "-r" | "--id-regexp" => id_regexp = arguments.next(),
            "-v" | "--invert-match" => invert_match = true,
            "-s" | "--sort-by-offset" => sort_by_offset = true,
            "--verbose-missing" => verbose_missing = true,
            "--temp-directory" => temp_directory = arguments.next(),
            "--sort-memory" => {
                let value = arguments
                    .next()
                    .ok_or_else(|| "--sort-memory needs a MiB value".to_string())?;
                sort_memory_mib = value
                    .parse()
                    .map_err(|_| format!("Invalid --sort-memory value: {value}"))?;
                if sort_memory_mib == 0 {
                    return Err("--sort-memory must be greater than zero".to_string());
                }
            }
            "--threads" => {
                let value = arguments
                    .next()
                    .ok_or_else(|| "--threads needs a value".to_string())?;
                threads = value
                    .parse()
                    .map_err(|_| format!("Invalid --threads value: {value}"))?;
                if threads == 0 {
                    return Err("--threads must be greater than zero".to_string());
                }
            }
            "--no-progress" => show_progress = false,
            "--fetch-threads" => {
                let value = arguments
                    .next()
                    .ok_or_else(|| "--fetch-threads needs a value".to_string())?;
                fetch_threads = value
                    .parse()
                    .ok()
                    .filter(|value: &usize| (1..=MAX_FETCH_THREADS).contains(value))
                    .ok_or_else(|| {
                        format!("--fetch-threads must be an integer from 1 to {MAX_FETCH_THREADS}")
                    })?;
            }
            "--full-id" => id_mode = IdMode::Exact,
            "-h" | "--help" => return Err(usage().to_string()),
            _ if argument.starts_with('-') => {
                return Err(format!("Unknown option: {argument}\n\n{}", usage()));
            }
            _ => ids.push(argument),
        }
    }

    let fasta = fasta.ok_or_else(|| format!("--fasta is required\n\n{}", usage()))?;
    let sort_memory_bytes = sort_memory_mib
        .checked_mul(1024 * 1024)
        .ok_or_else(|| "--sort-memory value is too large".to_string())?;
    Ok(FetchArgs {
        ffx: ffx.unwrap_or_else(|| default_ffx_path(&fasta)),
        fasta,
        ids_file,
        sort_by_offset,
        verbose_missing,
        id_mode,
        ignore_case,
        id_regexp,
        invert_match,
        temp_directory,
        sort_memory_bytes,
        threads,
        show_progress,
        fetch_threads,
        ids,
    })
}

fn parse_index_args() -> Result<IndexArgs, String> {
    let mut fasta = None;
    let mut fai = None;
    let mut gzi = None;
    let mut output = None;
    let mut temp_directory = None;
    let mut sort_memory_mib = DEFAULT_SORT_MEMORY_MIB;
    let mut threads = default_threads();
    let mut show_progress = true;
    let mut arguments = env::args().skip(2);

    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--fasta" => fasta = arguments.next(),
            "--fai" => fai = arguments.next(),
            "--gzi" => gzi = arguments.next(),
            "--output" => output = arguments.next(),
            "--temp-directory" => temp_directory = arguments.next(),
            "--sort-memory" => {
                let value = arguments
                    .next()
                    .ok_or_else(|| "--sort-memory needs a MiB value".to_string())?;
                sort_memory_mib = value
                    .parse()
                    .map_err(|_| format!("Invalid --sort-memory value: {value}"))?;
                if sort_memory_mib == 0 {
                    return Err("--sort-memory must be greater than zero".to_string());
                }
            }
            "--threads" => {
                let value = arguments
                    .next()
                    .ok_or_else(|| "--threads needs a value".to_string())?;
                threads = value
                    .parse()
                    .map_err(|_| format!("Invalid --threads value: {value}"))?;
                if threads == 0 {
                    return Err("--threads must be greater than zero".to_string());
                }
            }
            "--no-progress" => show_progress = false,
            "-h" | "--help" => {
                println!("{}", index_usage());
                std::process::exit(0);
            }
            _ => {
                return Err(format!(
                    "Unknown index option: {argument}\n\n{}",
                    index_usage()
                ));
            }
        }
    }

    let sort_memory_bytes = sort_memory_mib
        .checked_mul(1024 * 1024)
        .ok_or_else(|| "--sort-memory value is too large".to_string())?;
    Ok(IndexArgs {
        fasta,
        fai,
        gzi,
        output,
        temp_directory,
        sort_memory_bytes,
        threads,
        show_progress,
    })
}

fn parse_generate_args() -> Result<GenerateArgs, String> {
    let mut output = None;
    let mut number_sequences = None;
    let mut number_bases = None;
    let mut arguments = env::args().skip(2);

    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "-o" | "--output" => output = arguments.next(),
            "--ns" | "--number-seqs" => number_sequences = arguments.next(),
            "--nbp" | "--number-bases" => number_bases = arguments.next(),
            "-h" | "--help" => {
                println!("{}", generate_usage());
                std::process::exit(0);
            }
            _ => {
                return Err(format!(
                    "Unknown generate option: {argument}\n\n{}",
                    generate_usage()
                ));
            }
        }
    }

    Ok(GenerateArgs {
        output: output.ok_or_else(|| format!("--output is required\n\n{}", generate_usage()))?,
        number_sequences: parse_positive_count(
            number_sequences,
            "--number-seqs",
            generate_usage(),
        )?,
        number_bases: parse_positive_count(number_bases, "--number-bases", generate_usage())?,
    })
}

fn parse_positive_count(value: Option<String>, option: &str, help: &str) -> Result<u64, String> {
    let value = value.ok_or_else(|| format!("{option} is required\n\n{help}"))?;
    let value = value
        .parse::<u64>()
        .map_err(|_| format!("{option} must be a positive integer\n\n{help}"))?;
    if value == 0 {
        return Err(format!("{option} must be a positive integer\n\n{help}"));
    }
    Ok(value)
}

fn read_query_ids(arguments: &FetchArgs) -> io::Result<Vec<String>> {
    let mut ids = arguments
        .ids
        .iter()
        .map(|value| normalize_query(value, arguments.id_mode))
        .filter(|value| !value.is_empty())
        .collect::<Vec<_>>();

    if let Some(path) = &arguments.ids_file {
        let reader = BufReader::new(File::open(path)?);
        for line in reader.lines() {
            let id = normalize_query(&line?, arguments.id_mode);
            if !id.is_empty() {
                ids.push(id);
            }
        }
    }

    Ok(ids)
}

fn normalize_query(value: &str, mode: IdMode) -> String {
    let value = normalize_id(value);
    match mode {
        IdMode::Exact | IdMode::Prefix => value,
        IdMode::Partial => value.trim().to_string(),
    }
}

fn add_records(records: Vec<FfxRecord>, seen: &mut HashSet<u64>, output: &mut Vec<FfxRecord>) {
    for record in records {
        if seen.insert(record.record_id) {
            output.push(record);
        }
    }
}

fn record_bytes(reader: &mut FastaReader, record: &FfxRecord) -> io::Result<Vec<u8>> {
    reader.seek(record.virtual_offset)?;
    let sequence =
        reader.read_sequence(record.sequence_length, record.line_bases, record.line_width)?;
    let mut output = Vec::with_capacity(
        record
            .header
            .len()
            .saturating_add(sequence.len())
            .saturating_add(sequence.len().div_ceil(60))
            .saturating_add(2),
    );
    writeln!(output, ">{}", record.header)?;
    for line in sequence.chunks(60) {
        output.extend_from_slice(line);
        output.push(b'\n');
    }
    Ok(output)
}

fn write_record<W: Write>(
    reader: &mut FastaReader,
    output: &mut W,
    record: &FfxRecord,
) -> io::Result<()> {
    output.write_all(&record_bytes(reader, record)?)
}

fn write_records<W: Write>(
    fasta_path: &str,
    output: &mut W,
    records: &[FfxRecord],
    fetch_threads: usize,
) -> io::Result<()> {
    if fetch_threads == 1 || records.len() < 2 {
        let mut reader = FastaReader::new(fasta_path)?;
        for record in records {
            write_record(&mut reader, output, record)?;
        }
        return Ok(());
    }

    let batch_size = fetch_threads.saturating_mul(2);
    for batch in records.chunks(batch_size) {
        let worker_count = fetch_threads.min(batch.len());
        let fetched = std::thread::scope(|scope| -> io::Result<Vec<(usize, Vec<u8>)>> {
            let mut workers = Vec::with_capacity(worker_count);
            for worker_index in 0..worker_count {
                workers.push(scope.spawn(move || -> io::Result<Vec<(usize, Vec<u8>)>> {
                    let mut reader = FastaReader::new(fasta_path)?;
                    batch
                        .iter()
                        .enumerate()
                        .skip(worker_index)
                        .step_by(worker_count)
                        .map(|(record_index, record)| {
                            record_bytes(&mut reader, record).map(|bytes| (record_index, bytes))
                        })
                        .collect()
                }));
            }

            let mut fetched = Vec::with_capacity(batch.len());
            for worker in workers {
                fetched.extend(
                    worker
                        .join()
                        .map_err(|_| io::Error::other("FASTA fetch worker panicked"))??,
                );
            }
            Ok(fetched)
        })?;
        let mut fetched = fetched;
        fetched.sort_unstable_by_key(|(record_index, _)| *record_index);
        for (_, bytes) in fetched {
            output.write_all(&bytes)?;
        }
    }
    Ok(())
}

fn run_index_command() -> Result<(), Box<dyn Error>> {
    let arguments = parse_index_args().map_err(io::Error::other)?;
    let output = arguments.output.clone().unwrap_or_else(|| {
        arguments
            .fasta
            .as_ref()
            .map(|path| format!("{path}.ffx"))
            .unwrap_or_else(|| format!("{}.ffx", arguments.fai.as_deref().unwrap_or("fastars")))
    });

    match (arguments.fai.as_deref(), arguments.gzi.as_deref()) {
        (Some(fai), Some(gzi)) => {
            build_index_from_fai_gzi(
                fai,
                gzi,
                &output,
                arguments.temp_directory.as_deref(),
                arguments.sort_memory_bytes,
                arguments.threads,
                arguments.show_progress,
            )?;
        }
        (Some(_), None) | (None, Some(_)) => {
            return Err(io::Error::other("index requires both --fai and --gzi, or --fasta").into());
        }
        (None, None) => {
            let fasta = arguments
                .fasta
                .as_deref()
                .ok_or_else(|| io::Error::other("index requires --fasta, or --fai with --gzi"))?;
            build_index_from_fasta(
                fasta,
                &output,
                arguments.temp_directory.as_deref(),
                arguments.sort_memory_bytes,
                arguments.threads,
                arguments.show_progress,
            )?;
        }
    }

    eprintln!("Wrote {output}");
    Ok(())
}

fn run_generate_command() -> Result<(), Box<dyn Error>> {
    let arguments = parse_generate_args().map_err(io::Error::other)?;
    generate_fasta(
        &arguments.output,
        arguments.number_sequences,
        arguments.number_bases,
    )?;
    eprintln!("Wrote {}", arguments.output);
    Ok(())
}

fn generate_fasta(output: &str, number_sequences: u64, number_bases: u64) -> io::Result<()> {
    let mut writer = BgzfWriter::new(output)?;
    for sequence_number in 1..=number_sequences {
        writeln!(writer, ">seqid_{sequence_number}")?;
        write_generated_sequence(&mut writer, sequence_number, number_bases)?;
    }
    writer.finish()
}

fn write_generated_sequence(
    writer: &mut BgzfWriter,
    sequence_number: u64,
    number_bases: u64,
) -> io::Result<()> {
    let mut bases = [b'A'; 8_192];
    let mut state = sequence_number;
    let mut remaining = number_bases;
    while remaining > 0 {
        let count = remaining.min(bases.len() as u64) as usize;
        for base in &mut bases[..count] {
            *base = b"ACGT"[next_random(&mut state) as usize];
        }
        writer.write_all(&bases[..count])?;
        remaining -= count as u64;
    }
    writer.write_all(b"\n")
}

fn next_random(state: &mut u64) -> u8 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut value = *state;
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    ((value ^ (value >> 31)) & 3) as u8
}

fn run_fetch_command() -> Result<(), Box<dyn Error>> {
    let arguments = parse_fetch_args().map_err(io::Error::other)?;
    let query_ids = read_query_ids(&arguments)?;

    if query_ids.is_empty() && arguments.id_regexp.is_none() {
        return Err(io::Error::other("Provide at least one ID, --ids-file, or --id-regexp").into());
    }
    if arguments.invert_match && arguments.id_regexp.is_none() {
        return Err(io::Error::other("--invert-match requires --id-regexp").into());
    }

    if !Path::new(&arguments.ffx).exists() {
        eprintln!(
            "[INFO] {} not found; building it from {}",
            arguments.ffx, arguments.fasta
        );
        build_index_from_fasta(
            &arguments.fasta,
            &arguments.ffx,
            arguments.temp_directory.as_deref(),
            arguments.sort_memory_bytes,
            arguments.threads,
            arguments.show_progress,
        )?;
    }

    let mut index = FfxIndex::open(&arguments.ffx)?;
    let mut records = Vec::new();
    let mut seen = HashSet::new();
    let mut missing = Vec::new();

    let query_matches = match (arguments.id_mode, arguments.ignore_case) {
        (IdMode::Exact, false) => query_ids
            .iter()
            .map(|query| index.find_exact(query))
            .collect::<io::Result<Vec<_>>>()?,
        (IdMode::Prefix, false) => query_ids
            .iter()
            .map(|query| index.find_prefix(query))
            .collect::<io::Result<Vec<_>>>()?,
        (IdMode::Partial, ignore_case) => index.find_partial_many(&query_ids, ignore_case)?,
        (IdMode::Exact, true) => index.find_scan_many(&query_ids, |record, query| {
            if query.bytes().any(|byte| byte.is_ascii_whitespace()) {
                record.header.eq_ignore_ascii_case(query)
            } else {
                record.full_id.eq_ignore_ascii_case(query)
            }
        })?,
        (IdMode::Prefix, true) => index.find_scan_many(&query_ids, |record, query| {
            let candidate = if query.bytes().any(|byte| byte.is_ascii_whitespace()) {
                &record.header
            } else {
                &record.full_id
            };
            candidate
                .get(..query.len())
                .is_some_and(|prefix| prefix.eq_ignore_ascii_case(query))
        })?,
    };

    for (id, matches) in query_ids.into_iter().zip(query_matches) {
        if matches.is_empty() {
            missing.push(id);
        } else {
            add_records(matches, &mut seen, &mut records);
        }
    }

    if !missing.is_empty() {
        eprintln!("[WARN] {} ID queries had no matches", missing.len());
        if arguments.verbose_missing {
            for id in missing {
                eprintln!("[WARN] ID not found: {id}");
            }
        }
    }

    let regex = arguments
        .id_regexp
        .as_deref()
        .map(|pattern| {
            RegexBuilder::new(pattern)
                .case_insensitive(arguments.ignore_case)
                .build()
        })
        .transpose()?;
    if arguments.sort_by_offset {
        if let Some(regex) = &regex {
            let matches = index.find_regex(regex, arguments.invert_match)?;
            add_records(matches, &mut seen, &mut records);
        }
        records.sort_by_key(|record| record.virtual_offset);
    }

    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());
    write_records(
        &arguments.fasta,
        &mut output,
        &records,
        arguments.fetch_threads,
    )?;

    if !arguments.sort_by_offset
        && let Some(regex) = &regex
    {
        let mut reader = FastaReader::new(&arguments.fasta)?;
        index.for_each_regex(regex, arguments.invert_match, |record| {
            if seen.insert(record.record_id) {
                write_record(&mut reader, &mut output, &record)?;
            }
            Ok(())
        })?;
    }

    Ok(())
}

fn main() {
    let first_argument = env::args().nth(1);
    let second_argument = env::args().nth(2);
    if matches!(first_argument.as_deref(), Some("-h" | "--help"))
        || (first_argument.as_deref() == Some("fetch")
            && matches!(second_argument.as_deref(), Some("-h" | "--help")))
    {
        println!("{}", usage());
        return;
    }

    let result = match first_argument.as_deref() {
        Some("fetch") => run_fetch_command(),
        Some("index") => run_index_command(),
        Some("generate") => run_generate_command(),
        Some(command) => {
            Err(io::Error::other(format!("Unknown command: {command}\n\n{}", usage())).into())
        }
        None => Err(io::Error::other(format!("A command is required\n\n{}", usage())).into()),
    };

    if let Err(error) = result {
        eprintln!("[ERROR] {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_PATH: AtomicU64 = AtomicU64::new(0);

    fn test_path(suffix: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "fastars-generate-{}-{}-{suffix}",
            std::process::id(),
            NEXT_PATH.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn query_normalization_preserves_complete_headers() {
        assert_eq!(normalize_id("seq1 description"), "seq1 description");
        assert_eq!(normalize_id(">seq1 description\n"), "seq1 description");
        assert_eq!(normalize_id("  seq1  "), "seq1");
        assert_eq!(
            normalize_query("  >description with spaces  ", IdMode::Partial),
            "description with spaces"
        );
    }

    #[test]
    fn generated_fasta_is_indexable() {
        let fasta = test_path("fixture.fna.bgz");
        let index_path = test_path("fixture.fna.bgz.ffx");
        generate_fasta(fasta.to_str().unwrap(), 3, 65_281).unwrap();
        build_index_from_fasta(
            fasta.to_str().unwrap(),
            index_path.to_str().unwrap(),
            None,
            1024 * 1024,
            2,
            false,
        )
        .unwrap();

        let mut index = FfxIndex::open(index_path.to_str().unwrap()).unwrap();
        let record = index.find_exact("seqid_2").unwrap().pop().unwrap();
        assert_eq!(record.sequence_length, 65_281);

        let mut reader = FastaReader::new(&fasta).unwrap();
        reader.seek(record.virtual_offset).unwrap();
        let sequence = reader
            .read_sequence(record.sequence_length, record.line_bases, record.line_width)
            .unwrap();
        assert_eq!(sequence.len(), 65_281);
        for base in b"ACGT" {
            assert!(sequence.contains(base));
        }

        fs::remove_file(fasta).unwrap();
        fs::remove_file(index_path).unwrap();
    }

    #[test]
    fn parallel_fetch_preserves_requested_order() {
        let fasta = test_path("parallel-fetch.fna.bgz");
        let index_path = test_path("parallel-fetch.fna.bgz.ffx");
        generate_fasta(fasta.to_str().unwrap(), 8, 100).unwrap();
        build_index_from_fasta(
            fasta.to_str().unwrap(),
            index_path.to_str().unwrap(),
            None,
            1024 * 1024,
            2,
            false,
        )
        .unwrap();

        let mut index = FfxIndex::open(index_path.to_str().unwrap()).unwrap();
        let records = ["seqid_8", "seqid_2", "seqid_6", "seqid_1"]
            .iter()
            .map(|id| index.find_exact(id).unwrap().pop().unwrap())
            .collect::<Vec<_>>();
        let mut sequential = Vec::new();
        write_records(fasta.to_str().unwrap(), &mut sequential, &records, 1).unwrap();
        let mut parallel = Vec::new();
        write_records(fasta.to_str().unwrap(), &mut parallel, &records, 3).unwrap();
        assert_eq!(parallel, sequential);

        fs::remove_file(fasta).unwrap();
        fs::remove_file(index_path).unwrap();
    }
}
