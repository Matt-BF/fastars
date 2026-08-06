mod bgzf;
mod fasta;
mod ffx;
mod indexer;
mod progress;

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

fn usage() -> &'static str {
    r#"fastars — fast random retrieval from BGZF-compressed or plain FASTA files

USAGE:
    fastars fetch --fasta <FASTA.bgz> [OPTIONS] [ID ...]
    fastars index --fasta <FASTA.bgz> [OPTIONS]
    fastars index --fai <FASTA.bgz.fai> --gzi <FASTA.bgz.gzi> [OPTIONS]

COMMANDS:
    fetch    Fetch FASTA records by exact ID, prefix, partial text, or regular expression.
    index    Build a compressed, self-contained .ffx index.

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

fn write_record(
    reader: &mut FastaReader,
    output: &mut BufWriter<io::StdoutLock<'_>>,
    record: &FfxRecord,
) -> io::Result<()> {
    reader.seek(record.virtual_offset)?;
    let sequence =
        reader.read_sequence(record.sequence_length, record.line_bases, record.line_width)?;
    writeln!(output, ">{}", record.header)?;
    for line in sequence.chunks(60) {
        output.write_all(line)?;
        output.write_all(b"\n")?;
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

    let mut reader = FastaReader::new(&arguments.fasta)?;
    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());
    for record in &records {
        write_record(&mut reader, &mut output, record)?;
    }

    if !arguments.sort_by_offset
        && let Some(regex) = &regex
    {
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
}
