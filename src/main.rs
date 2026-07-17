mod bgzf;
mod ffx;
mod indexer;

use crate::bgzf::BgzfReader;
use crate::ffx::{FfxIndex, FfxRecord};
use crate::indexer::{build_index_from_fai_gzi, build_index_from_fasta};
use regex::Regex;
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
}

struct FetchArgs {
    fasta: String,
    ffx: String,
    ids_file: Option<String>,
    sort_by_offset: bool,
    verbose_missing: bool,
    id_mode: IdMode,
    id_regexp: Option<String>,
    invert_match: bool,
    temp_directory: Option<String>,
    ids: Vec<String>,
}

struct IndexArgs {
    fasta: Option<String>,
    fai: Option<String>,
    gzi: Option<String>,
    output: Option<String>,
    temp_directory: Option<String>,
}

fn usage() -> &'static str {
    r#"fastars — fast random retrieval from BGZF-compressed FASTA files

USAGE:
    fastars --fasta <FASTA.bgz> [OPTIONS] [ID ...]
    fastars index --fasta <FASTA.bgz> [OPTIONS]
    fastars index --fai <FASTA.bgz.fai> --gzi <FASTA.bgz.gzi> [OPTIONS]

COMMANDS:
    index    Build a self-contained .ffx index.

FETCH OPTIONS:
    --fasta <FILE>             BGZF-compressed FASTA to read. Required.
    --ffx <FILE>               Self-contained index. Default: <FASTA>.ffx
    -f, --ids-file <FILE>      Read one query ID per line.
    -m, --id-mode <MODE>       Query mode for IDs: exact or prefix. Default: exact
    -r, --id-regexp <REGEX>    Select indexed full IDs matching this regex.
    -v, --invert-match         With --id-regexp, select IDs that do not match.
    -s, --sort-by-offset       Fetch in FASTA order instead of request/index order.
    --verbose-missing          Print every exact/prefix query with no matches.
    --temp-directory <DIR>     Directory for temporary files if auto-indexing.
    -h, --help                 Print this help message.

INDEX OPTIONS:
    --fasta <FILE>             Build .ffx by scanning a BGZF FASTA directly.
    --fai <FILE>               Optional build source: existing FASTA .fai.
    --gzi <FILE>               Required with --fai to convert offsets once.
    --output <FILE>            Output index. Default: <FASTA>.ffx or <FAI>.ffx
    --temp-directory <DIR>     Directory for external sort temporary files.
    -h, --help                 Print index-specific help.

The .ffx stores full IDs, BGZF virtual offsets, sequence lengths, and line
layout. Fetching needs only the BGZF FASTA and .ffx; .fai/.gzi files are not
read during fetch. Use --id-mode prefix for prefix queries or --id-regexp for
regular-expression queries.
"#
}

fn index_usage() -> &'static str {
    r#"fastars index — build a self-contained lookup index

USAGE:
    fastars index --fasta <FASTA.bgz> [OPTIONS]
    fastars index --fai <FASTA.bgz.fai> --gzi <FASTA.bgz.gzi> [OPTIONS]

OPTIONS:
    --fasta <FILE>             Build by scanning a BGZF FASTA directly.
    --fai <FILE>               Build from an existing .fai file.
    --gzi <FILE>               BGZF .gzi file. Required with --fai.
    --output <FILE>            Output index. Default: <FASTA>.ffx or <FAI>.ffx
    --temp-directory <DIR>     Directory for external sort temporary files.
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
        .trim_start_matches('>')
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

fn parse_id_mode(value: &str) -> Result<IdMode, String> {
    match value {
        "exact" => Ok(IdMode::Exact),
        "prefix" => Ok(IdMode::Prefix),
        _ => Err(format!(
            "Invalid --id-mode value: {value}. Expected exact or prefix"
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
    let mut id_regexp = None;
    let mut invert_match = false;
    let mut temp_directory = None;
    let mut ids = Vec::new();
    let mut arguments = env::args().skip(1);

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
            "-r" | "--id-regexp" => id_regexp = arguments.next(),
            "-v" | "--invert-match" => invert_match = true,
            "-s" | "--sort-by-offset" => sort_by_offset = true,
            "--verbose-missing" => verbose_missing = true,
            "--temp-directory" => temp_directory = arguments.next(),
            "--full-id" => id_mode = IdMode::Exact,
            "-h" | "--help" => return Err(usage().to_string()),
            _ if argument.starts_with('-') => {
                return Err(format!("Unknown option: {argument}\n\n{}", usage()));
            }
            _ => ids.push(argument),
        }
    }

    let fasta = fasta.ok_or_else(|| format!("--fasta is required\n\n{}", usage()))?;
    Ok(FetchArgs {
        ffx: ffx.unwrap_or_else(|| default_ffx_path(&fasta)),
        fasta,
        ids_file,
        sort_by_offset,
        verbose_missing,
        id_mode,
        id_regexp,
        invert_match,
        temp_directory,
        ids,
    })
}

fn parse_index_args() -> Result<IndexArgs, String> {
    let mut fasta = None;
    let mut fai = None;
    let mut gzi = None;
    let mut output = None;
    let mut temp_directory = None;
    let mut arguments = env::args().skip(2);

    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--fasta" => fasta = arguments.next(),
            "--fai" => fai = arguments.next(),
            "--gzi" => gzi = arguments.next(),
            "--output" => output = arguments.next(),
            "--temp-directory" => temp_directory = arguments.next(),
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

    Ok(IndexArgs {
        fasta,
        fai,
        gzi,
        output,
        temp_directory,
    })
}

fn read_query_ids(arguments: &FetchArgs) -> io::Result<Vec<String>> {
    let mut ids = arguments
        .ids
        .iter()
        .map(|value| normalize_id(value))
        .filter(|value| !value.is_empty())
        .collect::<Vec<_>>();

    if let Some(path) = &arguments.ids_file {
        let reader = BufReader::new(File::open(path)?);
        for line in reader.lines() {
            let id = normalize_id(&line?);
            if !id.is_empty() {
                ids.push(id);
            }
        }
    }

    Ok(ids)
}

fn add_records(records: Vec<FfxRecord>, seen: &mut HashSet<u64>, output: &mut Vec<FfxRecord>) {
    for record in records {
        if seen.insert(record.record_id) {
            output.push(record);
        }
    }
}

fn write_record(
    reader: &mut BgzfReader,
    output: &mut BufWriter<io::StdoutLock<'_>>,
    record: &FfxRecord,
) -> io::Result<()> {
    reader.seek_virtual(record.virtual_offset)?;
    let sequence =
        reader.read_sequence(record.sequence_length, record.line_bases, record.line_width)?;
    writeln!(output, ">{}", record.full_id)?;
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
            build_index_from_fai_gzi(fai, gzi, &output, arguments.temp_directory.as_deref())?;
        }
        (Some(_), None) | (None, Some(_)) => {
            return Err(io::Error::other("index requires both --fai and --gzi, or --fasta").into());
        }
        (None, None) => {
            let fasta = arguments
                .fasta
                .as_deref()
                .ok_or_else(|| io::Error::other("index requires --fasta, or --fai with --gzi"))?;
            build_index_from_fasta(fasta, &output, arguments.temp_directory.as_deref())?;
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
        )?;
    }

    let mut index = FfxIndex::open(&arguments.ffx)?;
    let mut records = Vec::new();
    let mut seen = HashSet::new();
    let mut missing = Vec::new();

    for id in query_ids {
        let matches = match arguments.id_mode {
            IdMode::Exact => index.find_exact(&id)?,
            IdMode::Prefix => index.find_prefix(&id)?,
        };
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

    let regex = arguments.id_regexp.as_deref().map(Regex::new).transpose()?;
    if arguments.sort_by_offset {
        if let Some(regex) = &regex {
            let matches = index.find_regex(regex, arguments.invert_match)?;
            add_records(matches, &mut seen, &mut records);
        }
        records.sort_by_key(|record| record.virtual_offset);
    }

    let mut reader = BgzfReader::new(&arguments.fasta)?;
    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());
    for record in &records {
        write_record(&mut reader, &mut output, record)?;
    }

    if !arguments.sort_by_offset {
        if let Some(regex) = &regex {
            index.for_each_regex(regex, arguments.invert_match, |record| {
                if seen.insert(record.record_id) {
                    write_record(&mut reader, &mut output, &record)?;
                }
                Ok(())
            })?;
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
        run_fetch_command()
    };

    if let Err(error) = result {
        eprintln!("[ERROR] {error}");
        std::process::exit(1);
    }
}
