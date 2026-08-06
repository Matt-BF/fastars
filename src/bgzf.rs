mod parallel;
#[cfg(test)]
mod tests;

use parallel::{LoadedBlock, ParallelBlockLoader};
use std::ffi::c_void;
use std::fs::File;
use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::mem::{size_of, zeroed};
use std::os::raw::{c_char, c_int, c_uint, c_ulong};
use std::path::Path;

const Z_STREAM_END: c_int = 1;
const Z_FINISH: c_int = 4;
const Z_DEFLATED: c_int = 8;
const Z_DEFAULT_STRATEGY: c_int = 0;
const BGZF_BLOCK_PAYLOAD_SIZE: usize = 65_280;
const BGZF_EOF_BLOCK: [u8; 28] = [
    31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

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
    fn deflateInit2_(
        stream: *mut ZStream,
        level: c_int,
        method: c_int,
        window_bits: c_int,
        memory_level: c_int,
        strategy: c_int,
        version: *const c_char,
        stream_size: c_int,
    ) -> c_int;
    fn deflate(stream: *mut ZStream, flush: c_int) -> c_int;
    fn deflateEnd(stream: *mut ZStream) -> c_int;
    fn crc32(crc: c_ulong, buffer: *const u8, length: c_uint) -> c_ulong;
}

pub(crate) struct BgzfBlock {
    pub address: u64,
    pub bytes: Vec<u8>,
}

pub struct BgzfReader {
    loader: BlockLoader,
    file_len: u64,
    block_address: u64,
    block_size: u64,
    block: Vec<u8>,
    position: usize,
    started: bool,
}

pub(crate) struct BgzfWriter {
    file: BufWriter<File>,
    buffer: Vec<u8>,
}

impl BgzfWriter {
    pub(crate) fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            file: BufWriter::new(File::create(path)?),
            buffer: Vec::with_capacity(BGZF_BLOCK_PAYLOAD_SIZE),
        })
    }

    pub(crate) fn finish(mut self) -> io::Result<()> {
        self.flush_block()?;
        self.file.write_all(&BGZF_EOF_BLOCK)?;
        self.file.flush()
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let compressed = raw_deflate(&self.buffer)?;
        let block_size = 18_usize
            .checked_add(compressed.len())
            .and_then(|size| size.checked_add(8))
            .and_then(|size| u16::try_from(size).ok())
            .ok_or_else(|| io::Error::other("Generated BGZF block is too large"))?;
        let stored_size = (block_size - 1).to_le_bytes();
        self.file.write_all(&[
            31,
            139,
            8,
            4,
            0,
            0,
            0,
            0,
            0,
            255,
            6,
            0,
            66,
            67,
            2,
            0,
            stored_size[0],
            stored_size[1],
        ])?;
        self.file.write_all(&compressed)?;
        let checksum = unsafe { crc32(0, self.buffer.as_ptr(), self.buffer.len() as c_uint) };
        self.file.write_all(&(checksum as u32).to_le_bytes())?;
        self.file
            .write_all(&(self.buffer.len() as u32).to_le_bytes())?;
        self.buffer.clear();
        Ok(())
    }
}

impl Write for BgzfWriter {
    fn write(&mut self, mut bytes: &[u8]) -> io::Result<usize> {
        let input_len = bytes.len();
        while !bytes.is_empty() {
            let available = BGZF_BLOCK_PAYLOAD_SIZE - self.buffer.len();
            let count = available.min(bytes.len());
            self.buffer.extend_from_slice(&bytes[..count]);
            bytes = &bytes[count..];
            if self.buffer.len() == BGZF_BLOCK_PAYLOAD_SIZE {
                self.flush_block()?;
            }
        }
        Ok(input_len)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_block()?;
        self.file.flush()
    }
}

enum BlockLoader {
    Serial {
        file: File,
        spare_input: Vec<u8>,
        spare_output: Vec<u8>,
    },
    Parallel(ParallelBlockLoader),
}

impl BgzfReader {
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        let file_len = file.metadata()?.len();
        Ok(Self {
            loader: BlockLoader::Serial {
                file,
                spare_input: Vec::new(),
                spare_output: Vec::new(),
            },
            file_len,
            block_address: 0,
            block_size: 0,
            block: Vec::new(),
            position: 0,
            started: false,
        })
    }

    pub fn with_threads<P: AsRef<Path>>(path: P, threads: usize) -> io::Result<Self> {
        if threads == 0 {
            return Err(io::Error::other(
                "BGZF thread count must be greater than zero",
            ));
        }
        if threads == 1 {
            return Self::new(path);
        }
        let file = File::open(path)?;
        let file_len = file.metadata()?.len();
        Ok(Self {
            loader: BlockLoader::Parallel(ParallelBlockLoader::new(file, file_len, threads)),
            file_len,
            block_address: 0,
            block_size: 0,
            block: Vec::new(),
            position: 0,
            started: false,
        })
    }

    pub(crate) fn read_block(&mut self) -> io::Result<Option<BgzfBlock>> {
        let address = if self.started {
            self.block_address
                .checked_add(self.block_size)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "BGZF offset overflow"))?
        } else {
            0
        };
        if address >= self.file_len {
            return Ok(None);
        }

        self.load_block(address)?;
        self.started = true;
        Ok(Some(BgzfBlock {
            address,
            bytes: std::mem::take(&mut self.block),
        }))
    }

    pub(crate) fn total_bytes(&self) -> u64 {
        self.file_len
    }

    pub(crate) fn consumed_bytes(&self) -> u64 {
        if self.started {
            self.block_address
                .saturating_add(self.block_size)
                .min(self.file_len)
        } else {
            0
        }
    }

    pub(crate) fn recycle_block(&mut self, bytes: Vec<u8>) {
        match &mut self.loader {
            BlockLoader::Serial { spare_output, .. } => *spare_output = bytes,
            BlockLoader::Parallel(loader) => loader.recycle_output(bytes),
        }
    }

    pub fn seek_virtual(&mut self, offset: u64) -> io::Result<()> {
        if matches!(self.loader, BlockLoader::Parallel(_)) {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "Parallel BGZF readers only support sequential access",
            ));
        }
        self.load_block(offset >> 16)?;
        self.started = true;
        self.position = (offset & 0xffff) as usize;
        if self.position > self.block.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid BGZF virtual offset",
            ));
        }
        Ok(())
    }

    pub fn read_sequence(
        &mut self,
        sequence_length: u64,
        line_bases: u64,
        line_width: u64,
    ) -> io::Result<Vec<u8>> {
        if line_bases == 0 || line_width < line_bases {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid FASTA line layout in .ffx",
            ));
        }

        let mut sequence = Vec::with_capacity(sequence_length as usize);
        let mut remaining = sequence_length;
        while remaining > 0 {
            let bases = remaining.min(line_bases) as usize;
            sequence.extend_from_slice(&self.read_raw(bases)?);
            remaining -= bases as u64;
            if remaining > 0 {
                self.read_raw((line_width - line_bases) as usize)?;
            }
        }
        Ok(sequence)
    }

    fn read_raw(&mut self, count: usize) -> io::Result<Vec<u8>> {
        let mut result = Vec::with_capacity(count);
        while result.len() < count {
            if !self.ensure_block()? {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Unexpected end of BGZF stream",
                ));
            }
            let available = (count - result.len()).min(self.block.len() - self.position);
            if available == 0 {
                continue;
            }
            result.extend_from_slice(&self.block[self.position..self.position + available]);
            self.position += available;
        }
        Ok(result)
    }

    fn ensure_block(&mut self) -> io::Result<bool> {
        loop {
            if !self.started {
                if self.file_len == 0 {
                    return Ok(false);
                }
                self.load_block(0)?;
                self.started = true;
            } else if self.position >= self.block.len() {
                let next_address = self.block_address + self.block_size;
                if next_address >= self.file_len {
                    return Ok(false);
                }
                self.load_block(next_address)?;
            }

            if self.position < self.block.len() {
                return Ok(true);
            }
        }
    }

    fn load_block(&mut self, address: u64) -> io::Result<()> {
        let loaded = match &mut self.loader {
            BlockLoader::Serial {
                file,
                spare_input,
                spare_output,
            } => {
                let input = std::mem::take(spare_input);
                let (compressed, size) =
                    read_compressed_block(file, self.file_len, address, input)?;
                let output = std::mem::take(spare_output);
                let loaded = LoadedBlock {
                    address,
                    size,
                    bytes: gzip_decompress(&compressed, output)?,
                };
                *spare_input = compressed;
                loaded
            }
            BlockLoader::Parallel(loader) => loader.load(address)?,
        };
        self.block = loaded.bytes;
        self.block_address = loaded.address;
        self.block_size = loaded.size;
        self.position = 0;
        Ok(())
    }
}

fn read_compressed_block(
    file: &mut File,
    file_len: u64,
    address: u64,
    mut compressed: Vec<u8>,
) -> io::Result<(Vec<u8>, u64)> {
    if address >= file_len {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "BGZF block address is past end of file",
        ));
    }

    file.seek(SeekFrom::Start(address))?;
    compressed.resize(12, 0);
    file.read_exact(&mut compressed)?;
    if compressed[..4] != [31, 139, 8, 4] {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Input is not BGZF-compressed data",
        ));
    }

    let extra_length = u16::from_le_bytes(compressed[10..12].try_into().unwrap()) as usize;
    compressed.resize(12 + extra_length, 0);
    file.read_exact(&mut compressed[12..])?;
    let extra = &compressed[12..];

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
                u16::from_le_bytes(extra[position + 4..position + 6].try_into().unwrap()) as u64
                    + 1,
            );
            break;
        }
        position += 4 + field_length;
    }

    let block_size = block_size
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing BGZF BC field"))?;
    let prefix_len = compressed.len();
    let block_len = usize::try_from(block_size)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "BGZF block is too large"))?;
    if block_len < prefix_len || address.saturating_add(block_size) > file_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid BGZF block size",
        ));
    }
    compressed.resize(block_len, 0);
    file.read_exact(&mut compressed[prefix_len..])?;
    Ok((compressed, block_size))
}

fn gzip_decompress(block: &[u8], mut output: Vec<u8>) -> io::Result<Vec<u8>> {
    output.resize(65_536, 0);
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

fn raw_deflate(input: &[u8]) -> io::Result<Vec<u8>> {
    let mut output = vec![0_u8; input.len() + 64];
    let mut stream: ZStream = unsafe { zeroed() };
    stream.next_in = input.as_ptr() as *mut u8;
    stream.avail_in = input.len() as c_uint;
    stream.next_out = output.as_mut_ptr();
    stream.avail_out = output.len() as c_uint;

    let initialized = unsafe {
        deflateInit2_(
            &mut stream,
            6,
            Z_DEFLATED,
            -15,
            8,
            Z_DEFAULT_STRATEGY,
            zlibVersion(),
            size_of::<ZStream>() as c_int,
        )
    };
    if initialized != 0 {
        return Err(io::Error::other("Could not initialize zlib compression"));
    }

    let result = unsafe { deflate(&mut stream, Z_FINISH) };
    unsafe { deflateEnd(&mut stream) };
    if result != Z_STREAM_END {
        return Err(io::Error::other("Could not compress BGZF block"));
    }

    output.truncate(stream.total_out as usize);
    Ok(output)
}
