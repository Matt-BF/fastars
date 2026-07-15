use std::ffi::c_void;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::mem::{size_of, zeroed};
use std::os::raw::{c_char, c_int, c_uint, c_ulong};
use std::path::Path;

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

pub struct BgzfLine {
    pub start_virtual: u64,
    pub bytes: Vec<u8>,
}

pub struct BgzfReader {
    file: File,
    file_len: u64,
    block_address: u64,
    block_size: u64,
    block: Vec<u8>,
    position: usize,
    started: bool,
}

impl BgzfReader {
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        let file_len = file.metadata()?.len();
        Ok(Self {
            file,
            file_len,
            block_address: 0,
            block_size: 0,
            block: Vec::new(),
            position: 0,
            started: false,
        })
    }

    pub fn read_line(&mut self) -> io::Result<Option<BgzfLine>> {
        let Some((start_virtual, first)) = self.next_byte_with_offset()? else {
            return Ok(None);
        };
        let mut bytes = vec![first];
        while *bytes.last().unwrap() != b'\n' {
            let Some((_, byte)) = self.next_byte_with_offset()? else {
                break;
            };
            bytes.push(byte);
        }
        Ok(Some(BgzfLine {
            start_virtual,
            bytes,
        }))
    }

    pub fn seek_virtual(&mut self, offset: u64) -> io::Result<()> {
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

    fn next_byte_with_offset(&mut self) -> io::Result<Option<(u64, u8)>> {
        if !self.ensure_block()? {
            return Ok(None);
        }

        let virtual_offset = (self.block_address << 16) | self.position as u64;
        let byte = self.block[self.position];
        self.position += 1;
        Ok(Some((virtual_offset, byte)))
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
        if address >= self.file_len {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "BGZF block address is past end of file",
            ));
        }

        self.file.seek(SeekFrom::Start(address))?;
        let mut header = [0_u8; 12];
        self.file.read_exact(&mut header)?;
        if header[..4] != [31, 139, 8, 4] {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Input is not BGZF-compressed data",
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
