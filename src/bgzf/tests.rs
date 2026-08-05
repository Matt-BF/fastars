use super::*;
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_PATH: AtomicU64 = AtomicU64::new(0);

fn test_path() -> PathBuf {
    std::env::temp_dir().join(format!(
        "fastars-bgzf-{}-{}.bgz",
        std::process::id(),
        NEXT_PATH.fetch_add(1, Ordering::Relaxed)
    ))
}

fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = u32::MAX;
    for byte in bytes {
        crc ^= u32::from(*byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xedb8_8320 & (0_u32.wrapping_sub(crc & 1)));
        }
    }
    !crc
}

fn uncompressed_bgzf_block(payload: &[u8]) -> Vec<u8> {
    assert!(payload.len() <= u16::MAX as usize);
    let block_size = 18 + 5 + payload.len() + 8;
    let bsize = (block_size - 1) as u16;
    let length = payload.len() as u16;
    let mut block = vec![
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
        b'B',
        b'C',
        2,
        0,
        bsize as u8,
        (bsize >> 8) as u8,
        1,
        length as u8,
        (length >> 8) as u8,
        !length as u8,
        (!length >> 8) as u8,
    ];
    block.extend_from_slice(payload);
    block.extend_from_slice(&crc32(payload).to_le_bytes());
    block.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    assert_eq!(block.len(), block_size);
    block
}

fn read_all_blocks(mut reader: BgzfReader) -> io::Result<Vec<(u64, Vec<u8>)>> {
    let mut blocks = Vec::new();
    while let Some(block) = reader.read_block()? {
        blocks.push((block.address, block.bytes));
    }
    Ok(blocks)
}

#[test]
fn parallel_reader_matches_serial_blocks() {
    let path = test_path();
    let mut contents = uncompressed_bgzf_block(b">alpha\nACGT\n>be");
    contents.extend_from_slice(&uncompressed_bgzf_block(b"ta\nTT\n"));
    fs::write(&path, contents).unwrap();

    let serial = read_all_blocks(BgzfReader::new(&path).unwrap()).unwrap();
    let single_threaded = read_all_blocks(BgzfReader::with_threads(&path, 1).unwrap()).unwrap();
    let parallel = read_all_blocks(BgzfReader::with_threads(&path, 4).unwrap()).unwrap();
    assert_eq!(single_threaded, serial);
    assert_eq!(parallel, serial);
    assert_eq!(
        parallel
            .iter()
            .map(|(_, bytes)| bytes.as_slice())
            .collect::<Vec<_>>(),
        vec![b">alpha\nACGT\n>be".as_slice(), b"ta\nTT\n".as_slice(),]
    );

    fs::remove_file(path).unwrap();
}
