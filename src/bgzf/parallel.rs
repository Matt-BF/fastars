use super::{gzip_decompress, read_compressed_block};
use std::fs::File;
use std::io;
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread::{self, JoinHandle};

pub(super) struct LoadedBlock {
    pub address: u64,
    pub size: u64,
    pub bytes: Vec<u8>,
}

struct DecompressJob {
    sequence: usize,
    address: u64,
    size: u64,
    compressed: Vec<u8>,
    output: Vec<u8>,
}

struct DecompressResult {
    sequence: usize,
    input: Vec<u8>,
    block: io::Result<LoadedBlock>,
}

pub(super) struct ParallelBlockLoader {
    file: File,
    file_len: u64,
    senders: Vec<Sender<DecompressJob>>,
    results: Receiver<DecompressResult>,
    workers: Vec<JoinHandle<()>>,
    next_worker: usize,
    next_read_address: u64,
    expected_address: u64,
    submitted: usize,
    delivered: usize,
    ready: Vec<DecompressResult>,
    spare_inputs: Vec<Vec<u8>>,
    spare_outputs: Vec<Vec<u8>>,
    max_in_flight: usize,
}

impl ParallelBlockLoader {
    pub fn new(file: File, file_len: u64, threads: usize) -> Self {
        let (result_sender, results) = mpsc::channel();
        let mut senders = Vec::with_capacity(threads);
        let mut workers = Vec::with_capacity(threads);
        for _ in 0..threads {
            let (sender, receiver) = mpsc::channel::<DecompressJob>();
            let result_sender = result_sender.clone();
            workers.push(thread::spawn(move || {
                while let Ok(job) = receiver.recv() {
                    let input = job.compressed;
                    let block = gzip_decompress(&input, job.output).map(|bytes| LoadedBlock {
                        address: job.address,
                        size: job.size,
                        bytes,
                    });
                    if result_sender
                        .send(DecompressResult {
                            sequence: job.sequence,
                            input,
                            block,
                        })
                        .is_err()
                    {
                        break;
                    }
                }
            }));
            senders.push(sender);
        }
        drop(result_sender);
        Self {
            file,
            file_len,
            senders,
            results,
            workers,
            next_worker: 0,
            next_read_address: 0,
            expected_address: 0,
            submitted: 0,
            delivered: 0,
            ready: Vec::new(),
            spare_inputs: Vec::new(),
            spare_outputs: Vec::new(),
            max_in_flight: threads.saturating_mul(2),
        }
    }

    pub fn load(&mut self, address: u64) -> io::Result<LoadedBlock> {
        if address != self.expected_address {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "Parallel BGZF reader received a non-sequential offset",
            ));
        }
        self.fill_queue()?;
        let result = self.receive_next()?;
        self.delivered += 1;
        self.spare_inputs.push(result.input);
        let block = result.block?;
        if block.address != address {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Parallel BGZF blocks were delivered out of order",
            ));
        }
        self.expected_address = self
            .expected_address
            .checked_add(block.size)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "BGZF offset overflow"))?;
        Ok(block)
    }

    pub fn recycle_output(&mut self, output: Vec<u8>) {
        self.spare_outputs.push(output);
    }

    fn fill_queue(&mut self) -> io::Result<()> {
        while self.submitted - self.delivered < self.max_in_flight
            && self.next_read_address < self.file_len
        {
            let address = self.next_read_address;
            let input = self.spare_inputs.pop().unwrap_or_default();
            let (compressed, size) =
                read_compressed_block(&mut self.file, self.file_len, address, input)?;
            let output = self.spare_outputs.pop().unwrap_or_default();
            self.senders[self.next_worker]
                .send(DecompressJob {
                    sequence: self.submitted,
                    address,
                    size,
                    compressed,
                    output,
                })
                .map_err(|_| io::Error::other("Parallel BGZF worker stopped unexpectedly"))?;
            self.next_worker = (self.next_worker + 1) % self.senders.len();
            self.submitted += 1;
            self.next_read_address = self.next_read_address.checked_add(size).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "BGZF offset overflow")
            })?;
        }
        Ok(())
    }

    fn receive_next(&mut self) -> io::Result<DecompressResult> {
        loop {
            if let Some(position) = self
                .ready
                .iter()
                .position(|result| result.sequence == self.delivered)
            {
                return Ok(self.ready.swap_remove(position));
            }
            let result = self
                .results
                .recv()
                .map_err(|_| io::Error::other("Parallel BGZF worker stopped unexpectedly"))?;
            if result.sequence == self.delivered {
                return Ok(result);
            }
            self.ready.push(result);
        }
    }
}

impl Drop for ParallelBlockLoader {
    fn drop(&mut self) {
        self.senders.clear();
        for worker in self.workers.drain(..) {
            let _ = worker.join();
        }
    }
}
