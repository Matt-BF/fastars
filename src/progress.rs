use std::io::{self, IsTerminal, Write};
use std::time::{Duration, Instant};

const REFRESH_INTERVAL: Duration = Duration::from_secs(1);
const LOG_INTERVAL: Duration = Duration::from_secs(30);
const LOG_PERCENT_STEP: u64 = 10;

#[derive(Clone, Copy)]
pub(crate) enum Unit {
    Bytes,
    Records,
    Runs,
}

pub(crate) struct ProgressBar {
    enabled: bool,
    terminal: bool,
    phase: String,
    total: u64,
    unit: Unit,
    started: Instant,
    last_draw: Instant,
    last_log: Instant,
    last_sample: Instant,
    last_value: u64,
    rate: Option<f64>,
    next_log_percent: u64,
    rendered_width: usize,
    line_active: bool,
}

impl ProgressBar {
    pub(crate) fn new(enabled: bool, phase: impl Into<String>, total: u64, unit: Unit) -> Self {
        let now = Instant::now();
        let mut progress = Self {
            enabled,
            terminal: io::stderr().is_terminal(),
            phase: phase.into(),
            total,
            unit,
            started: now,
            last_draw: now,
            last_log: now,
            last_sample: now,
            last_value: 0,
            rate: None,
            next_log_percent: LOG_PERCENT_STEP,
            rendered_width: 0,
            line_active: false,
        };
        if progress.enabled {
            progress.draw(0, false);
        }
        progress
    }

    pub(crate) fn update(&mut self, value: u64) {
        if !self.enabled || value >= self.total {
            return;
        }
        let now = Instant::now();
        let percent = percentage(value, self.total) as u64;
        let due = if self.terminal {
            now.duration_since(self.last_draw) >= REFRESH_INTERVAL
        } else {
            percent >= self.next_log_percent || now.duration_since(self.last_log) >= LOG_INTERVAL
        };
        if due {
            self.sample_rate(value, now);
            self.draw(value, false);
        }
    }

    pub(crate) fn finish(&mut self, value: u64) {
        if !self.enabled {
            return;
        }
        let now = Instant::now();
        self.sample_rate(value, now);
        self.draw(value, true);
        self.line_active = false;
    }

    fn sample_rate(&mut self, value: u64, now: Instant) {
        let elapsed = now.duration_since(self.last_sample).as_secs_f64();
        if elapsed > 0.0 && value >= self.last_value {
            let current = (value - self.last_value) as f64 / elapsed;
            self.rate = Some(match self.rate {
                Some(previous) => previous * 0.7 + current * 0.3,
                None => current,
            });
        }
        self.last_sample = now;
        self.last_value = value;
    }

    fn draw(&mut self, value: u64, finished: bool) {
        let value = value.min(self.total);
        let elapsed = self.started.elapsed();
        let mut line = format!(
            "[{}] {} / {} ({:.1}%)",
            self.phase,
            format_value(value, self.unit),
            format_value(self.total, self.unit),
            percentage(value, self.total)
        );
        if finished {
            line.push_str(&format!(" — done in {}", format_duration(elapsed)));
        } else if let Some(rate) = self.rate.filter(|rate| *rate > 0.0) {
            line.push_str(&format!(" — {}/s", format_rate(rate, self.unit)));
            let remaining = self.total.saturating_sub(value) as f64 / rate;
            if remaining.is_finite() {
                line.push_str(&format!(
                    " — ETA {}",
                    format_duration(Duration::from_secs_f64(remaining))
                ));
            }
        }

        if self.terminal {
            let padding = self.rendered_width.saturating_sub(line.len());
            eprint!("\r{line}{}", " ".repeat(padding));
            if finished {
                eprintln!();
            } else {
                let _ = io::stderr().flush();
                self.line_active = true;
            }
            self.rendered_width = line.len();
            self.last_draw = Instant::now();
        } else {
            eprintln!("{line}");
            self.last_log = Instant::now();
            let percent = percentage(value, self.total) as u64;
            self.next_log_percent = ((percent / LOG_PERCENT_STEP) + 1) * LOG_PERCENT_STEP;
        }
    }
}

impl Drop for ProgressBar {
    fn drop(&mut self) {
        if self.enabled && self.terminal && self.line_active {
            eprintln!();
        }
    }
}

pub(crate) fn status(enabled: bool, phase: &str, message: &str) {
    if enabled {
        eprintln!("[{phase}] {message}");
    }
}

pub(crate) fn format_duration(duration: Duration) -> String {
    let seconds = duration.as_secs();
    let hours = seconds / 3600;
    let minutes = seconds % 3600 / 60;
    let seconds = seconds % 60;
    if hours > 0 {
        format!("{hours}h {minutes}m {seconds}s")
    } else if minutes > 0 {
        format!("{minutes}m {seconds}s")
    } else {
        format!("{seconds}s")
    }
}

fn percentage(value: u64, total: u64) -> f64 {
    if total == 0 {
        100.0
    } else {
        value as f64 / total as f64 * 100.0
    }
}

fn format_value(value: u64, unit: Unit) -> String {
    match unit {
        Unit::Bytes => format_bytes(value as f64),
        Unit::Records => format!("{} records", separated(value)),
        Unit::Runs => format!("{} runs", separated(value)),
    }
}

fn format_rate(rate: f64, unit: Unit) -> String {
    match unit {
        Unit::Bytes => format_bytes(rate),
        Unit::Records => format!("{} records", separated(rate.round() as u64)),
        Unit::Runs => format!("{rate:.1} runs"),
    }
}

fn format_bytes(bytes: f64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes;
    let mut unit = 0;
    while value >= 1024.0 && unit < UNITS.len() - 1 {
        value /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{value:.0} {}", UNITS[unit])
    } else {
        format!("{value:.1} {}", UNITS[unit])
    }
}

fn separated(value: u64) -> String {
    let digits = value.to_string();
    let mut output = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, byte) in digits.bytes().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            output.push(',');
        }
        output.push(byte as char);
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_byte_units() {
        assert_eq!(format_bytes(512.0), "512 B");
        assert_eq!(format_bytes(1536.0), "1.5 KiB");
        assert_eq!(format_bytes(2.0 * 1024.0 * 1024.0), "2.0 MiB");
    }

    #[test]
    fn formats_counts_with_separators() {
        assert_eq!(separated(42), "42");
        assert_eq!(separated(1_234_567), "1,234,567");
    }

    #[test]
    fn formats_durations() {
        assert_eq!(format_duration(Duration::from_secs(9)), "9s");
        assert_eq!(format_duration(Duration::from_secs(271)), "4m 31s");
        assert_eq!(format_duration(Duration::from_secs(7_325)), "2h 2m 5s");
    }
}
