//! Verbose logging utilities.
//!
//! When enabled (the default), log messages are printed to **stderr** and
//! buffered in a shared ring buffer so the GUI can display them in a
//! scrollable log panel.
//!
//! Toggle with [`set_verbose`].  The flag is a global [`AtomicBool`] so it
//! is safe to read from any thread.

use std::collections::VecDeque;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

/// Maximum number of log lines retained in the GUI buffer.
const MAX_LOG_LINES: usize = 500;

/// Global verbose flag — **on** by default.
static VERBOSE: AtomicBool = AtomicBool::new(true);

/// Returns `true` when verbose logging is enabled.
pub fn is_verbose() -> bool {
    VERBOSE.load(Ordering::Relaxed)
}

/// Enable or disable verbose logging.
pub fn set_verbose(on: bool) {
    VERBOSE.store(on, Ordering::Relaxed);
}

/// Thread-safe ring buffer of log messages for GUI display.
#[derive(Clone, Debug)]
pub struct LogBuffer {
    inner: Arc<Mutex<VecDeque<String>>>,
}

impl Default for LogBuffer {
    fn default() -> Self {
        Self::new()
    }
}

impl LogBuffer {
    pub fn new() -> Self {
        Self {
            inner: Arc::new(Mutex::new(VecDeque::with_capacity(MAX_LOG_LINES))),
        }
    }

    /// Append a message to the buffer (and print to stderr when verbose).
    pub fn log(&self, msg: &str) {
        if is_verbose() {
            eprintln!("[viewer] {msg}");
        }
        if let Ok(mut buf) = self.inner.lock() {
            if buf.len() >= MAX_LOG_LINES {
                buf.pop_front();
            }
            buf.push_back(msg.to_string());
        }
    }

    /// Return a snapshot of all buffered lines.
    pub fn lines(&self) -> Vec<String> {
        self.inner
            .lock()
            .map(|buf| buf.iter().cloned().collect())
            .unwrap_or_default()
    }

    /// Number of buffered lines.
    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.inner.lock().map(|buf| buf.len()).unwrap_or(0)
    }
}

/// Convenience: log a message through a [`LogBuffer`] reference.
///
/// Usage: `vlog!(buf, "loaded {} reads", n);`
macro_rules! vlog {
    ($buf:expr, $($arg:tt)*) => {
        $buf.log(&format!($($arg)*))
    };
}
pub(crate) use vlog;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_verbose_default_is_on() {
        assert!(is_verbose());
    }

    #[test]
    fn test_set_verbose_off_and_on() {
        // Temporarily set to off.  Restore at end so other tests aren't affected.
        set_verbose(false);
        assert!(!is_verbose());
        set_verbose(true);
        assert!(is_verbose());
    }

    #[test]
    fn test_log_buffer_stores_messages() {
        let buf = LogBuffer::new();
        buf.log("hello");
        buf.log("world");
        let lines = buf.lines();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "hello");
        assert_eq!(lines[1], "world");
    }

    #[test]
    fn test_log_buffer_ring_eviction() {
        let buf = LogBuffer::new();
        for i in 0..MAX_LOG_LINES + 10 {
            buf.log(&format!("line {i}"));
        }
        assert_eq!(buf.len(), MAX_LOG_LINES);
        // First line should be "line 10" (first 10 were evicted).
        let lines = buf.lines();
        assert_eq!(lines[0], "line 10");
    }

    #[test]
    fn test_vlog_macro() {
        let buf = LogBuffer::new();
        vlog!(buf, "count = {}", 42);
        assert_eq!(buf.lines(), vec!["count = 42"]);
    }

    #[test]
    fn test_log_buffer_clone_shares_state() {
        let buf = LogBuffer::new();
        let clone = buf.clone();
        buf.log("from original");
        assert_eq!(clone.lines(), vec!["from original"]);
    }
}
