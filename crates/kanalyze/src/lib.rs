//! Rust port of the KAnalyze subset used by Kestrel.

/// Sequence batch storage helpers.
pub mod batch;
/// Sequence source and reader components.
pub mod comp;
/// Small concurrency primitives used by pipeline code.
pub mod concurrent;
/// Condition reporting support.
pub mod condition;
/// Binary KAnalyze file formats.
pub mod io;
/// Runnable KAnalyze modules.
pub mod module;
/// Shared KAnalyze utility types.
pub mod util;

pub use util::Base;

/// Maximum Java-compatible array size.
pub const MAX_ARRAY_SIZE: usize = i32::MAX as usize - 8;

/// Common lifecycle contract for KAnalyze-style runnable components.
pub trait KAnalyzeRunnable {
    /// Error type returned by the runnable.
    type Error;

    /// Runs the component to completion.
    fn run(&mut self) -> Result<(), Self::Error>;
    /// Requests that the component stop as soon as practical.
    fn abort(&mut self);
    /// Returns whether abort has been requested.
    fn is_aborted(&self) -> bool;
}

#[cfg(test)]
mod tests {
    #[test]
    fn crate_loads() {
        assert_eq!(env!("CARGO_PKG_NAME"), "kanalyze");
    }
}
