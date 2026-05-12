//! Rust port of the Kestrel variant caller.

/// Active-region detection and haplotype data structures.
pub mod activeregion;
/// Alignment scoring, trace, and k-mer alignment types.
pub mod align;
/// Command-line parsing for the Rust Kestrel binary.
pub mod cli;
/// Kestrel constants mirrored from the Java implementation.
pub mod constants;
/// Count map implementations used by the runner.
pub mod counter;
/// Haplotype writer registry and writer implementations.
pub mod hapwriter;
/// Genomic interval parsing and storage.
pub mod interval;
/// Shared input/output value types.
pub mod io;
/// Logging level parsing and tracing integration.
pub mod log_level;
/// Reference sequence reader and region types.
pub mod refreader;
/// Kestrel runner configuration and orchestration.
pub mod runner;
/// Kestrel utility modules.
pub mod util;
/// Variant filter registry and filter implementations.
pub mod varfilter;
/// Variant call types and caller.
pub mod variant;
/// Variant writer registry and writer implementations.
pub mod writer;

pub use constants::*;
pub use log_level::{LogLevel, LogLevelError};

#[cfg(test)]
mod tests {
    #[test]
    fn crate_loads() {
        assert_eq!(env!("CARGO_PKG_NAME"), "kestrel");
    }
}
