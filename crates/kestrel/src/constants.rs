/// Program name used in command-line help and diagnostics.
pub const PROG_NAME: &str = "kestrel";

/// Major version component matching the Java Kestrel release.
pub const VERSION_MAJOR: u32 = 1;
/// Minor version component matching the Java Kestrel release.
pub const VERSION_MINOR: u32 = 0;
/// Revision version component matching the Java Kestrel release.
pub const VERSION_REVISION: u32 = 2;
/// Development version component.
pub const VERSION_DEV: u32 = 0;
/// Full semantic version string.
pub const VERSION: &str = "1.0.2";

/// Maximum array-like allocation size mirrored from Java.
pub const MAX_ARRAY_SIZE: usize = i32::MAX as usize - 8;
/// Minimum k-mer size supported by Kestrel.
pub const MIN_KMER_SIZE: usize = 4;
/// Default array expansion factor.
pub const ARRAY_EXPAND_FACTOR: f32 = 1.5;

/// No-error process status.
pub const ERR_NONE: i32 = 0;
/// Usage-error process status.
pub const ERR_USAGE: i32 = 1;
/// I/O-error process status.
pub const ERR_IO: i32 = 2;
/// Security-error process status.
pub const ERR_SECURITY: i32 = 3;
/// File-not-found process status.
pub const ERR_FILE_NOT_FOUND: i32 = 4;
/// Data-format-error process status.
pub const ERR_DATA_FORMAT: i32 = 5;
/// Analysis-error process status.
pub const ERR_ANALYSIS: i32 = 6;
/// Interrupted process status.
pub const ERR_INTERRUPTED: i32 = 7;
/// Resource-limit process status.
pub const ERR_LIMITS: i32 = 8;
/// Abort process status.
pub const ERR_ABORT: i32 = 98;
/// System-error process status.
pub const ERR_SYSTEM: i32 = 99;

/// Java resource root path retained for compatibility.
pub const RESOURCE_ROOT: &str = "edu/gatech/kestrel";
/// Java test resource path retained for compatibility.
pub const RESOURCE_TEST: &str = "edu/gatech/kestrel/test";
/// Valid format-type identifier pattern.
pub const FORMAT_TYPE_PATTERN: &str = "[A-Za-z][A-Za-z0-9_-]*";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn program_name_is_kestrel() {
        assert_eq!(PROG_NAME, "kestrel");
    }

    #[test]
    fn version_parts_match_version() {
        assert_eq!(VERSION_MAJOR, 1);
        assert_eq!(VERSION_MINOR, 0);
        assert_eq!(VERSION_REVISION, 2);
        assert_eq!(VERSION_DEV, 0);
        assert_eq!(VERSION, "1.0.2");
    }

    #[test]
    fn limits_match_java() {
        assert_eq!(MIN_KMER_SIZE, 4);
        assert_eq!(ARRAY_EXPAND_FACTOR, 1.5);
        assert_eq!(MAX_ARRAY_SIZE, i32::MAX as usize - 8);
    }

    #[test]
    fn error_codes_are_distinct() {
        let codes = [
            ERR_NONE,
            ERR_USAGE,
            ERR_IO,
            ERR_SECURITY,
            ERR_FILE_NOT_FOUND,
            ERR_DATA_FORMAT,
            ERR_ANALYSIS,
            ERR_INTERRUPTED,
            ERR_LIMITS,
            ERR_ABORT,
            ERR_SYSTEM,
        ];

        for (left_index, left) in codes.iter().enumerate() {
            for right in codes.iter().skip(left_index + 1) {
                assert_ne!(left, right);
            }
        }
        assert_eq!(ERR_NONE, 0);
    }

    #[test]
    fn resource_paths_match_java() {
        assert_eq!(RESOURCE_ROOT, "edu/gatech/kestrel");
        assert_eq!(RESOURCE_TEST, "edu/gatech/kestrel/test");
    }
}
