use thiserror::Error;

/// Floating-point tolerance treated as zero.
pub const ZERO_RANGE: f32 = 0.0001;
/// Negative floating-point tolerance treated as zero.
pub const N_ZERO_RANGE: f32 = -ZERO_RANGE;

/// Errors returned while computing quantiles.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum QuantileError {
    /// Quantile was outside the supported range.
    #[error("q is out of range (0, 1): {0}")]
    OutOfRange(String),
}

/// Returns true when a value is within the zero tolerance.
#[must_use]
pub fn is_zero(value: f32) -> bool {
    (N_ZERO_RANGE..=ZERO_RANGE).contains(&value)
}

/// Returns true when a value is outside the zero tolerance.
#[must_use]
pub fn is_not_zero(value: f32) -> bool {
    !(N_ZERO_RANGE..=ZERO_RANGE).contains(&value)
}

/// Computes a quantile over adjacent absolute count differences.
pub fn count_diff_quantile(count: &[i32], quantile: f64) -> Result<i32, QuantileError> {
    if quantile <= 0.0 || quantile >= 1.0 {
        return Err(QuantileError::OutOfRange(quantile.to_string()));
    }

    match count.len() {
        0 | 1 => Ok(0),
        2 => Ok((count[1] - count[0]).abs()),
        len => {
            let mut diffs = count
                .windows(2)
                .map(|window| (window[0] - window[1]).abs())
                .collect::<Vec<_>>();
            diffs.sort_unstable();

            let n_count = len - 2;
            let raw = n_count as f64 * quantile;
            let loc = raw as usize;
            let offset = raw - loc as f64;
            let value = diffs[loc] as f64 * (1.0 - offset) + diffs[loc + 1] as f64 * offset;

            Ok(value as i32)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_range_constants_match_java() {
        assert_eq!(ZERO_RANGE, 0.0001);
        assert_eq!(N_ZERO_RANGE, -0.0001);
    }

    #[test]
    fn zero_checks_match_java_boundaries() {
        assert!(is_zero(0.0));
        assert!(is_zero(0.00001));
        assert!(is_zero(-0.00001));
        assert!(is_zero(ZERO_RANGE));
        assert!(is_zero(N_ZERO_RANGE));
        assert!(!is_zero(0.001));
        assert!(!is_zero(-0.001));

        assert!(!is_not_zero(ZERO_RANGE));
        assert!(!is_not_zero(N_ZERO_RANGE));
        assert!(is_not_zero(0.001));
        assert!(is_not_zero(-0.001));
    }

    #[test]
    fn count_diff_quantile_matches_java_examples() {
        assert_eq!(count_diff_quantile(&[], 0.5).unwrap(), 0);
        assert_eq!(count_diff_quantile(&[42], 0.5).unwrap(), 0);
        assert_eq!(count_diff_quantile(&[10, 5], 0.5).unwrap(), 5);
        assert_eq!(count_diff_quantile(&[5, 10], 0.5).unwrap(), 5);
        assert_eq!(count_diff_quantile(&[10, 5, 20, 7, 30], 0.5).unwrap(), 14);
    }

    #[test]
    fn count_diff_quantile_rejects_out_of_range_q() {
        assert_eq!(
            count_diff_quantile(&[1, 2, 3], 0.0),
            Err(QuantileError::OutOfRange("0".to_owned()))
        );
        assert!(matches!(
            count_diff_quantile(&[1, 2, 3], 1.0),
            Err(QuantileError::OutOfRange(_))
        ));
    }
}
