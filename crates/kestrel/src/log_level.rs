use std::fmt;
use std::str::FromStr;

use thiserror::Error;
use tracing::Level;

/// Kestrel log levels accepted by the command-line interface.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum LogLevel {
    /// All logging, mapped to trace.
    All,
    /// Trace logging.
    Trace,
    /// Debug logging.
    Debug,
    /// Informational logging.
    Info,
    /// Warning logging.
    Warn,
    /// Error logging.
    Error,
    /// Disable logging.
    Off,
}

/// Errors returned while parsing log levels.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum LogLevelError {
    /// Log level string was empty.
    #[error("log level name is empty")]
    Empty,
    /// Log level string was not recognized.
    #[error("unknown log level name: {0}")]
    Unknown(String),
}

impl LogLevel {
    /// All supported log levels in display order.
    pub const ALL: [Self; 7] = [
        Self::All,
        Self::Trace,
        Self::Debug,
        Self::Info,
        Self::Warn,
        Self::Error,
        Self::Off,
    ];

    /// Converts to a `tracing` level, or `None` for disabled logging.
    #[must_use]
    pub fn tracing_level(self) -> Option<Level> {
        match self {
            Self::All | Self::Trace => Some(Level::TRACE),
            Self::Debug => Some(Level::DEBUG),
            Self::Info => Some(Level::INFO),
            Self::Warn => Some(Level::WARN),
            Self::Error => Some(Level::ERROR),
            Self::Off => None,
        }
    }

    /// Returns a comma-separated list of log level names.
    #[must_use]
    pub fn level_list() -> String {
        Self::level_list_with_and(false)
    }

    /// Returns a comma-separated list of log level names, optionally adding `and`.
    #[must_use]
    pub fn level_list_with_and(with_and: bool) -> String {
        let names = Self::ALL.map(|level| level.to_string());

        if with_and {
            let (last, rest) = names.split_last().expect("LogLevel::ALL is nonempty");
            format!("{}, and {last}", rest.join(", "))
        } else {
            names.join(", ")
        }
    }
}

impl FromStr for LogLevel {
    type Err = LogLevelError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let trimmed = value.trim();

        if trimmed.is_empty() {
            return Err(LogLevelError::Empty);
        }

        match trimmed.to_ascii_uppercase().as_str() {
            "ALL" => Ok(Self::All),
            "TRACE" => Ok(Self::Trace),
            "DEBUG" => Ok(Self::Debug),
            "INFO" => Ok(Self::Info),
            "WARN" => Ok(Self::Warn),
            "ERROR" => Ok(Self::Error),
            "OFF" => Ok(Self::Off),
            _ => Err(LogLevelError::Unknown(trimmed.to_owned())),
        }
    }
}

impl fmt::Display for LogLevel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::All => "ALL",
            Self::Trace => "TRACE",
            Self::Debug => "DEBUG",
            Self::Info => "INFO",
            Self::Warn => "WARN",
            Self::Error => "ERROR",
            Self::Off => "OFF",
        };

        f.write_str(name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levels_map_to_tracing() {
        assert_eq!(LogLevel::All.tracing_level(), Some(Level::TRACE));
        assert_eq!(LogLevel::Trace.tracing_level(), Some(Level::TRACE));
        assert_eq!(LogLevel::Debug.tracing_level(), Some(Level::DEBUG));
        assert_eq!(LogLevel::Info.tracing_level(), Some(Level::INFO));
        assert_eq!(LogLevel::Warn.tracing_level(), Some(Level::WARN));
        assert_eq!(LogLevel::Error.tracing_level(), Some(Level::ERROR));
        assert_eq!(LogLevel::Off.tracing_level(), None);
    }

    #[test]
    fn parses_level_names_like_java() {
        assert_eq!("INFO".parse::<LogLevel>().unwrap(), LogLevel::Info);
        assert_eq!("Debug".parse::<LogLevel>().unwrap(), LogLevel::Debug);
        assert_eq!("  warn  ".parse::<LogLevel>().unwrap(), LogLevel::Warn);
    }

    #[test]
    fn rejects_empty_or_unknown_names() {
        assert_eq!("".parse::<LogLevel>(), Err(LogLevelError::Empty));
        assert_eq!(
            "MAYBE".parse::<LogLevel>(),
            Err(LogLevelError::Unknown("MAYBE".to_owned()))
        );
    }

    #[test]
    fn level_list_mentions_all_levels() {
        let list = LogLevel::level_list();

        for level in LogLevel::ALL {
            assert!(list.contains(&level.to_string()));
        }
        assert!(!list.contains(", and "));
        assert!(LogLevel::level_list_with_and(true).contains(", and "));
    }
}
