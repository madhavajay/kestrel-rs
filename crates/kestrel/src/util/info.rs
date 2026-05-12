use std::str::FromStr;

use crate::constants::{VERSION, VERSION_DEV, VERSION_MAJOR, VERSION_MINOR};

/// Version information fields.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InfoMode {
    /// Major version.
    VersionMajor,
    /// Minor version.
    VersionMinor,
    /// Development version.
    VersionDev,
    /// Full version string.
    Version,
}

impl InfoMode {
    /// All supported info modes.
    pub const ALL: [Self; 4] = [
        Self::VersionMajor,
        Self::VersionMinor,
        Self::VersionDev,
        Self::Version,
    ];

    /// Returns the dotted mode name.
    #[must_use]
    pub fn mode_name(self) -> &'static str {
        match self {
            Self::VersionMajor => "version.major",
            Self::VersionMinor => "version.minor",
            Self::VersionDev => "version.dev",
            Self::Version => "version",
        }
    }

    /// Returns the value for this info mode.
    #[must_use]
    pub fn value(self) -> String {
        match self {
            Self::VersionMajor => VERSION_MAJOR.to_string(),
            Self::VersionMinor => VERSION_MINOR.to_string(),
            Self::VersionDev => VERSION_DEV.to_string(),
            Self::Version => VERSION.to_owned(),
        }
    }
}

impl FromStr for InfoMode {
    type Err = ();

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let value = value.trim().to_ascii_lowercase();

        Self::ALL
            .into_iter()
            .find(|mode| mode.mode_name() == value)
            .ok_or(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mode_names_are_lowercase_dotted() {
        assert_eq!(InfoMode::VersionMajor.mode_name(), "version.major");
        assert_eq!(InfoMode::VersionMinor.mode_name(), "version.minor");
        assert_eq!(InfoMode::VersionDev.mode_name(), "version.dev");
        assert_eq!(InfoMode::Version.mode_name(), "version");
    }

    #[test]
    fn modes_parse_case_insensitively() {
        assert_eq!("version".parse::<InfoMode>(), Ok(InfoMode::Version));
        assert_eq!(
            "Version.Major".parse::<InfoMode>(),
            Ok(InfoMode::VersionMajor)
        );
        assert!("garbage".parse::<InfoMode>().is_err());
    }

    #[test]
    fn mode_values_match_constants() {
        assert_eq!(InfoMode::VersionMajor.value(), "1");
        assert_eq!(InfoMode::VersionMinor.value(), "0");
        assert_eq!(InfoMode::VersionDev.value(), "0");
        assert_eq!(InfoMode::Version.value(), "1.0.2");
    }
}
