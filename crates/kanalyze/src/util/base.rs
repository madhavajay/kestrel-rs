use std::fmt;

/// DNA base encoded with Java KAnalyze numeric values.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[repr(u8)]
pub enum Base {
    /// Adenine.
    A = 0,
    /// Cytosine.
    C = 1,
    /// Guanine.
    G = 2,
    /// Thymine.
    T = 3,
}

/// Error returned when a value cannot be parsed as a DNA base.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ParseBaseError;

impl Base {
    /// All supported DNA bases in Java value order.
    pub const ALL: [Self; 4] = [Self::A, Self::C, Self::G, Self::T];

    /// Returns the uppercase ASCII base character.
    #[must_use]
    pub const fn as_char(self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::G => 'G',
            Self::T => 'T',
        }
    }

    /// Returns the uppercase ASCII base byte.
    #[must_use]
    pub const fn as_byte(self) -> u8 {
        self.as_char() as u8
    }

    /// Returns the Java numeric base value.
    #[must_use]
    pub const fn value(self) -> u8 {
        self as u8
    }

    /// Returns the Watson-Crick complement.
    #[must_use]
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
        }
    }

    /// Converts a Java numeric base value into a base.
    #[must_use]
    pub const fn from_value(value: u8) -> Option<Self> {
        match value {
            0 => Some(Self::A),
            1 => Some(Self::C),
            2 => Some(Self::G),
            3 => Some(Self::T),
            _ => None,
        }
    }

    /// Converts an uppercase DNA base character into a base.
    #[must_use]
    pub const fn from_char(base: char) -> Option<Self> {
        match base {
            'A' => Some(Self::A),
            'C' => Some(Self::C),
            'G' => Some(Self::G),
            'T' => Some(Self::T),
            _ => None,
        }
    }
}

impl TryFrom<u8> for Base {
    type Error = ParseBaseError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        Self::from_value(value).ok_or(ParseBaseError)
    }
}

impl TryFrom<char> for Base {
    type Error = ParseBaseError;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        Self::from_char(value).ok_or(ParseBaseError)
    }
}

impl From<Base> for char {
    fn from(value: Base) -> Self {
        value.as_char()
    }
}

impl From<Base> for u8 {
    fn from(value: Base) -> Self {
        value.value()
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.as_char().to_string())
    }
}

impl fmt::Display for ParseBaseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid DNA base")
    }
}

impl std::error::Error for ParseBaseError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn java_values_are_preserved() {
        let cases = [
            (Base::A, 'A', 0, Base::T),
            (Base::C, 'C', 1, Base::G),
            (Base::G, 'G', 2, Base::C),
            (Base::T, 'T', 3, Base::A),
        ];

        for (base, ch, value, complement) in cases {
            assert_eq!(base.as_char(), ch);
            assert_eq!(base.as_byte(), ch as u8);
            assert_eq!(base.value(), value);
            assert_eq!(base.complement(), complement);
            assert_eq!(Base::from_value(value), Some(base));
            assert_eq!(Base::from_char(ch), Some(base));
        }
    }

    #[test]
    fn invalid_values_are_missing() {
        assert_eq!(Base::from_value(4), None);
        assert_eq!(Base::from_char('N'), None);
        assert_eq!(Base::try_from(4), Err(ParseBaseError));
        assert_eq!(Base::try_from('a'), Err(ParseBaseError));
    }
}
