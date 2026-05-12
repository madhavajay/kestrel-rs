/// Converts a string to Java-style name case.
#[must_use]
pub fn to_name_case(value: impl AsRef<str>) -> String {
    let value = value.as_ref();

    if value.is_empty() {
        return String::new();
    }

    let mut chars = value.chars();
    let Some(first) = chars.next() else {
        return String::new();
    };

    first
        .to_uppercase()
        .chain(chars.flat_map(char::to_lowercase))
        .collect()
}

/// Returns a printable description of a Unicode code point.
#[must_use]
pub fn char_description(codepoint: u32) -> String {
    if (0xD800..=0xDBFF).contains(&codepoint) {
        return format!("<HIGH_SURROGATE> (U+{codepoint:04X}, HIGH SURROGATE)");
    }

    if (0xDC00..=0xDFFF).contains(&codepoint) {
        return format!("<LOW_SURROGATE> (U+{codepoint:04X}, LOW SURROGATE)");
    }

    let Some(ch) = char::from_u32(codepoint) else {
        return format!("<ILLEGAL_CODEPOINT> (0x{codepoint:X}, {codepoint})");
    };

    let print_char = if ch.is_whitespace() {
        "<WHITESPACE_CHAR>".to_owned()
    } else if ch.is_control() {
        "<ISO_CTRL_CHAR>".to_owned()
    } else {
        ch.to_string()
    };

    let name = control_name(codepoint)
        .map(str::to_owned)
        .or_else(|| unicode_names2::name(ch).map(|name| name.to_string()))
        .unwrap_or_else(|| "UNKNOWN CHARACTER NAME".to_owned());

    format!("{print_char} (U+{codepoint:04X}, {name})")
}

fn control_name(codepoint: u32) -> Option<&'static str> {
    let name = match codepoint {
        0x00 => "NULL",
        0x01 => "START OF HEADING",
        0x02 => "START OF TEXT",
        0x03 => "END OF TEXT",
        0x04 => "END OF TRANSMISSION",
        0x05 => "ENQUIRY",
        0x06 => "ACKNOWLEDGE",
        0x07 => "BELL",
        0x08 => "BACKSPACE",
        0x09 => "CHARACTER TABULATION",
        0x0A => "LINE FEED",
        0x0B => "LINE TABULATION",
        0x0C => "FORM FEED",
        0x0D => "CARRIAGE RETURN",
        0x0E => "SHIFT OUT",
        0x0F => "SHIFT IN",
        0x10 => "DATA LINK ESCAPE",
        0x11 => "DEVICE CONTROL ONE",
        0x12 => "DEVICE CONTROL TWO",
        0x13 => "DEVICE CONTROL THREE",
        0x14 => "DEVICE CONTROL FOUR",
        0x15 => "NEGATIVE ACKNOWLEDGE",
        0x16 => "SYNCHRONOUS IDLE",
        0x17 => "END OF TRANSMISSION BLOCK",
        0x18 => "CANCEL",
        0x19 => "END OF MEDIUM",
        0x1A => "SUBSTITUTE",
        0x1B => "ESCAPE",
        0x1C => "INFORMATION SEPARATOR FOUR",
        0x1D => "INFORMATION SEPARATOR THREE",
        0x1E => "INFORMATION SEPARATOR TWO",
        0x1F => "INFORMATION SEPARATOR ONE",
        0x7F => "DELETE",
        0x80 => "PADDING CHARACTER",
        0x81 => "HIGH OCTET PRESET",
        0x82 => "BREAK PERMITTED HERE",
        0x83 => "NO BREAK HERE",
        0x84 => "INDEX",
        0x85 => "NEXT LINE",
        0x86 => "START OF SELECTED AREA",
        0x87 => "END OF SELECTED AREA",
        0x88 => "CHARACTER TABULATION SET",
        0x89 => "CHARACTER TABULATION WITH JUSTIFICATION",
        0x8A => "LINE TABULATION SET",
        0x8B => "PARTIAL LINE FORWARD",
        0x8C => "PARTIAL LINE BACKWARD",
        0x8D => "REVERSE LINE FEED",
        0x8E => "SINGLE SHIFT TWO",
        0x8F => "SINGLE SHIFT THREE",
        0x90 => "DEVICE CONTROL STRING",
        0x91 => "PRIVATE USE ONE",
        0x92 => "PRIVATE USE TWO",
        0x93 => "SET TRANSMIT STATE",
        0x94 => "CANCEL CHARACTER",
        0x95 => "MESSAGE WAITING",
        0x96 => "START OF GUARDED AREA",
        0x97 => "END OF GUARDED AREA",
        0x98 => "START OF STRING",
        0x99 => "SINGLE GRAPHIC CHARACTER INTRODUCER",
        0x9A => "SINGLE CHARACTER INTRODUCER",
        0x9B => "CONTROL SEQUENCE INTRODUCER",
        0x9C => "STRING TERMINATOR",
        0x9D => "OPERATING SYSTEM COMMAND",
        0x9E => "PRIVACY MESSAGE",
        0x9F => "APPLICATION PROGRAM COMMAND",
        _ => return None,
    };

    Some(name)
}

/// Returns Unicode scalar values for a string.
#[must_use]
pub fn code_point_array(value: impl AsRef<str>) -> Vec<u32> {
    value.as_ref().chars().map(u32::from).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name_case_matches_java_examples() {
        let cases = [
            ("one", "One"),
            ("1one", "1one"),
            ("", ""),
            ("&", "&"),
            ("One!andtwo", "One!andtwo"),
            ("one!andtwo", "One!andtwo"),
            ("thisisit", "Thisisit"),
        ];

        for (input, expected) in cases {
            assert_eq!(to_name_case(input), expected);
        }
    }

    #[test]
    fn code_points_preserve_non_bmp_scalars() {
        assert_eq!(code_point_array("A\u{1F9EC}T"), vec![0x41, 0x1F9EC, 0x54]);
    }

    #[test]
    fn character_descriptions_classify_special_values() {
        assert_eq!(char_description(0x20), "<WHITESPACE_CHAR> (U+0020, SPACE)");
        assert_eq!(char_description(0x00), "<ISO_CTRL_CHAR> (U+0000, NULL)");
        assert_eq!(
            char_description(0x11_0000),
            "<ILLEGAL_CODEPOINT> (0x110000, 1114112)"
        );
    }
}
