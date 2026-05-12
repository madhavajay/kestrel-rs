use std::any::type_name_of_val;
use std::path::PathBuf;

use thiserror::Error;

/// Errors from resource lookup helpers.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum ResourceError {
    /// Resource name was empty after trimming.
    #[error("resource name is empty")]
    EmptyName,
    /// Resource could not be found in known search roots.
    #[error("resource file does not exist: {0}")]
    NotFound(String),
}

/// Resolves a resource name to an existing filesystem path.
pub fn get_file_by_resource(resource_name: impl AsRef<str>) -> Result<PathBuf, ResourceError> {
    let resource_name = resource_name.as_ref().trim();

    if resource_name.is_empty() {
        return Err(ResourceError::EmptyName);
    }

    let candidates = [
        PathBuf::from(resource_name),
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(resource_name),
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join(resource_name),
    ];

    candidates
        .into_iter()
        .find(|path| path.exists())
        .map(|path| path.canonicalize().unwrap_or(path))
        .ok_or_else(|| ResourceError::NotFound(resource_name.to_owned()))
}

/// Formats an object identity string similar to Java's `Object.toString`.
#[must_use]
pub fn object_to_string<T: ?Sized>(value: Option<&T>) -> String {
    match value {
        Some(value) => format!(
            "{}@{:x}",
            type_name_of_val(value),
            value as *const T as *const () as usize
        ),
        None => "null".to_owned(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_workspace_resources() {
        let path = get_file_by_resource("fixtures/refreader/general.us-ascii.fasta").unwrap();

        assert!(path.ends_with("fixtures/refreader/general.us-ascii.fasta"));
    }

    #[test]
    fn rejects_missing_or_empty_resource_names() {
        assert_eq!(get_file_by_resource(" "), Err(ResourceError::EmptyName));
        assert_eq!(
            get_file_by_resource("fixtures/refreader/missing.fasta"),
            Err(ResourceError::NotFound(
                "fixtures/refreader/missing.fasta".to_owned()
            ))
        );
    }

    #[test]
    fn formats_object_identity() {
        let value = 7_u32;
        let rendered = object_to_string(Some(&value));

        assert!(rendered.starts_with("u32@"));
        assert_eq!(object_to_string::<u32>(None), "null");
    }
}
