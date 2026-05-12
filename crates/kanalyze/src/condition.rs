use std::fmt;
use std::sync::{Arc, Mutex};

/// Severity of a reported condition.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ConditionLevel {
    /// Non-fatal warning condition.
    Warning,
    /// Error condition.
    Error,
}

/// Condition notification emitted by a component.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConditionEvent {
    /// Condition severity.
    pub level: ConditionLevel,
    /// Human-readable condition message.
    pub message: String,
}

impl ConditionEvent {
    /// Creates a warning event.
    #[must_use]
    pub fn warning(message: impl Into<String>) -> Self {
        Self {
            level: ConditionLevel::Warning,
            message: message.into(),
        }
    }

    /// Creates an error event.
    #[must_use]
    pub fn error(message: impl Into<String>) -> Self {
        Self {
            level: ConditionLevel::Error,
            message: message.into(),
        }
    }
}

/// Listener for component condition events.
pub trait ConditionListener: Send + Sync {
    /// Handles a condition event.
    fn on_condition(&self, event: &ConditionEvent);
}

impl<F> ConditionListener for F
where
    F: Fn(&ConditionEvent) + Send + Sync,
{
    fn on_condition(&self, event: &ConditionEvent) {
        self(event);
    }
}

/// In-memory condition listener used by tests and simple integrations.
#[derive(Clone, Default)]
pub struct StreamConditionListener {
    events: Arc<Mutex<Vec<ConditionEvent>>>,
}

impl StreamConditionListener {
    /// Creates an empty stream listener.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns all collected events.
    #[must_use]
    pub fn events(&self) -> Vec<ConditionEvent> {
        self.events
            .lock()
            .expect("condition mutex poisoned")
            .clone()
    }
}

impl ConditionListener for StreamConditionListener {
    fn on_condition(&self, event: &ConditionEvent) {
        self.events
            .lock()
            .expect("condition mutex poisoned")
            .push(event.clone());
    }
}

impl fmt::Debug for StreamConditionListener {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("StreamConditionListener")
            .field("events", &self.events())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn listener_collects_events() {
        let listener = StreamConditionListener::new();

        listener.on_condition(&ConditionEvent::warning("low coverage"));
        listener.on_condition(&ConditionEvent::error("bad input"));

        assert_eq!(
            listener.events(),
            vec![
                ConditionEvent::warning("low coverage"),
                ConditionEvent::error("bad input")
            ]
        );
    }
}
