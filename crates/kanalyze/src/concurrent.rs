use crossbeam_channel::{Receiver, RecvError, SendError, Sender, bounded};

/// Cloneable bounded queue backed by a crossbeam channel.
#[derive(Debug)]
pub struct BoundedQueue<T> {
    sender: Sender<T>,
    receiver: Receiver<T>,
}

impl<T> BoundedQueue<T> {
    /// Creates a queue with the given channel capacity.
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        let (sender, receiver) = bounded(capacity);
        Self { sender, receiver }
    }

    /// Sends a value to the queue.
    pub fn push(&self, value: T) -> Result<(), SendError<T>> {
        self.sender.send(value)
    }

    /// Receives the next value from the queue.
    pub fn pop(&self) -> Result<T, RecvError> {
        self.receiver.recv()
    }

    /// Returns a cloned sender handle.
    #[must_use]
    pub fn sender(&self) -> Sender<T> {
        self.sender.clone()
    }

    /// Returns a cloned receiver handle.
    #[must_use]
    pub fn receiver(&self) -> Receiver<T> {
        self.receiver.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sends_and_receives_values() {
        let queue = BoundedQueue::new(1);

        queue.push(7).unwrap();

        assert_eq!(queue.pop().unwrap(), 7);
    }

    #[test]
    fn cloned_sender_reports_shutdown() {
        let queue = BoundedQueue::new(1);
        let sender = queue.sender();
        drop(queue);

        assert_eq!(sender.send(1).unwrap_err().0, 1);
    }
}
