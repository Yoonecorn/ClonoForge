//! Bounded ordered parallel map. Results are committed in input order.
//!
//! A permit is returned only after the corresponding result is consumed, so
//! queued jobs, running jobs and out-of-order results together never exceed
//! `window`. `window` is a record count, not a byte budget: callers must also
//! enforce an admitted bound on every job and result before enqueueing.
use crate::model::invalid;
use std::{
    collections::BTreeMap,
    io,
    panic::{catch_unwind, AssertUnwindSafe},
    sync::{
        atomic::{AtomicBool, Ordering},
        mpsc, Arc, Mutex,
    },
    thread,
};

enum Event<O> {
    Item(u64, io::Result<O>),
    End(u64),
    WorkerExit(Option<u64>),
}

/// Apply independent worker operations and consume the results in input order.
/// A terminal error stops admission, drops channel endpoints, and joins threads.
/// It cannot interrupt an operating-system read or worker closure already blocked
/// inside external code; the caller must supply finite/cancellable operations.
pub fn ordered_map<I, T, O, F, C>(
    input: I,
    workers: usize,
    window: usize,
    work: F,
    consume: C,
) -> io::Result<()>
where
    I: Iterator<Item = io::Result<T>> + Send,
    T: Send,
    O: Send,
    F: Fn(T) -> io::Result<O> + Sync,
    C: FnMut(O) -> io::Result<()>,
{
    ordered_map_inner(input, workers, window, work, consume, 0)
}
pub(crate) fn ordered_map_inner<I, T, O, F, C>(
    input: I,
    workers: usize,
    window: usize,
    work: F,
    mut consume: C,
    fault: u8,
) -> io::Result<()>
where
    I: Iterator<Item = io::Result<T>> + Send,
    T: Send,
    O: Send,
    F: Fn(T) -> io::Result<O> + Sync,
    C: FnMut(O) -> io::Result<()>,
{
    if workers == 0 || window == 0 {
        return Err(invalid("executor: zero workers or window"));
    }
    let (task_tx, task_rx) = mpsc::sync_channel::<(u64, T)>(window);
    let task_rx = Arc::new(Mutex::new(task_rx));
    let (event_tx, event_rx) = mpsc::sync_channel::<Event<O>>(window);
    let (permit_tx, permit_rx) = mpsc::sync_channel::<()>(window);
    for _ in 0..window {
        permit_tx
            .send(())
            .map_err(|_| invalid("executor: permit initialization"))?;
    }
    let cancelled = Arc::new(AtomicBool::new(false));
    thread::scope(|scope| {
        let producer_cancel = Arc::clone(&cancelled);
        let producer_events = event_tx.clone();
        let producer = thread::Builder::new().spawn_scoped(scope, move || {
            let mut count = 0u64;
            let mut input = input;
            loop {
                if producer_cancel.load(Ordering::Acquire) || permit_rx.recv().is_err() {
                    break;
                }
                if producer_cancel.load(Ordering::Acquire) {
                    break;
                }
                let item = catch_unwind(AssertUnwindSafe(|| input.next()));
                match item {
                    Ok(Some(Ok(job))) => {
                        let Some(next) = count.checked_add(1) else {
                            producer_cancel.store(true, Ordering::Release);
                            let _ = producer_events.send(Event::Item(
                                count,
                                Err(invalid("executor: ordinal overflow")),
                            ));
                            return;
                        };
                        if task_tx.send((count, job)).is_err() {
                            break;
                        }
                        count = next;
                    }
                    Ok(None) => {
                        let _ = producer_events.send(Event::End(count));
                        return;
                    }
                    Ok(Some(Err(e))) => {
                        producer_cancel.store(true, Ordering::Release);
                        let _ = producer_events.send(Event::Item(count, Err(e)));
                        return;
                    }
                    Err(_) => {
                        producer_cancel.store(true, Ordering::Release);
                        let _ = producer_events
                            .send(Event::Item(count, Err(invalid("executor: input panicked"))));
                        return;
                    }
                }
            }
            let _ = producer_events.send(Event::End(count));
        })?;
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            let rx = Arc::clone(&task_rx);
            let events = event_tx.clone();
            let worker_cancel = Arc::clone(&cancelled);
            let work = &work;
            let handle = thread::Builder::new().spawn_scoped(scope, move || {
                let mut active = None;
                let _ = catch_unwind(AssertUnwindSafe(|| loop {
                    let item = match rx.lock() {
                        Ok(rx) => rx.recv(),
                        Err(_) => break,
                    };
                    let Ok((ordinal, job)) = item else {
                        break;
                    };
                    active = Some(ordinal);
                    if fault == 1 || (fault == 2 && ordinal == 1) {
                        break;
                    }
                    let result = match catch_unwind(AssertUnwindSafe(|| work(job))) {
                        Ok(r) => r,
                        Err(_) => Err(invalid("executor: worker panicked")),
                    };
                    if result.is_err() {
                        worker_cancel.store(true, Ordering::Release);
                    }
                    if events.send(Event::Item(ordinal, result)).is_err() {
                        break;
                    }
                    active = None;
                }));
                let _ = events.send(Event::WorkerExit(active));
            });
            match handle {
                Ok(handle) => handles.push(handle),
                Err(error) => {
                    cancelled.store(true, Ordering::Release);
                    drop(event_rx);
                    drop(permit_tx);
                    drop(task_rx);
                    drop(event_tx);
                    let _ = producer.join();
                    for h in handles {
                        let _ = h.join();
                    }
                    return Err(error);
                }
            }
        }
        drop(event_tx);
        drop(task_rx);
        let mut next = 0u64;
        let mut pending = BTreeMap::new();
        let mut end = None;
        let mut exited = 0usize;
        let result = loop {
            if end == Some(next) {
                break Ok(());
            }
            let event = match event_rx.recv() {
                Ok(e) => e,
                Err(_) => break Err(invalid("executor: result channel closed early")),
            };
            match event {
                Event::End(total) => {
                    if end.replace(total).is_some() || total < next {
                        break Err(invalid("executor: invalid producer completion"));
                    }
                }
                Event::Item(index, result) => {
                    if index < next || pending.insert(index, result).is_some() {
                        break Err(invalid("executor: duplicate result ordinal"));
                    }
                }
                Event::WorkerExit(abandoned) => {
                    exited += 1;
                    if let Some(index) = abandoned {
                        if index < next
                            || pending
                                .insert(index, Err(invalid("executor: worker exited unexpectedly")))
                                .is_some()
                        {
                            break Err(invalid("executor: invalid worker exit ordinal"));
                        }
                    }
                }
            }
            while let Some(result) = pending.remove(&next) {
                let result = result.and_then(|value| {
                    catch_unwind(AssertUnwindSafe(|| consume(value)))
                        .unwrap_or_else(|_| Err(invalid("executor: consumer panicked")))
                });
                match result {
                    Ok(()) => {}
                    Err(e) => {
                        cancelled.store(true, Ordering::Release);
                        drop(event_rx);
                        drop(permit_tx);
                        let _ = producer.join();
                        for h in handles {
                            let _ = h.join();
                        }
                        return Err(e);
                    }
                }
                next = match next.checked_add(1) {
                    Some(n) => n,
                    None => break,
                };
                // A closed permit channel means production has ended; no refill is needed.
                let _ = permit_tx.send(());
            }
            if exited == workers && end != Some(next) {
                break Err(invalid(
                    "executor: all workers exited before results completed",
                ));
            }
        };
        cancelled.store(true, Ordering::Release);
        drop(event_rx);
        drop(permit_tx);
        let p = producer.join();
        let mut panicked = p.is_err();
        for h in handles {
            panicked |= h.join().is_err();
        }
        if panicked && result.is_ok() {
            Err(invalid("executor: thread panicked"))
        } else {
            result
        }
    })
}

#[cfg(test)]
mod closure_cancellation_tests {
    use super::*;
    use std::{sync::atomic::AtomicUsize, time::Duration};
    #[test]
    fn all_workers_exit_early_closes_queue_without_deadlock() {
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let result = ordered_map_inner((0..100).map(Ok), 4, 4, Ok, |_: i32| Ok(()), 1);
            tx.send(result.map_err(|e| e.to_string())).unwrap();
        });
        let result = rx
            .recv_timeout(Duration::from_secs(2))
            .expect("worker-exit queue closure deadlocked");
        assert!(result.is_err());
    }
    #[test]
    fn one_worker_exits_with_a_job_without_deadlock() {
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let result = ordered_map_inner((0..100).map(Ok), 4, 4, Ok, |_: i32| Ok(()), 2);
            tx.send(result.map_err(|e| e.to_string())).unwrap();
        });
        let result = rx
            .recv_timeout(Duration::from_secs(2))
            .expect("single-worker exit deadlocked");
        assert!(result.unwrap_err().contains("worker exited"));
    }
    #[test]
    fn consumer_cancellation_joins_all_active_workers() {
        struct Active<'a>(&'a AtomicUsize);
        impl Drop for Active<'_> {
            fn drop(&mut self) {
                self.0.fetch_sub(1, Ordering::SeqCst);
            }
        }
        let active = AtomicUsize::new(0);
        let result = ordered_map(
            (0..100).map(Ok),
            4,
            8,
            |n| {
                active.fetch_add(1, Ordering::SeqCst);
                let _active = Active(&active);
                std::thread::sleep(Duration::from_millis(5));
                Ok(n)
            },
            |_| Err(invalid("consumer cancellation")),
        );
        assert!(result.is_err());
        assert_eq!(active.load(Ordering::SeqCst), 0);
    }
}
