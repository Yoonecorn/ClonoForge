use clonoforge_core::ordered_executor::ordered_map;
use std::{
    io,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
    thread,
    time::Duration,
};
fn err() -> io::Error {
    io::Error::other("sentinel")
}
#[test]
fn ordered_bounded_backpressure() {
    for workers in [1, 2, 8] {
        for _ in 0..4 {
            let live = Arc::new(AtomicUsize::new(0));
            let peak = Arc::new(AtomicUsize::new(0));
            let input_live = live.clone();
            let input_peak = peak.clone();
            let input = (0..300).map(move |n| {
                let x = input_live.fetch_add(1, Ordering::SeqCst) + 1;
                input_peak.fetch_max(x, Ordering::SeqCst);
                Ok(n)
            });
            let mut out = Vec::new();
            ordered_map(
                input,
                workers,
                7,
                |n| {
                    if n % 19 == 0 {
                        thread::sleep(Duration::from_millis(1));
                    }
                    Ok(n * 2)
                },
                |n| {
                    out.push(n);
                    live.fetch_sub(1, Ordering::SeqCst);
                    Ok(())
                },
            )
            .unwrap();
            assert_eq!(out, (0..300).map(|x| x * 2).collect::<Vec<_>>());
            assert!(peak.load(Ordering::SeqCst) <= 7);
        }
    }
}
#[test]
fn errors_and_panics_terminate() {
    for workers in [1, 2, 8] {
        for kind in 0..5 {
            let input = (0..200).map(move |n| {
                if kind == 0 && n == 21 {
                    Err(err())
                } else {
                    Ok(n)
                }
            });
            let result = ordered_map(
                input,
                workers,
                4,
                move |n| {
                    if kind == 1 && n == 20 {
                        return Err(err());
                    }
                    if kind == 2 && n == 20 {
                        panic!("worker");
                    }
                    Ok(n)
                },
                move |n| {
                    if kind == 3 && n == 20 {
                        return Err(err());
                    }
                    if kind == 4 && n == 20 {
                        panic!("consumer");
                    }
                    Ok(())
                },
            );
            assert!(result.is_err(), "{workers} {kind}");
        }
    }
}
#[test]
fn invalid_configuration() {
    assert!(ordered_map(std::iter::empty::<io::Result<u32>>(), 0, 1, Ok, |_| Ok(())).is_err());
    assert!(ordered_map(std::iter::empty::<io::Result<u32>>(), 1, 0, Ok, |_| Ok(())).is_err());
}

#[test]
fn earlier_submitted_failure_wins_over_faster_later_failure() {
    for _ in 0..20 {
        let result = ordered_map(
            (0..20).map(Ok),
            4,
            8,
            |n| {
                if n == 1 {
                    thread::sleep(Duration::from_millis(5));
                    return Err(io::Error::other("earlier"));
                }
                if n == 3 {
                    return Err(io::Error::other("later"));
                }
                Ok(n)
            },
            |_| Ok(()),
        );
        assert_eq!(result.unwrap_err().to_string(), "earlier");
    }
}
