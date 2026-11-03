use clonoforge_core::admission::Request;
fn request() -> Request {
    Request {
        budget: u64::MAX,
        workers: 2,
        window: 4,
        merge_fan_in: 2,
        max_record_bytes: 4096,
        max_group_records: 50,
        sort_chunk_bytes: 1024 * 1024,
    }
}
#[test]
fn deterministic_exact_budget_boundary() {
    let mut r = request();
    let p = r.plan().unwrap();
    assert_eq!(p, r.plan().unwrap());
    r.budget = p.total;
    assert_eq!(r.plan().unwrap(), p);
    r.budget = p.total - 1;
    assert!(r.plan().is_err());
}
#[test]
fn rejects_invalid_and_overflow() {
    let mut r = request();
    r.workers = 0;
    assert!(r.plan().is_err());
    r = request();
    r.merge_fan_in = 1;
    assert!(r.plan().is_err());
    r = request();
    r.max_group_records = 51;
    assert!(r.plan().is_err());
    r = request();
    r.window = u64::MAX;
    assert!(r.plan().is_err());
    r = request();
    r.sort_chunk_bytes = u64::MAX;
    assert!(r.plan().is_err());
}
#[test]
fn increasing_resources_never_reduces_model() {
    let r = request();
    let a = r.plan().unwrap();
    let mut b = r;
    b.workers += 1;
    assert!(b.plan().unwrap().total > a.total);
    b = r;
    b.window += 1;
    assert!(b.plan().unwrap().total > a.total);
    b = r;
    b.merge_fan_in += 1;
    assert!(b.plan().unwrap().total > a.total);
    b = r;
    b.max_record_bytes += 1;
    assert!(b.plan().unwrap().total > a.total);
}
#[test]
fn deterministic_resource_selection() {
    let (r, p) = Request::select(u64::MAX, 8, 4096, 50).unwrap();
    assert_eq!(r.workers, 8);
    assert_eq!(r.window, 8);
    assert_eq!(p.max_line_bytes, 4099);
    let (r2, p2) = Request::select(p.total, 8, 4096, 50).unwrap();
    assert_eq!(p, p2);
    assert_eq!(r.workers, r2.workers);
    let (r3, _) = Request::select(p2.total - 1, 8, 4096, 50).unwrap();
    assert_eq!(r3.workers, 8);
    assert_eq!(r3.merge_fan_in, 2);
    let mut minimum_request = r3;
    minimum_request.workers = 1;
    minimum_request.window = 1;
    minimum_request.budget = u64::MAX;
    let minimum = minimum_request.plan().unwrap();
    assert!(Request::select(minimum.total - 1, 8, 4096, 50).is_err());
    let (boundary, plan) = Request::select(minimum.total, 1, 4096, 50).unwrap();
    assert_eq!(boundary.workers, 1);
    assert_eq!(plan.total, minimum.total);
}
