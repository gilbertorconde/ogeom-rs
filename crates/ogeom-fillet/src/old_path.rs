use std::cell::Cell;
use std::io::Write as _;
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Clone, Copy)]
struct State {
    request: u64,
    api: &'static str,
    public_depth: usize,
    entry_depth: usize,
}

thread_local! {
    static STATE: Cell<State> = const { Cell::new(State {
        request: 0, api: "-", public_depth: 0, entry_depth: 0,
    }) };
}
static NEXT_REQUEST: AtomicU64 = AtomicU64::new(1);

#[must_use]
pub(crate) struct PublicRequest;

pub(crate) fn public_request(api: &'static str) -> PublicRequest {
    STATE.with(|state| {
        let mut current = state.get();
        if current.public_depth == 0 {
            current.request = NEXT_REQUEST.fetch_add(1, Ordering::Relaxed);
            current.api = api;
        }
        current.public_depth += 1;
        state.set(current);
    });
    PublicRequest
}

impl Drop for PublicRequest {
    fn drop(&mut self) {
        STATE.with(|state| {
            let mut current = state.get();
            current.public_depth -= 1;
            if current.public_depth == 0 {
                current.request = 0;
                current.api = "-";
            }
            state.set(current);
        });
    }
}

#[must_use]
pub(crate) struct Entry;

impl Drop for Entry {
    fn drop(&mut self) {
        STATE.with(|state| {
            let mut current = state.get();
            current.entry_depth -= 1;
            state.set(current);
        });
    }
}

pub(crate) fn entry(entry: &str) -> Entry {
    let current = STATE.with(|state| {
        let current = state.get();
        state.set(State {
            entry_depth: current.entry_depth + 1,
            ..current
        });
        current
    });
    let guard = Entry;
    if let Some(path) = std::env::var_os("OGEOM_OLD_PATH_LOG") {
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .unwrap_or_else(|error| panic!("cannot open old-path log: {error}"));
        let test = std::thread::current().name().unwrap_or("?").to_owned();
        let executable = std::env::current_exe()
            .unwrap_or_else(|error| panic!("cannot identify old-path caller: {error}"));
        // One append keeps concurrent test names on separate lines.
        let line = format!(
            "{entry}\t{test}\t{}\t{}\t{}\t{}\n",
            executable.display(),
            current.request,
            current.api,
            current.entry_depth,
        );
        file.write_all(line.as_bytes())
            .unwrap_or_else(|error| panic!("cannot write old-path log: {error}"));
    }
    guard
}
