//! Processor time for the timing tests.
//!
//! The machine the tests run on is shared with other work, and wall time
//! stretches with however busy it is. On Linux [`cpu_time`] runs its
//! closure with the kernel's parallel stages on the calling thread and
//! reads that thread's processor clock (`CLOCK_THREAD_CPUTIME_ID`) around
//! it: the reading is the work the closure did, and time the thread spent
//! waiting for a core is not in it. Elsewhere it reads wall time and
//! leaves the thread count alone.

use std::time::Duration;

/// Run `f`, returning its result and the processor time it took on this
/// thread (wall time off Linux).
pub fn cpu_time<R>(f: impl FnOnce() -> R) -> (R, Duration) {
    #[cfg(target_os = "linux")]
    {
        let _one = OneThread::hold();
        let start = thread_clock();
        let out = f();
        let took = thread_clock().saturating_sub(start);
        (out, took)
    }
    #[cfg(not(target_os = "linux"))]
    {
        let start = ogeom_core::clock::Instant::now();
        let out = f();
        (out, start.elapsed())
    }
}

/// The processor time this thread has run.
#[cfg(target_os = "linux")]
fn thread_clock() -> Duration {
    use rustix::time::{ClockId, clock_gettime};
    Duration::try_from(clock_gettime(ClockId::ThreadCPUTime)).unwrap_or_default()
}

/// The kernel's parallel stages run on their caller while any measurement
/// in this test binary holds one of these. Tests run side by side in one
/// process and the thread count is process-wide, so the last holder to
/// let go restores the machine default. The answers do not depend on the
/// thread count; only the time the other tests take does.
#[cfg(target_os = "linux")]
struct OneThread;

#[cfg(target_os = "linux")]
static HOLDERS: std::sync::Mutex<usize> = std::sync::Mutex::new(0);

#[cfg(target_os = "linux")]
impl OneThread {
    fn hold() -> Self {
        let mut holders = HOLDERS
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if *holders == 0 {
            ogeom_core::parallel::set_threads(1);
        }
        *holders += 1;
        Self
    }
}

#[cfg(target_os = "linux")]
impl Drop for OneThread {
    fn drop(&mut self) {
        let mut holders = HOLDERS
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        *holders -= 1;
        if *holders == 0 {
            ogeom_core::parallel::set_threads(0);
        }
    }
}
