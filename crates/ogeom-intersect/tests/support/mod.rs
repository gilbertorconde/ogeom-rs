//! The intersector's instruments.
//!
//! These are not kernel capability and they are not tests either. They are the
//! things a test asserts *about*: [`benchmark`] scores accuracy — every point
//! of every reported curve against both surfaces, with ground truth being the
//! surfaces themselves — and [`coverage`] scores completeness, asking the
//! surfaces by signed distance where the intersection must be and checking
//! that some branch reached there.
//!
//! They live under `tests/` because that is what they are for. Shipping them
//! as public API would advertise them as something a caller of the kernel
//! wants, and no caller does; what a caller wants is an intersector that has
//! been held to them. Both have negative controls in the suites that use them:
//! they demonstrably fail when something is genuinely missing.
//!
//! Each test binary that includes this module uses some of it and not the
//! rest, hence the blanket `dead_code` allowance — the alternative is a
//! per-item annotation that says nothing.
#![allow(dead_code, reason = "each test binary uses a different part")]

pub mod benchmark;
pub mod coverage;
