//! Diagnostics: the kernel's debug output, by area, off unless asked for.
//!
//! Each area of the kernel writes to its own [`Channel`]. Every channel is off
//! by default, and a probe on an off channel costs one atomic load. Channels
//! are switched on in one of two ways:
//!
//! - the `OGEOM_DEBUG` environment variable, read once on the first probe:
//!   channel names separated by commas or spaces (`OGEOM_DEBUG=sew,wire`), or
//!   `all`;
//! - [`set_enabled`], from the host or a test, which overrides the variable
//!   for that channel.
//!
//! Lines go to a sink. The default one writes each line to stderr, which on
//! `wasm32-unknown-unknown` goes nowhere. A host that wants the lines
//! elsewhere (a log, a console panel) installs its own with [`set_sink`].
//!
//! Probes write with the [`diag!`](crate::diag!) macro, which checks the
//! channel before it formats anything:
//!
//! ```
//! use ogeom_core::diag::{self, Channel};
//!
//! let gap = 1.5e-4;
//! ogeom_core::diag!(Channel::Sew, "gap {gap:.2e} left open");
//! if diag::enabled(Channel::Sew) {
//!     // Work done only to describe a failure goes behind the same check.
//! }
//! ```

use std::fmt;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Arc, PoisonError, RwLock};

/// The environment variable that names the channels to switch on.
pub const ENV: &str = "OGEOM_DEBUG";

/// An area of the kernel with debug output of its own.
///
/// Each has a stable name, the one `OGEOM_DEBUG` takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Channel {
    /// `arrange`: the boolean's per-face arrangement, its pieces and their
    /// states, and the open edges of a result shell that does not close.
    Arrange,
    /// `bool-audit`: fill every face pair of a boolean a second time without
    /// the bounding-box filter and assert both fills keep the same section
    /// material. A check rather than output: it doubles the fill's cost and
    /// panics on a disagreement, so `all` leaves it off.
    BoolAudit,
    /// `defeature`: the wire a face removal rebuilds, and the bridge it could
    /// not finish.
    Defeature,
    /// `domain`: a curve evaluated outside its parameter domain.
    Domain,
    /// `edge`: an edge whose end vertex lies off its curve.
    Edge,
    /// `fillet`: where a rolling-ball blend runs out, and why.
    Fillet,
    /// `mass`: how the mass properties walk faces and orient shells.
    Mass,
    /// `mesh`: faces the tessellator drew again finer, refused constraint
    /// segments, and the slow faces' timings.
    Mesh,
    /// `sew`: pcurve refits and the wires and faces the sewing could not
    /// build.
    Sew,
    /// `strands`: the boolean's strands per face, snapping and dust.
    Strands,
    /// `sweep`: a sweep or ring whose sewn shell does not close.
    Sweep,
    /// `wire`: the boolean's pave, junction and wire-building decisions, and
    /// gaps between consecutive edges of a wire.
    Wire,
}

impl Channel {
    /// Every channel, in name order.
    pub const ALL: &'static [Self] = &[
        Self::Arrange,
        Self::BoolAudit,
        Self::Defeature,
        Self::Domain,
        Self::Edge,
        Self::Fillet,
        Self::Mass,
        Self::Mesh,
        Self::Sew,
        Self::Strands,
        Self::Sweep,
        Self::Wire,
    ];

    /// The channel's name, as `OGEOM_DEBUG` takes it.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Arrange => "arrange",
            Self::BoolAudit => "bool-audit",
            Self::Defeature => "defeature",
            Self::Domain => "domain",
            Self::Edge => "edge",
            Self::Fillet => "fillet",
            Self::Mass => "mass",
            Self::Mesh => "mesh",
            Self::Sew => "sew",
            Self::Strands => "strands",
            Self::Sweep => "sweep",
            Self::Wire => "wire",
        }
    }

    /// The channel with this name, ignoring ASCII case.
    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL
            .iter()
            .copied()
            .find(|c| c.name().eq_ignore_ascii_case(name))
    }

    const fn bit(self) -> u32 {
        1 << self as u32
    }

    /// Whether `all` switches the channel on: every channel that only
    /// writes.
    const fn in_all(self) -> bool {
        !matches!(self, Self::BoolAudit)
    }
}

impl fmt::Display for Channel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// Set while the environment has not been read yet. No channel uses it.
const UNREAD: u32 = 1 << 31;

/// The enabled channels, one bit each.
static MASK: AtomicU32 = AtomicU32::new(UNREAD);

/// Whether `channel` is on.
#[inline]
#[must_use]
pub fn enabled(channel: Channel) -> bool {
    let mut mask = MASK.load(Ordering::Relaxed);
    if mask & UNREAD != 0 {
        mask = read_environment();
    }
    mask & channel.bit() != 0
}

/// Switch `channel` on or off for the whole process. Overrides `OGEOM_DEBUG`
/// for that channel and leaves the others as they are.
pub fn set_enabled(channel: Channel, on: bool) {
    if MASK.load(Ordering::Acquire) & UNREAD != 0 {
        read_environment();
    }
    if on {
        MASK.fetch_or(channel.bit(), Ordering::AcqRel);
    } else {
        MASK.fetch_and(!channel.bit(), Ordering::AcqRel);
    }
}

/// Read `OGEOM_DEBUG` into the mask, once: the thread that stores it wins,
/// and any other sees its answer.
#[cold]
fn read_environment() -> u32 {
    let spec = std::env::var(ENV).unwrap_or_default();
    let (mask, unknown) = parse(&spec);
    match MASK.compare_exchange(UNREAD, mask, Ordering::AcqRel, Ordering::Acquire) {
        Ok(_) => {
            if !unknown.is_empty() {
                let known: Vec<&str> = Channel::ALL.iter().map(|c| c.name()).collect();
                write_stderr(format_args!(
                    "{ENV}: no channel named {}; the channels are all, {}",
                    unknown.join(", "),
                    known.join(", ")
                ));
            }
            mask
        }
        Err(current) => current,
    }
}

/// The mask a channel list asks for, and the names it does not know.
fn parse(spec: &str) -> (u32, Vec<&str>) {
    let mut mask = 0;
    let mut unknown = Vec::new();
    for name in spec
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
    {
        if name.eq_ignore_ascii_case("all") {
            for c in Channel::ALL.iter().filter(|c| c.in_all()) {
                mask |= c.bit();
            }
        } else if let Some(c) = Channel::from_name(name) {
            mask |= c.bit();
        } else {
            unknown.push(name);
        }
    }
    (mask, unknown)
}

/// A sink: hears each line with the channel it was written to.
type Sink = Arc<dyn Fn(Channel, fmt::Arguments<'_>) + Send + Sync>;

/// The host's sink, when it installed one.
static SINK: RwLock<Option<Sink>> = RwLock::new(None);

/// Send every diagnostic line to `sink` from now on, in place of stderr.
///
/// One sink serves the whole process. It runs on whichever thread writes the
/// line, possibly several at once, and only for channels that are on.
pub fn set_sink(sink: impl Fn(Channel, fmt::Arguments<'_>) + Send + Sync + 'static) {
    *SINK.write().unwrap_or_else(PoisonError::into_inner) = Some(Arc::new(sink));
}

/// Send diagnostic lines to stderr again, as before any [`set_sink`].
pub fn reset_sink() {
    *SINK.write().unwrap_or_else(PoisonError::into_inner) = None;
}

/// Write one line to `channel`'s sink, whether or not the channel is on.
/// Probes go through [`diag!`](crate::diag!), which checks first.
pub fn emit(channel: Channel, line: fmt::Arguments<'_>) {
    let sink = SINK.read().unwrap_or_else(PoisonError::into_inner).clone();
    match sink {
        Some(sink) => sink(channel, line),
        None => write_stderr(line),
    }
}

/// The default sink: the line on stderr.
#[allow(clippy::print_stderr)]
fn write_stderr(line: fmt::Arguments<'_>) {
    eprintln!("{line}");
}

/// Write a formatted line to a diagnostics channel, formatting nothing when
/// the channel is off.
///
/// ```
/// use ogeom_core::diag::Channel;
///
/// let edges = 3;
/// ogeom_core::diag!(Channel::Wire, "{edges} edges left over");
/// ```
#[macro_export]
macro_rules! diag {
    ($channel:expr, $($arg:tt)+) => {{
        let channel: $crate::diag::Channel = $channel;
        if $crate::diag::enabled(channel) {
            $crate::diag::emit(channel, format_args!($($arg)+));
        }
    }};
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    #[test]
    fn every_channel_has_its_own_bit_below_the_unread_flag() {
        let mut seen = 0;
        for c in Channel::ALL {
            assert_eq!(seen & c.bit(), 0, "{c} shares a bit");
            assert_eq!(c.bit() & UNREAD, 0, "{c} uses the unread flag");
            seen |= c.bit();
        }
    }

    #[test]
    fn names_round_trip() {
        for &c in Channel::ALL {
            assert_eq!(Channel::from_name(c.name()), Some(c));
            assert_eq!(Channel::from_name(&c.name().to_uppercase()), Some(c));
        }
        assert_eq!(Channel::from_name("nothing"), None);
    }

    #[test]
    fn a_list_parses_by_comma_or_space_and_reports_what_it_does_not_know() {
        let (mask, unknown) = parse("sew, wire mass,,typo");
        assert_eq!(
            mask,
            Channel::Sew.bit() | Channel::Wire.bit() | Channel::Mass.bit()
        );
        assert_eq!(unknown, ["typo"]);
        assert_eq!(parse(""), (0, vec![]));
    }

    #[test]
    fn all_leaves_the_audit_to_be_named() {
        let (mask, _) = parse("all");
        assert_eq!(mask & Channel::BoolAudit.bit(), 0);
        assert_ne!(mask & Channel::Arrange.bit(), 0);
        let (mask, _) = parse("all,bool-audit");
        assert_ne!(mask & Channel::BoolAudit.bit(), 0);
    }

    // The one test that touches the process-wide switches and sink, so no
    // other test sees them change.
    #[test]
    fn a_probe_reaches_the_sink_only_while_its_channel_is_on() {
        let heard: Arc<Mutex<Vec<String>>> = Arc::default();
        let into = Arc::clone(&heard);
        set_sink(move |c, line| into.lock().unwrap().push(format!("{c}: {line}")));

        set_enabled(Channel::Sweep, false);
        let mut formatted = 0;
        let mut count = || {
            formatted += 1;
            formatted
        };
        crate::diag!(Channel::Sweep, "off {}", count());
        set_enabled(Channel::Sweep, true);
        crate::diag!(Channel::Sweep, "on {}", count());
        set_enabled(Channel::Sweep, false);
        reset_sink();

        assert_eq!(*heard.lock().unwrap(), ["sweep: on 1"]);
    }
}
