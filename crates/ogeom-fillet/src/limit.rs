//! The largest size a blend builds at, read from its seat and confirmed
//! against the build itself (ADR 0001, amendment "largest-size readings
//! confirmed by the build").
//!
//! The seat reading (where a contact first leaves its face) is fast and
//! closed-form, but a request can be refused below it for a reason the
//! seat does not see (a corner that does not close, a bevel reaching a
//! concave floor, a result failing its check), and can build above it
//! (a bevel rolling over its face's far side as a notch). So the seat
//! reading only starts the search: the size is then confirmed by asking
//! the same call the request makes, on a copy of the model, at sizes
//! bracketing the boundary between the last size that built and the first
//! that was refused, bisected to [`LIMIT_TOLERANCE`].

use crate::{Chamfer, ChamferCorner, FilletCorner, RoundProfile};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_topo::{Model, Shape};

/// The relative tolerance a confirmed limit is bisected to: the smallest
/// size found refused is at most this share above the size returned.
pub const LIMIT_TOLERANCE: f64 = 1e-4;

/// The most builds one reading probes. The bracket and the bisection to
/// [`LIMIT_TOLERANCE`] take at most about 45 from any start in a range of
/// a billion; the seat reading usually brackets the limit in two, and
/// [`LIMIT_SAMPLES`] each side are spent looking for a size that builds
/// again.
pub const LIMIT_PROBES: u32 = 64;

/// How many sizes below the limit, and how many above its first refusal,
/// a reading probes for a range that is not one interval: a refusal below
/// a size that builds, or a size that builds above a refusal.
pub const LIMIT_SAMPLES: u32 = 3;

/// What a largest-size reading found: the size and how it was confirmed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BlendLimit {
    /// The largest size confirmed to build, at most the bound asked: the
    /// size of the first boundary from below between sizes that build and
    /// sizes that are refused. None where no size probed above the
    /// confusion distance builds.
    pub size: Option<f64>,
    /// The smallest size probed above `size` that was refused, within
    /// [`LIMIT_TOLERANCE`] of it unless [`LIMIT_PROBES`] ran out first. None
    /// where the bound asked builds.
    pub refused: Option<f64>,
    /// A size above `refused` found to build again: the sizes that build
    /// are not one interval, and `size` is where they are first broken.
    /// None where none of the sizes probed for one builds, which does not
    /// prove there is none.
    pub builds_again: Option<f64>,
    /// The seat reading the search started from, where the seat could be
    /// read: the size where a contact first leaves its face.
    pub seat: Option<f64>,
    /// How many builds were probed.
    pub probes: u32,
}

/// A request whose largest size [`blend_limit`] reads, with the options
/// its blend call takes. The size scales what the request names.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum BlendSizing<'a> {
    /// [`fillet_edges_corner`](crate::fillet_edges_corner) on `edges`
    /// closed as `corner` says; the size is the radius.
    Rounds {
        /// The edges rounded together.
        edges: &'a [Shape],
        /// How a vertex where rounds meet is closed.
        corner: FilletCorner,
    },
    /// [`fillet_edge_distances`](crate::fillet_edge_distances) on `edge`,
    /// meeting `face` `size · on_face` from the edge and the other face
    /// `size · on_other` from it, its section `profile`; the size is the
    /// common scale.
    RoundDistances {
        /// The edge rounded.
        edge: &'a Shape,
        /// The face the first distance runs along.
        face: &'a Shape,
        /// The distance along `face` at a scale of one.
        on_face: f64,
        /// The distance along the other face at a scale of one.
        on_other: f64,
        /// The section between the contacts.
        profile: RoundProfile,
    },
    /// [`chamfer_edges_cornered`](crate::chamfer_edges_cornered) on `edge`
    /// alone, closed as `corner` says; the size is the measure
    /// [`largest_chamfer_distance`](crate::largest_chamfer_distance) reads
    /// for `spec`'s spelling, the measure `spec` itself names not read.
    /// One edge has no vertex where three of its bevels meet, so either
    /// corner builds alike there; the corner is the call's own option.
    Chamfer {
        /// The edge bevelled.
        edge: &'a Shape,
        /// The bevel's spelling.
        spec: &'a Chamfer,
        /// How a vertex where three or more bevels meet is closed.
        corner: ChamferCorner,
    },
}

/// The largest size up to `up_to` at which `sizing`'s blend call builds on
/// `solid`, with how it was confirmed: the reading behind
/// [`largest_fillet_radius_all`](crate::largest_fillet_radius_all),
/// [`largest_fillet_distances`](crate::largest_fillet_distances) and
/// [`largest_chamfer_distance`](crate::largest_chamfer_distance), which
/// return its [`BlendLimit::size`], with the options their calls take.
///
/// The seat is read first, without building, as those readings did; then
/// the call is asked on a copy of the model, first at the seat reading
/// (at `up_to` where the seat cannot be read) and a [`LIMIT_TOLERANCE`]
/// past it, then doubling or halving the size until one builds and one
/// is refused, and bisecting between the largest size that built and the
/// smallest refused until they stand within [`LIMIT_TOLERANCE`] of each
/// other. [`LIMIT_SAMPLES`] sizes evenly below the result, and as many
/// spread geometrically from the refusal to `up_to`, are then probed: a
/// refusal below moves the answer down to the boundary under it, a build
/// above is reported as [`BlendLimit::builds_again`]. At most
/// [`LIMIT_PROBES`] builds are probed.
///
/// A reading where the seat brackets the limit costs `2 + 2 ·`
/// [`LIMIT_SAMPLES`] builds of the request; one bisected from a factor of
/// two costs about 14 more.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `up_to`, an edge list or a distance is not usable; where the seat
/// itself is refused and no size builds, that refusal.
pub fn blend_limit(
    model: &Model,
    solid: &Shape,
    sizing: &BlendSizing<'_>,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<BlendLimit> {
    let conf = tol.confusion();
    if !up_to.is_finite() || up_to <= conf {
        ogeom_bail!(
            Construction,
            "a largest size up to {up_to} asks for nothing"
        );
    }
    let (seat, floor) = match *sizing {
        BlendSizing::Rounds { edges, .. } => {
            if edges.is_empty() {
                ogeom_bail!(Construction, "a request of no edges has no largest radius");
            }
            (
                crate::fits::seat_fillet_radius(model, solid, edges, up_to, tol),
                conf,
            )
        }
        BlendSizing::RoundDistances {
            edge,
            face,
            on_face,
            on_other,
            ..
        } => {
            for d in [on_face, on_other] {
                if !d.is_finite() || d <= conf {
                    ogeom_bail!(
                        Construction,
                        "a round set back {d} from its edge rounds nothing"
                    );
                }
            }
            (
                crate::asym::seat_fillet_distances(
                    model, solid, edge, face, on_face, on_other, up_to, tol,
                ),
                conf / on_face.min(on_other),
            )
        }
        BlendSizing::Chamfer { edge, spec, .. } => {
            let floor = match spec {
                Chamfer::Distances {
                    on_face, on_other, ..
                } => {
                    for d in [*on_face, *on_other] {
                        if !d.is_finite() || d <= conf {
                            ogeom_bail!(
                                Construction,
                                "a bevel set back {d} from its edge bevels nothing"
                            );
                        }
                    }
                    conf / on_face.min(*on_other)
                }
                Chamfer::Symmetric(_) | Chamfer::Angle { .. } => conf,
            };
            (
                crate::chamfer::seat_chamfer_distance(model, solid, edge, spec, up_to, tol),
                floor,
            )
        }
    };
    let builds = |size: f64| probe(model, solid, sizing, size, tol);
    match seat {
        Ok(seat) => Ok(confirm(Some(seat), up_to, floor, builds)),
        Err(refusal) => {
            let found = confirm(None, up_to, floor, builds);
            if found.size.is_none() {
                return Err(refusal);
            }
            Ok(found)
        }
    }
}

/// Whether `sizing`'s blend call builds at `size`, asked on a copy of the
/// model as the request asks it.
fn probe(
    model: &Model,
    solid: &Shape,
    sizing: &BlendSizing<'_>,
    size: f64,
    tol: Tolerances,
) -> bool {
    let mut scratch = model.clone();
    match *sizing {
        BlendSizing::Rounds { edges, corner } => {
            crate::fillet_edges_corner(&mut scratch, solid, edges, size, corner, tol).is_ok()
        }
        BlendSizing::RoundDistances {
            edge,
            face,
            on_face,
            on_other,
            profile,
        } => crate::fillet_edge_distances(
            &mut scratch,
            solid,
            edge,
            face,
            size * on_face,
            size * on_other,
            profile,
            tol,
        )
        .is_ok(),
        BlendSizing::Chamfer { edge, spec, corner } => {
            let specs = [(edge.clone(), scaled(spec, size))];
            crate::chamfer_edges_cornered(&mut scratch, solid, &specs, corner, tol).is_ok()
        }
    }
}

/// `spec` at the measure `size`, as
/// [`largest_chamfer_distance`](crate::largest_chamfer_distance) reads it:
/// the distance, the scale of both distances, or the distance at the
/// angle.
pub(crate) fn scaled(spec: &Chamfer, size: f64) -> Chamfer {
    match spec {
        Chamfer::Symmetric(_) => Chamfer::Symmetric(size),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => Chamfer::Distances {
            face: face.clone(),
            on_face: size * on_face,
            on_other: size * on_other,
        },
        Chamfer::Angle { face, angle, .. } => Chamfer::Angle {
            face: face.clone(),
            distance: size,
            angle: *angle,
        },
    }
}

/// The probes one reading has asked, each size once.
struct Probes<F> {
    builds: F,
    asked: Vec<(f64, bool)>,
}

impl<F: FnMut(f64) -> bool> Probes<F> {
    /// Whether the call builds at `size`; None once [`LIMIT_PROBES`] are
    /// spent.
    fn ask(&mut self, size: f64) -> Option<bool> {
        if let Some(&(_, built)) = self.asked.iter().find(|(s, _)| *s == size) {
            return Some(built);
        }
        if self.spent() >= LIMIT_PROBES {
            return None;
        }
        let built = (self.builds)(size);
        self.asked.push((size, built));
        Some(built)
    }

    fn spent(&self) -> u32 {
        u32::try_from(self.asked.len()).unwrap_or(u32::MAX)
    }
}

/// Whether `hi` stands within [`LIMIT_TOLERANCE`] of `lo`.
fn close(lo: f64, hi: f64) -> bool {
    hi <= lo * (1.0 + LIMIT_TOLERANCE)
}

/// Bisect between `lo`, which built, and `hi`, which was refused, until
/// they are [`close`] or the probes are spent.
fn bisect<F: FnMut(f64) -> bool>(probes: &mut Probes<F>, mut lo: f64, mut hi: f64) -> (f64, f64) {
    while !close(lo, hi) {
        let mid = f64::midpoint(lo, hi);
        match probes.ask(mid) {
            Some(true) => lo = mid,
            Some(false) => hi = mid,
            None => break,
        }
    }
    (lo, hi)
}

/// Search down from `hi`, refused, halving until a size above `floor`
/// builds, then bisect: the boundary below `hi`, or None with the smallest
/// size refused where none builds.
fn below<F: FnMut(f64) -> bool>(probes: &mut Probes<F>, hi: f64, floor: f64) -> (Option<f64>, f64) {
    let mut hi = hi;
    let mut t = hi * 0.5;
    loop {
        if t <= floor {
            return (None, hi);
        }
        match probes.ask(t) {
            Some(true) => break,
            Some(false) => {
                hi = t;
                t *= 0.5;
            }
            None => return (None, hi),
        }
    }
    let (lo, hi) = bisect(probes, t, hi);
    (Some(lo), hi)
}

/// The shared confirmation every largest-size reading ends in: the
/// boundary between sizes that build and sizes refused, by `builds`,
/// started from the `seat` reading where there is one, sizes at or below
/// `floor` never asked. See [`blend_limit`] for the search.
pub(crate) fn confirm(
    seat: Option<Option<f64>>,
    up_to: f64,
    floor: f64,
    builds: impl FnMut(f64) -> bool,
) -> BlendLimit {
    let mut probes = Probes {
        builds,
        asked: Vec::new(),
    };
    let seat_size = seat.flatten();
    // Where the seat holds no size, the smallest worth asking.
    let start = match seat {
        Some(None) => (floor * 2.0).min(up_to),
        Some(Some(s)) => s.min(up_to).max(floor * 2.0).min(up_to),
        None => up_to,
    };
    let (mut lo, mut hi): (Option<f64>, Option<f64>) = match probes.ask(start) {
        Some(true) => {
            // The seat can stand short of what builds (a bevel notching
            // past its face's far side): step just past it, then double.
            let mut lo = start;
            let mut hi = None;
            let mut t = (start * (1.0 + LIMIT_TOLERANCE)).min(up_to);
            while lo < up_to {
                match probes.ask(t) {
                    Some(true) => {
                        lo = t;
                        t = (t * 2.0).min(up_to);
                    }
                    Some(false) => {
                        hi = Some(t);
                        break;
                    }
                    None => break,
                }
            }
            match hi {
                Some(h) => {
                    let (l, h) = bisect(&mut probes, lo, h);
                    (Some(l), Some(h))
                }
                None => (Some(lo), None),
            }
        }
        Some(false) => {
            let (l, h) = below(&mut probes, start, floor);
            (l, Some(h))
        }
        None => (None, None),
    };
    let mut builds_again = None;
    if let Some(l) = lo {
        // A refusal below the answer: the sizes that build break there.
        let samples: Vec<f64> = (1..=LIMIT_SAMPLES)
            .map(|k| l * f64::from(k) / f64::from(LIMIT_SAMPLES + 1))
            .filter(|s| *s > floor)
            .collect();
        let mut first_refused = None;
        let mut built_below = None;
        for s in samples {
            match probes.ask(s) {
                Some(true) if first_refused.is_none() => built_below = Some(s),
                Some(false) if first_refused.is_none() => first_refused = Some(s),
                _ => {}
            }
        }
        if let Some(r) = first_refused {
            builds_again = Some(l);
            let (l2, h2) = match built_below {
                Some(b) => {
                    let (l2, h2) = bisect(&mut probes, b, r);
                    (Some(l2), h2)
                }
                None => below(&mut probes, r, floor),
            };
            lo = l2;
            hi = Some(h2);
        }
    }
    if let (Some(h), None) = (hi, builds_again) {
        // A size building again above the first refusal.
        if h < up_to {
            for k in 1..=LIMIT_SAMPLES {
                let s = h * (up_to / h).powf(f64::from(k) / f64::from(LIMIT_SAMPLES));
                if probes.ask(s) == Some(true) {
                    builds_again = Some(s);
                    break;
                }
            }
        }
    }
    BlendLimit {
        size: lo,
        refused: hi,
        builds_again,
        seat: seat_size,
        probes: probes.spent(),
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    fn run(seat: Option<Option<f64>>, up_to: f64, builds: impl Fn(f64) -> bool) -> BlendLimit {
        confirm(seat, up_to, 1e-7, builds)
    }

    #[test]
    fn an_exact_seat_is_confirmed_in_two_probes_and_the_samples() {
        let found = run(Some(Some(2.0)), 10.0, |s| s <= 2.0);
        assert_eq!(found.size, Some(2.0));
        assert!(close(2.0, found.refused.unwrap()));
        assert_eq!(found.builds_again, None);
        assert_eq!(found.probes, 2 + 2 * LIMIT_SAMPLES);
    }

    #[test]
    fn a_seat_short_of_what_builds_is_searched_past() {
        let found = run(Some(Some(2.0)), 100.0, |s| s <= 7.3);
        let size = found.size.unwrap();
        assert!(size <= 7.3 && close(size, 7.3), "{size}");
        assert!(found.refused.unwrap() > 7.3);
    }

    #[test]
    fn a_seat_past_what_builds_is_searched_below() {
        let found = run(Some(Some(2.0)), 100.0, |s| s <= 0.31);
        let size = found.size.unwrap();
        assert!(size <= 0.31 && close(size, 0.31), "{size}");
    }

    #[test]
    fn the_bound_building_answers_the_bound() {
        let found = run(Some(Some(0.5)), 0.5, |_| true);
        assert_eq!(found.size, Some(0.5));
        assert_eq!(found.refused, None);
    }

    #[test]
    fn nothing_building_answers_none_within_the_cap() {
        let found = run(None, 100.0, |_| false);
        assert_eq!(found.size, None);
        assert!(found.probes <= LIMIT_PROBES);
    }

    #[test]
    fn a_size_building_again_above_a_refusal_is_reported() {
        // Refused only between 2 and 3: the answer is the first boundary.
        let found = run(Some(Some(2.0)), 10.0, |s| !(s > 2.0 && s < 3.0));
        assert_eq!(found.size, Some(2.0));
        assert!(found.builds_again.unwrap() >= 3.0);
    }

    #[test]
    fn a_refusal_below_the_seat_s_boundary_moves_the_answer_down() {
        // Refused between 0.9 and 1.1, building again up to 2: the answer
        // is 0.9, and 2 is reported as building again.
        let found = run(Some(Some(2.0)), 10.0, |s| {
            s <= 0.9 || (1.1..=2.0).contains(&s)
        });
        let size = found.size.unwrap();
        assert!(size <= 0.9 && close(size, 0.9), "{size}");
        assert_eq!(found.builds_again, Some(2.0));
    }
}
