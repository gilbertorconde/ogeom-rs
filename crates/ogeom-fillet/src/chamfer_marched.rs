//! The chamfer on a seat only the marching machinery speaks: a crease
//! between two faces with no closed form (a bore's rim where it leaves a
//! drum, a spline's crease along a wall).
//!
//! No ball is rolled. At each station along the crease the section plane
//! stands square to it, and each face meets that plane in a section curve
//! leaving the crease. A leg runs its distance along that curve, the
//! geodesic of the face in the section: the curve is walked from the
//! crease point in short steps, each step put back exactly on the face and
//! the plane, the chords summed. `Angle` walks the named face's leg, then
//! the other face's section until it crosses the line leaving the named
//! contact at the angle. The bevel is the ruled band through the chords
//! between the contacts, fitted round the loop, and the wedge is the
//! marched fillet's with the chord in place of the ball's arc.

use crate::chamfer::Chamfer;
use crate::march::{BlendStop, MarchedBlend, Sides};
use crate::marched::{Ball, Chord, MarchedSeat, chart_of, crease_convexity, marched_seat};
use crate::support::keeping_face_meshes;
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Point, Point2, Vector};
use ogeom_topo::{Model, Shape};

/// How many steps each leg's walk along its section curve takes.
const WALK_STEPS: u32 = 64;

/// How the legs are given: a distance along each face, or a distance along
/// one face and the bevel's angle off it.
#[derive(Clone, Copy)]
enum Law {
    /// The legs along the first host and the second.
    Distances(f64, f64),
    /// The leg along the host `named` (0 or 1), and the angle the bevel
    /// leaves that face at in the section.
    Angle {
        named: usize,
        distance: f64,
        angle: f64,
    },
}

/// One host as the walks read it: its surface and outward sign.
struct Host<'a> {
    surface: &'a SurfaceGeometry,
    sign: f64,
}

impl Host<'_> {
    /// The outward unit normal at `(u, v)`.
    fn normal(&self, u: f64, v: f64, tol: Tolerances) -> Option<Vector> {
        crate::march::unit_normal(self.surface, u, v, tol).map(|n| n * self.sign)
    }

    /// The point of the face on the section plane through `g` square to
    /// `along`, a distance `step` past `from` along `dir`, solved from
    /// `(u, v)`: the plane and the advance along `dir` held exactly.
    fn advance(
        &self,
        (g, along): (Point, Vector),
        from: Point,
        dir: Vector,
        step: f64,
        (mut u, mut v): (f64, f64),
        tol: Tolerances,
    ) -> Option<(f64, f64)> {
        let target = from + dir * step;
        for _ in 0..40 {
            let p = self.surface.point_at(u, v, tol).ok()?;
            let (su, sv) = self.surface.d1_at(u, v, tol).ok()?;
            let f = [(p - g).dot(along), (p - target).dot(dir)];
            if f[0].abs() <= tol.confusion() * 1e-3 && f[1].abs() <= tol.confusion() * 1e-3 {
                return Some((u, v));
            }
            let (a, b, c, d) = (su.dot(along), sv.dot(along), su.dot(dir), sv.dot(dir));
            let det = a.mul_add(d, -(b * c));
            if det.abs() <= f64::EPSILON {
                return None;
            }
            u -= d.mul_add(f[0], -(b * f[1])) / det;
            v -= a.mul_add(f[1], -(c * f[0])) / det;
        }
        None
    }
}

/// A walk along one host's section curve from the crease.
struct Walk<'a, 'h> {
    host: &'a Host<'h>,
    plane: (Point, Vector),
    at: Point,
    uv: (f64, f64),
    dir: Vector,
    walked: f64,
}

impl Walk<'_, '_> {
    /// Step `step` on along the section, keeping the walk's heading.
    fn step(&mut self, step: f64, tol: Tolerances) -> Option<()> {
        let uv = self
            .host
            .advance(self.plane, self.at, self.dir, step, self.uv, tol)?;
        let p = self.host.surface.point_at(uv.0, uv.1, tol).ok()?;
        self.walked += p.distance(self.at);
        let normal = self.host.normal(uv.0, uv.1, tol)?;
        let heading = self.plane.1.cross(normal);
        let m = heading.magnitude();
        if m <= tol.angular() {
            return None;
        }
        let heading = heading / m;
        self.dir = if heading.dot(self.dir) >= 0.0 {
            heading
        } else {
            -heading
        };
        self.at = p;
        self.uv = uv;
        Some(())
    }

    /// Walk on to `length` from the crease, the last step trimmed so the
    /// summed chords come to it.
    fn to(&mut self, length: f64, tol: Tolerances) -> Option<()> {
        let stride = length / f64::from(WALK_STEPS);
        while length - self.walked > stride * 0.5 {
            self.step(stride.min(length - self.walked), tol)?;
        }
        let rest = length - self.walked;
        if rest.abs() > tol.confusion() * 1e-3 {
            self.step(rest, tol)?;
        }
        Some(())
    }
}

/// The band of `spec`'s legs on a marched seat, crossing each station by
/// `chord`: an asymmetric round's conic. With `laws`, a distance law for
/// `spec`'s named face and for the other, each leg runs its law's distance
/// at its station, the law laid along the crease by arc length from the
/// edge's start as a variable radius is: held level past an open edge's
/// ends, round a closed one. `mates`, where given, are the request's
/// edges by their ends, the edge's own at `index`: a mate continuing the
/// crease at an end stops the band flush there.
#[allow(
    clippy::too_many_arguments,
    reason = "the seat, the legs and their laws"
)]
pub(crate) fn marched_bevel_with(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    spec: &Chamfer,
    laws: Option<&[crate::law::Profile; 2]>,
    given: &[Shape],
    mates: Option<(usize, &[crate::fillet::Mate])>,
    chord: Chord,
    tol: Tolerances,
) -> OgeomResult<Built> {
    // What the refusals name: the bevel, or the round its chord stands for.
    let what = "an asymmetric round";
    let reach = match spec {
        Chamfer::Symmetric(d) => *d,
        Chamfer::Distances {
            on_face, on_other, ..
        } => on_face.max(*on_other),
        Chamfer::Angle { distance, .. } => *distance,
    };
    if !reach.is_finite() || reach <= tol.confusion() {
        ogeom_bail!(Construction, "a chamfer of {reach} cuts nothing");
    }
    let MarchedSeat {
        faces,
        hosts: [(first, sign_first), (second, sign_second)],
        guide,
        guide_range,
        edge_range,
        closed,
        loops,
        ..
    } = marched_seat(model, solid, edge, reach, tol)?;
    // An open seat runs out past the edge's ends, along the crease
    // continued where it ends with the edge.
    let (guide, guide_range) = if loops {
        (guide, guide_range)
    } else {
        crate::marched::run_out_guide(guide, guide_range, closed, [&first, &second], reach, tol)?
    };
    let side = |face: &Shape| -> OgeomResult<usize> {
        if faces[0].is_same(face) {
            Ok(0)
        } else if faces[1].is_same(face) {
            Ok(1)
        } else {
            ogeom_bail!(
                Construction,
                "the named face does not meet the edge being chamfered"
            )
        }
    };
    let law = match spec {
        Chamfer::Symmetric(d) => Law::Distances(*d, *d),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => {
            if side(face)? == 0 {
                Law::Distances(*on_face, *on_other)
            } else {
                Law::Distances(*on_other, *on_face)
            }
        }
        Chamfer::Angle {
            face,
            distance,
            angle,
        } => {
            if !angle.is_finite() || *angle <= tol.angular() {
                ogeom_bail!(
                    Construction,
                    "a chamfer at an angle of {angle} cuts nothing"
                );
            }
            Law::Angle {
                named: side(face)?,
                distance: *distance,
                angle: *angle,
            }
        }
    };
    if let Law::Distances(a, b) = law {
        for d in [a, b] {
            if !d.is_finite() || d <= tol.confusion() {
                ogeom_bail!(Construction, "a chamfer of {d} cuts nothing");
            }
        }
    }
    let convex = crease_convexity(
        model,
        &faces[0],
        [(&first, sign_first), (&second, sign_second)],
        &guide,
        edge_range,
        reach,
        tol,
    )?;
    let hosts = [
        Host {
            surface: &first,
            sign: sign_first,
        },
        Host {
            surface: &second,
            sign: sign_second,
        },
    ];
    // Each law in the seat's face order, laid along the crease.
    let along = match laws {
        Some(laws) => {
            let named = match spec {
                Chamfer::Distances { face, .. } => side(face)?,
                _ => 0,
            };
            let ordered = if named == 0 {
                [&laws[0], &laws[1]]
            } else {
                [&laws[1], &laws[0]]
            };
            Some(laid_along(
                model,
                edge,
                &guide,
                guide_range,
                edge_range,
                loops,
                ordered,
                tol,
            )?)
        }
        None => None,
    };
    let station = |w: f64, near: Option<[f64; 4]>| -> Option<[f64; 5]> {
        let law = match &along {
            Some([a, b]) => Law::Distances(a.at(w, tol).ok()?.0, b.at(w, tol).ok()?.0),
            None => law,
        };
        station_at(&hosts, &guide, law, convex, w, near, tol)
    };

    // Stations evenly over the crease, as many as the band's fit needs:
    // round the whole loop, or over the run around the edge on an open
    // seat, which stops where a station no longer seats.
    let aim = crate::marched::band_fit_target(tol);
    let (lo, hi) = guide_range;
    let mid = f64::midpoint(edge_range.0, edge_range.1);
    let mut count: u32 = 64;
    let blend = loop {
        let mut found: Vec<Option<[f64; 5]>> = Vec::with_capacity(count as usize);
        let mut near = None;
        for k in 0..count {
            let w = (hi - lo).mul_add(f64::from(k) / f64::from(count), lo);
            let x = station(w, near);
            if let Some(x) = x {
                near = Some([x[0], x[1], x[2], x[3]]);
            }
            found.push(x);
        }
        if !loops {
            // The run that holds the edge's middle, the rest dropped.
            let at = |k: usize| {
                let k = u32::try_from(k).unwrap_or(u32::MAX);
                (hi - lo).mul_add(f64::from(k) / f64::from(count), lo)
            };
            let centre = (0..found.len())
                .min_by(|a, b| (at(*a) - mid).abs().total_cmp(&(at(*b) - mid).abs()))
                .unwrap_or(0);
            let mut first_kept = centre;
            while first_kept > 0 && found[first_kept - 1].is_some() {
                first_kept -= 1;
            }
            let mut last_kept = centre;
            while last_kept + 1 < found.len() && found[last_kept + 1].is_some() {
                last_kept += 1;
            }
            for (k, x) in found.iter_mut().enumerate() {
                if k < first_kept || k > last_kept {
                    *x = None;
                }
            }
        }
        let mut blend = MarchedBlend {
            spine: Vec::new(),
            on_first: Vec::new(),
            on_second: Vec::new(),
            touch_first: Vec::new(),
            touch_second: Vec::new(),
            along: Vec::new(),
            sides: Sides {
                first: 1,
                second: 1,
            },
            stopped: if loops {
                BlendStop::Closed
            } else {
                BlendStop::RanPastTheGuide
            },
        };
        for (k, x) in found.into_iter().enumerate() {
            let k = u32::try_from(k).unwrap_or(u32::MAX);
            let w = (hi - lo).mul_add(f64::from(k) / f64::from(count), lo);
            let Some(x) = x else {
                if loops {
                    ogeom_bail!(
                        Construction,
                        "{what}'s leg runs off its face's section at {:?}",
                        guide.point_at(w, tol)?
                    );
                }
                continue;
            };
            let p1 = first.point_at(x[0], x[1], tol)?;
            let p2 = second.point_at(x[2], x[3], tol)?;
            blend.spine.push(p1.midpoint(p2));
            blend.on_first.push((x[0], x[1]));
            blend.on_second.push((x[2], x[3]));
            blend.touch_first.push(p1);
            blend.touch_second.push(p2);
            blend.along.push(w);
        }
        if blend.len() < 8 {
            ogeom_bail!(Construction, "{what}'s sections do not seat along the edge");
        }
        let met = match chord {
            Chord::Conic(weight) => crate::marched::conic_fit_met(
                &blend,
                [&first, &second],
                &guide,
                weight,
                loops,
                aim,
                tol,
            )?,
        };
        if met || count >= 1024 {
            break blend;
        }
        count *= 2;
    };
    if !loops && (blend.along[0] > edge_range.0 || blend.along[blend.len() - 1] < edge_range.1) {
        ogeom_bail!(
            Construction,
            "{what}'s sections do not seat all along the edge"
        );
    }

    // Each contact over the edge's own stretch stands on its face, or the
    // wedge would cut past the face; past an open edge's ends the run-out
    // leaves the faces by design. A contact on a face as the request named
    // it (`given`), where an earlier bevel of the request has cut that face
    // back, stands in what that bevel took, as two bevels' mitre does.
    keeping_face_meshes(|| -> OgeomResult<()> {
        let deflection = ogeom_mesh::Deflection {
            chord: reach * 1e-2,
            ..ogeom_mesh::Deflection::default()
        };
        for (touches, face) in [
            (&blend.touch_first, &faces[0]),
            (&blend.touch_second, &faces[1]),
        ] {
            for (p, w) in touches.iter().zip(&blend.along) {
                if !loops && (*w < edge_range.0 || *w > edge_range.1) {
                    continue;
                }
                let off = |face: &Shape| -> OgeomResult<bool> {
                    Ok(
                        ogeom_algo::classify_on_face(model, face, *p, deflection, tol)?
                            == ogeom_algo::Containment::Out,
                    )
                };
                let mut on_given = false;
                for named in given {
                    if !off(named)? {
                        on_given = true;
                        break;
                    }
                }
                if off(face)? && !on_given {
                    ogeom_bail!(Construction, "{what}'s contact at {p:?} runs past its face");
                }
            }
        }
        Ok(())
    })?;

    let ruled = |w: f64, near: [f64; 4]| station(w, Some(near));
    if !loops {
        return crate::marched::open_ruled_wedge(
            model,
            solid,
            edge,
            blend,
            &guide,
            edge_range,
            [(&first, sign_first), (&second, sign_second)],
            [&faces[0], &faces[1]],
            mates,
            reach,
            &ruled,
            chord,
            convex,
            tol,
        );
    }
    crate::marched::closed_band_wedge(
        model,
        solid,
        Some(edge),
        blend,
        &guide,
        guide_range,
        loops,
        [(&first, sign_first), (&second, sign_second)],
        Ball::Ruled(&ruled, chord),
        convex,
        None,
        tol,
    )
}

/// `laws` laid along a marched seat's crease `guide` by arc length from
/// `edge`'s start, the way the edge runs: over the edge's own window of the
/// crease on an open seat, held level past it, or round the whole loop the
/// crease closes, from the edge's start.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where an open law is laid over one arc of a crease that closes through
/// other edges, and as reading the edge and measuring the crease.
#[allow(clippy::too_many_arguments, reason = "the seat's crease and the laws")]
fn laid_along(
    model: &Model,
    edge: &Shape,
    guide: &Curve,
    guide_range: (f64, f64),
    edge_range: (f64, f64),
    loops: bool,
    laws: [&crate::law::Profile; 2],
    tol: Tolerances,
) -> OgeomResult<[crate::law::RadiusAlong; 2]> {
    use crate::law::{RadiusAlong, Ruler};
    use ogeom_topo::Orientation;
    if laws.iter().any(|law| loops && !law.closed()) {
        ogeom_bail!(
            Construction,
            "an asymmetric round whose distances vary over one arc of a crease that \
             closes through other edges is still owed; see docs/CAPABILITIES.md, \
             fillet.asymmetric-rounds"
        );
    }
    let reversed = edge.orientation() == Orientation::Reversed;
    let (ruler, origin, reversed) = if loops {
        let ruler = Ruler::new(guide.clone(), guide_range, true, tol)?;
        let (stored, (lo, hi)) = crate::support::edge_curve(model, edge, tol)?;
        let (start, heading) = if reversed {
            (stored.point_at(hi, tol)?, -stored.d1_at(hi, tol)?)
        } else {
            (stored.point_at(lo, tol)?, stored.d1_at(lo, tol)?)
        };
        let at = ogeom_algo::project_on_curve(guide, start, 256, tol)?.parameter;
        let against = guide.d1_at(at, tol)?.dot(heading) < 0.0;
        let origin = ruler.fraction(at, tol)?.0;
        (ruler, origin, against)
    } else {
        let ruler = Ruler::new(guide.clone(), edge_range, false, tol)?;
        (ruler, if reversed { 1.0 } else { 0.0 }, reversed)
    };
    Ok(laws.map(|law| RadiusAlong::new(law.clone(), ruler.clone(), origin, reversed)))
}

/// A chamfer's contacts on a marched seat at stations over the edge's own
/// stretch, as the marched bevel reads them: each station's crease
/// parameter and its contact on each of the seat's two faces.
pub(crate) struct CreaseContacts {
    /// The two faces either side of the edge.
    pub faces: [Shape; 2],
    /// Whether the edge closes on itself: the last station is the first
    /// again.
    pub closed: bool,
    /// The stations' crease parameters, increasing.
    pub along: Vec<f64>,
    /// Each station's contact on the first face and on the second.
    pub points: [Vec<Point>; 2],
    /// Whether each station's faces, read in its section plane, leave
    /// their contacts along tangents that meet in front of the crease: the
    /// corner a conic section between the contacts stands in.
    pub shouldered: Vec<bool>,
    /// Each station's point on the crease.
    pub crease: Vec<Point>,
}

/// [`CreaseContacts`] of `spec` on `edge` of `solid` at `count + 1`
/// stations evenly over the edge's stretch of its crease, the legs walked
/// as [`marched_bevel_with`] walks them; with `whole`, once round the whole loop
/// of a crease that comes back to its start, from its point nearest
/// `whole`, read as closed (the pieces of one loop the boolean left, taken
/// as one edge), the stations' parameters run on past the crease's own join
/// rather than wrapped. `None` where a station does not seat, or `whole` is
/// asked of a crease that does not loop.
///
/// # Errors
///
/// As [`marched_bevel_with`] reading the seat and the law.
pub(crate) fn crease_contacts(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    spec: &Chamfer,
    count: u32,
    whole: Option<Point>,
    tol: Tolerances,
) -> OgeomResult<Option<CreaseContacts>> {
    let reach = match spec {
        Chamfer::Symmetric(d) => *d,
        Chamfer::Distances {
            on_face, on_other, ..
        } => on_face.max(*on_other),
        Chamfer::Angle { distance, .. } => *distance,
    };
    if !reach.is_finite() || reach <= tol.confusion() || count == 0 {
        return Ok(None);
    }
    let MarchedSeat {
        faces,
        hosts: [(first, sign_first), (second, sign_second)],
        guide,
        guide_range,
        edge_range,
        closed,
        loops,
        ..
    } = marched_seat(model, solid, edge, reach, tol)?;
    let period = guide_range.1 - guide_range.0;
    // The piece's own stretch, where the crease's convexity is read.
    let own_range = edge_range;
    let (edge_range, closed) = match whole {
        Some(start) => {
            if !loops {
                return Ok(None);
            }
            let w0 = ogeom_algo::project_on_curve(&guide, start, 256, tol)?.parameter;
            ((w0, w0 + period), true)
        }
        None => (edge_range, closed),
    };
    let wrapped = |w: f64| {
        if whole.is_some() {
            guide_range.0 + (w - guide_range.0).rem_euclid(period)
        } else {
            w
        }
    };

    let side = |face: &Shape| -> Option<usize> {
        if faces[0].is_same(face) {
            Some(0)
        } else if faces[1].is_same(face) {
            Some(1)
        } else {
            None
        }
    };
    let law = match spec {
        Chamfer::Symmetric(d) => Law::Distances(*d, *d),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => match side(face) {
            Some(0) => Law::Distances(*on_face, *on_other),
            Some(_) => Law::Distances(*on_other, *on_face),
            None => return Ok(None),
        },
        Chamfer::Angle {
            face,
            distance,
            angle,
        } => {
            let Some(named) = side(face) else {
                return Ok(None);
            };
            if !angle.is_finite() || *angle <= tol.angular() {
                return Ok(None);
            }
            Law::Angle {
                named,
                distance: *distance,
                angle: *angle,
            }
        }
    };
    let convex = crease_convexity(
        model,
        &faces[0],
        [(&first, sign_first), (&second, sign_second)],
        &guide,
        own_range,
        reach,
        tol,
    )?;
    let hosts = [
        Host {
            surface: &first,
            sign: sign_first,
        },
        Host {
            surface: &second,
            sign: sign_second,
        },
    ];
    let (lo, hi) = edge_range;
    let mut along = Vec::with_capacity(count as usize + 1);
    let mut points = [Vec::new(), Vec::new()];
    let mut shouldered = Vec::with_capacity(count as usize + 1);
    let mut crease = Vec::with_capacity(count as usize + 1);
    let mut near = None;
    for k in 0..=count {
        let w = (hi - lo).mul_add(f64::from(k) / f64::from(count), lo);
        let found = station_at(&hosts, &guide, law, convex, wrapped(w), near, tol);
        let Some(x) = found else {
            return Ok(None);
        };
        near = Some([x[0], x[1], x[2], x[3]]);
        along.push(w);
        let p = [
            first.point_at(x[0], x[1], tol)?,
            second.point_at(x[2], x[3], tol)?,
        ];
        let n = [
            hosts[0].normal(x[0], x[1], tol),
            hosts[1].normal(x[2], x[3], tol),
        ];
        let on_crease = guide.point_at(wrapped(w), tol)?;
        let tangent = guide.d1_at(wrapped(w), tol)?;
        shouldered.push(match n {
            [Some(n0), Some(n1)] => shoulder_in_front(p, [n0, n1], on_crease, tangent, tol),
            _ => false,
        });
        points[0].push(p[0]);
        points[1].push(p[1]);
        crease.push(on_crease);
    }
    Ok(Some(CreaseContacts {
        faces,
        closed,
        along,
        points,
        shouldered,
        crease,
    }))
}

/// Whether the lines leaving the contacts `p` square to the crease's
/// `tangent` along their faces (of normals `n`), each back towards the
/// `crease` point, meet in front of both contacts and are not parallel:
/// the shoulder a conic section's middle control point stands on. A leg
/// turning a quarter round a drum leaves along a line parallel to the
/// other face's, its shoulder gone to infinity.
fn shoulder_in_front(
    p: [Point; 2],
    n: [Vector; 2],
    crease: Point,
    tangent: Vector,
    tol: Tolerances,
) -> bool {
    let length = tangent.magnitude();
    if length <= tol.confusion() {
        return false;
    }
    let t = tangent / length;
    let mut d = [Vector::ZERO; 2];
    for i in 0..2 {
        let along = t.cross(n[i]);
        let m = along.magnitude();
        if m <= tol.angular() {
            return false;
        }
        d[i] = if along.dot(crease - p[i]) < 0.0 {
            -along / m
        } else {
            along / m
        };
    }
    // p0 + a d0 = p1 + b d1, in the section plane.
    let cross = d[0].cross(d[1]);
    let sin = cross.magnitude();
    if sin <= tol.angular() {
        return false;
    }
    let gap = p[1] - p[0];
    let a = gap.cross(d[1]).dot(cross) / (sin * sin);
    let b = gap.cross(d[0]).dot(cross) / (sin * sin);
    a > tol.confusion() && b > tol.confusion()
}

/// The station at guide parameter `w`: both contacts' parameters and `w`.
fn station_at(
    hosts: &[Host<'_>; 2],
    guide: &Curve,
    law: Law,
    convex: bool,
    w: f64,
    near: Option<[f64; 4]>,
    tol: Tolerances,
) -> Option<[f64; 5]> {
    let g = guide.point_at(w, tol).ok()?;
    let along = guide.d1_at(w, tol).ok()?;
    let m = along.magnitude();
    if m <= tol.confusion() {
        return None;
    }
    let along = along / m;
    // Each walk starts at the crease point, its chart position read near
    // the last station's contact.
    let start = |i: usize| -> Option<(f64, f64)> {
        let prev = near.map(|n| Point2::new(n[2 * i], n[2 * i + 1]));
        let uv = chart_of(hosts[i].surface, g, prev, tol).ok()?;
        Some((uv.x, uv.y))
    };
    let normals = [
        {
            let (u, v) = start(0)?;
            hosts[0].normal(u, v, tol)?
        },
        {
            let (u, v) = start(1)?;
            hosts[1].normal(u, v, tol)?
        },
    ];
    // Into face `i` from the crease, square to it in the face: away from
    // the other face's tangent plane at a convex crease, toward it at a
    // concave one.
    let walk = |i: usize| -> Option<Walk<'_, '_>> {
        let heading = along.cross(normals[i]);
        let hm = heading.magnitude();
        if hm <= tol.angular() {
            return None;
        }
        let heading = heading / hm;
        let lean = heading.dot(normals[1 - i]);
        let into = if (lean < 0.0) == convex {
            heading
        } else {
            -heading
        };
        Some(Walk {
            host: &hosts[i],
            plane: (g, along),
            at: g,
            uv: start(i)?,
            dir: into,
            walked: 0.0,
        })
    };
    let (uv0, uv1) = match law {
        Law::Distances(a, b) => {
            let mut w0 = walk(0)?;
            w0.to(a, tol)?;
            let mut w1 = walk(1)?;
            w1.to(b, tol)?;
            (w0.uv, w1.uv)
        }
        Law::Angle {
            named,
            distance,
            angle,
        } => {
            let mut leg = walk(named)?;
            leg.to(distance, tol)?;
            // Back toward the crease along the section, turned by the angle
            // into the wedge's side, both in the section plane: the
            // material at a convex crease, the open notch at a concave one.
            let normal = hosts[named].normal(leg.uv.0, leg.uv.1, tol)?;
            let square = along.cross(leg.dir);
            let square = square / square.magnitude();
            let wedge_side = if convex { -normal } else { normal };
            let side = if square.dot(wedge_side) >= 0.0 {
                square
            } else {
                -square
            };
            let line = -leg.dir * angle.cos() + side * angle.sin();
            let across = along.cross(line);
            let offset = |p: Point| (p - leg.at).dot(across);
            let mut other = walk(1 - named)?;
            let stride = distance / f64::from(WALK_STEPS);
            let bound = {
                let (u, v) = other.host.surface.domain();
                let corner = |a: f64, b: f64| other.host.surface.point_at(a, b, tol).ok();
                match (corner(u.0, v.0), corner(u.1, v.1)) {
                    (Some(p), Some(q)) => p.distance(q).max(distance) * 4.0,
                    _ => distance * 64.0,
                }
            };
            let mut before = offset(other.at);
            loop {
                let (at, uv, dir, walked) = (other.at, other.uv, other.dir, other.walked);
                other.step(stride, tol)?;
                let now = offset(other.at);
                if now == 0.0 || now.signum() != before.signum() {
                    // The crossing lies within this step: bisect its
                    // length from where it began.
                    let (mut lo, mut hi) = (0.0, stride);
                    for _ in 0..60 {
                        let mid = f64::midpoint(lo, hi);
                        other.at = at;
                        other.uv = uv;
                        other.dir = dir;
                        other.walked = walked;
                        other.step(mid, tol)?;
                        if offset(other.at).signum() == before.signum() {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                    }
                    break;
                }
                before = now;
                if other.walked > bound {
                    return None;
                }
            }
            if named == 0 {
                (leg.uv, other.uv)
            } else {
                (other.uv, leg.uv)
            }
        }
    };
    Some([uv0.0, uv0.1, uv1.0, uv1.1, w])
}
