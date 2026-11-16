//! Rounds of straight convex creases meeting at inside corners, carried on
//! past each other (ADR 0001, amendment "M5e, rounds carried on through an
//! inside corner").
//!
//! Where two rounded top creases meet at a corner that is reflex in plan
//! (an inside corner of a notch, a pocket's rim), each band runs on past
//! the corner as its own crease's round, and the solid keeps whatever lies
//! under either band: a point goes only where it lies above both bands
//! carried on, or above one and past the other's wall, where the corner is
//! open. Off square in plan, or with walls leaning, the bands then meet
//! along their own intersection and the corner's cell past both walls,
//! which no flush cap or patch at the vertex holds.
//!
//! Each crease's band carried on is its spandrel prism: its two faces'
//! wedge, far longer than the solid, less that wedge rounded by the engine
//! (an exact cylinder). Each crease's tool is its spandrel where, at each
//! end shared with another asked crease, it also lies past that crease's
//! wall or within its spandrel; the solid is cut by the tools (the one
//! boolean of amendment "M5a"). A crease end with no asked neighbour must
//! carry its band out of the solid there, so the band ends flush on the
//! faces it meets.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Axis, Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::curved::{Curved, read};
use super::fallback::{Sense, Tool, polytope};
use super::trimmed::{faces_of, outward};
use super::{Id, id};

/// One asked crease: its edge, its ends, its two faces with their outward
/// normals and offsets at its start, its radius, and its kind.
struct Crease {
    shape: Shape,
    ends: (Id, Id),
    a: Point,
    b: Point,
    faces: [usize; 2],
    planes: [(Vector, f64); 2],
    /// The largest radius it is rounded at.
    radius: f64,
    kind: Kind,
    /// A varying round's law over the crease, from its first vertex.
    law: Option<crate::law::Profile>,
}

/// A crease's kind: straight between two planes, or a rim, an arc about an
/// axis between a plane square to it and a coaxial cylinder or cone.
#[derive(Clone, Copy)]
enum Kind {
    Straight,
    /// The arc's circle (its curve), its range on it, and which face is
    /// the plane (the other is the wall).
    Rim {
        range: (f64, f64),
        top: usize,
    },
}

impl Crease {
    /// The direction the crease leaves vertex `v` in.
    fn leaving(&self, curved: &Curved, v: Id, tol: Tolerances) -> OgeomResult<Vector> {
        let start = self.ends.0 == v;
        let d = match self.kind {
            Kind::Straight => {
                if start {
                    self.b - self.a
                } else {
                    self.a - self.b
                }
            }
            Kind::Rim { range, .. } => {
                let curve = &curved.edges[&id(&self.shape)].curve;
                if start {
                    curve.d1_at(range.0, tol)?
                } else {
                    -curve.d1_at(range.1, tol)?
                }
            }
        };
        Ok(d / d.magnitude())
    }
}

/// What a crease is rounded with: a radius, or a law laid over it from the
/// end at the point given to the other.
#[derive(Debug, Clone)]
pub(crate) enum Spec {
    Round(f64),
    Law(crate::law::Profile, Point),
}

/// The rounds of `asked` on `solid` where every asked edge is a straight
/// convex crease between planes, two or more of them meet at inside
/// corners of a face they share, and every crease end met by no other
/// asked crease runs out of the solid; `None` otherwise, or where a tool
/// or the boolean will not build.
///
/// # Errors
///
/// As reading the solid and building the tools, other than refusals,
/// which decline.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    asked: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let specs: Vec<(Shape, Spec)> = asked
        .iter()
        .map(|(e, r)| (e.clone(), Spec::Round(*r)))
        .collect();
    rounds_with(model, solid, &specs, tol)
}

/// [`rounds`] with each crease's own [`Spec`]: a varying round's band
/// follows its law over its crease and holds level past its ends.
///
/// # Errors
///
/// As [`rounds`].
pub(crate) fn rounds_with(
    model: &mut Model,
    solid: &Shape,
    asked: &[(Shape, Spec)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let ang = tol.angular();
    if asked.len() < 2 {
        return Ok(None);
    }
    let varying = asked.iter().any(|(_, s)| matches!(s, Spec::Law(..)));
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    // In one order whatever the request's, so its booleans run alike.
    let mut asked = asked.to_vec();
    asked.sort_by_key(|(edge, _)| id(edge));
    let mut creases = Vec::with_capacity(asked.len());
    for (edge, spec) in &asked {
        let radius = match spec {
            Spec::Round(r) => *r,
            Spec::Law(profile, _) => profile.largest(),
        };
        let Some(mut crease) = crease_of(model, &curved, id(edge), edge, radius, tol)? else {
            return Ok(None);
        };
        // A straight crease's ball seated on both its faces: one running
        // off a face narrower than its setback is cut off there, not
        // carried on.
        if matches!(crease.kind, Kind::Straight) && !seated_along(model, solid, &crease, tol)? {
            return Ok(None);
        }
        if let Spec::Law(profile, from) = spec {
            // Laid from the crease's own first vertex.
            let reversed = from.distance(crease.b) < from.distance(crease.a);
            crease.law = Some(profile.scaled_by(1.0, reversed));
        }
        creases.push(crease);
    }
    // The extent of the solid, for the wedges and the walls' half-spaces.
    let (mut lo, mut hi) = (
        Vector::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        Vector::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in curved.points.values() {
        let v = *p - Point::ORIGIN;
        lo = Vector::new(lo.x.min(v.x), lo.y.min(v.y), lo.z.min(v.z));
        hi = Vector::new(hi.x.max(v.x), hi.y.max(v.y), hi.z.max(v.z));
    }
    let reach = (hi - lo).magnitude() + 4.0 * creases.iter().map(|c| c.radius).fold(0.0, f64::max);
    let mid = (lo + hi) * 0.5;
    let boxed = |planes: &mut Vec<(Vector, f64)>| {
        for axis in [Vector::X, Vector::Y, Vector::Z] {
            let m = axis.dot(mid);
            planes.push((axis, m + reach));
            planes.push((-axis, -(m - reach)));
        }
    };
    // Each crease's neighbour at each end: the asked crease sharing that
    // vertex, with the face both walk and the neighbour's other face (its
    // wall); the corner reflex in plan.
    let mut neighbours: Vec<[Option<(usize, usize, bool)>; 2]> = vec![[None, None]; creases.len()];
    // Each corner whose rails meet only past the face both creases share:
    // the creases, that face and the vertex (amendment "corners where a
    // round's rail misses").
    let mut rolled: Vec<(usize, usize, usize, Id)> = Vec::new();
    // The planes a crease's band ends flush on, at free ends past which the
    // material goes on.
    // With each, the half of space beyond the crease's middle from that end
    // (the plane there square to a straight crease, the meridian plane for a
    // rim), where the clip does not reach: the band carried on past the
    // other end is that end's to keep.
    let mut ends_at: Vec<Vec<[(Vector, f64); 2]>> = vec![Vec::new(); creases.len()];
    let mut met = false;
    for (i, c) in creases.iter().enumerate() {
        for (k, v) in [c.ends.0, c.ends.1].into_iter().enumerate() {
            let others: Vec<usize> = (0..creases.len())
                .filter(|j| *j != i && (creases[*j].ends.0 == v || creases[*j].ends.1 == v))
                .collect();
            match others[..] {
                [] => {
                    // Carried on past a free end, the band leaves the
                    // solid: a point just past the vertex, a little inside
                    // the wedge, stands outside it.
                    let at = if k == 0 { c.a } else { c.b };
                    let away = -c.leaving(&curved, v, tol)?;
                    let (Some(n0), Some(n1)) = (
                        outward(&curved, c.faces[0], at),
                        outward(&curved, c.faces[1], at),
                    ) else {
                        return Ok(None);
                    };
                    let inward = -(n0 + n1);
                    let inward = inward / inward.magnitude();
                    let probe = at + away * (c.radius * 0.5) + inward * (c.radius * 0.1);
                    if ogeom_algo::classify_in_solid_exact(model, solid, probe, tol)?
                        != ogeom_algo::Containment::Out
                    {
                        // Into the material: the band ends flush on the one
                        // plane the crease meets there, cut to the crease's
                        // side of it.
                        let at_v: Vec<usize> = (0..curved.faces.len())
                            .filter(|f| {
                                !c.faces.contains(f)
                                    && curved.faces[*f].loops.iter().flatten().any(|(x, _)| {
                                        let g = &curved.edges[x];
                                        g.ends.0 == v || g.ends.1 == v
                                    })
                            })
                            .collect();
                        let planes: Vec<Vector> = at_v
                            .iter()
                            .filter_map(|f| match curved.faces[*f].surface {
                                SurfaceGeometry::Plane(_) => outward(&curved, *f, at),
                                _ => None,
                            })
                            .collect();
                        let Some(&n) = planes.first() else {
                            return Ok(None);
                        };
                        if planes.len() != at_v.len()
                            || planes
                                .iter()
                                .any(|m| m.cross(n).magnitude() > ang || m.dot(n) < 0.0)
                        {
                            return Ok(None);
                        }
                        let d = n.dot(at - Point::ORIGIN);
                        let inside = c.leaving(&curved, v, tol)?;
                        let plane = if inside.dot(n) <= 0.0 {
                            (n, d)
                        } else {
                            (-n, -d)
                        };
                        let (middle, split) = match c.kind {
                            Kind::Straight => (c.a.midpoint(c.b), c.b - c.a),
                            Kind::Rim { range, .. } => {
                                let e = &curved.edges[&id(&c.shape)];
                                let t = 0.5 * (range.0 + range.1);
                                (e.curve.point_at(t, tol)?, e.curve.d1_at(t, tol)?)
                            }
                        };
                        let toward = if (at - middle).dot(split) >= 0.0 {
                            split
                        } else {
                            -split
                        };
                        let toward = toward / toward.magnitude();
                        ends_at[i].push([plane, (toward, toward.dot(middle - Point::ORIGIN))]);
                    }
                }
                [j] => {
                    let n = &creases[j];
                    let shared: Vec<usize> = c
                        .faces
                        .iter()
                        .copied()
                        .filter(|f| n.faces.contains(f))
                        .collect();
                    let [top] = shared[..] else {
                        return Ok(None);
                    };
                    let (Some(wall), Some(own)) = (
                        n.faces.iter().copied().find(|f| *f != top),
                        c.faces.iter().copied().find(|f| *f != top),
                    ) else {
                        return Ok(None);
                    };
                    let at = curved.points[&v];
                    let (Some(n_wall), Some(own_wall)) =
                        (outward(&curved, wall, at), outward(&curved, own, at))
                    else {
                        return Ok(None);
                    };
                    let (go, back) = (
                        c.leaving(&curved, v, tol)?.dot(n_wall),
                        n.leaving(&curved, v, tol)?.dot(own_wall),
                    );
                    // An inside corner: each crease, leaving the vertex,
                    // runs on the open side of the other's wall, and its band
                    // is kept only past that wall or within the other's
                    // band. A convex turn: each runs out past the vertex as
                    // it does alone. Anything else is not built here. Only
                    // an inside corner takes a constant request here:
                    // with convex turns alone the other routes close it, as
                    // one piece with no boolean where they can. A varying
                    // chain they do not take closes here at either.
                    if go > ang && back > ang {
                        let joined =
                            joined_at(model, solid, &curved, [c, n], [own, wall], top, v, tol)?;
                        if !joined {
                            if !rolls_round(&curved, [c, n], top, v, tol) {
                                return Ok(None);
                            }
                            if i < j {
                                rolled.push((i, j, top, v));
                            }
                        }
                        neighbours[i][k] = Some((j, wall, joined));
                        met = true;
                    } else if go < -ang && back < -ang {
                        met |= varying;
                    } else {
                        return Ok(None);
                    }
                }
                _ => return Ok(None),
            }
        }
    }
    if !met {
        return Ok(None);
    }
    let mut scratch = model.clone();
    // Each crease's spandrel carried on, its band's face named.
    let mut spandrels: Vec<(Shape, Vec<Shape>)> = Vec::with_capacity(creases.len());
    for c in &creases {
        let Some(spandrel) = spandrel(&mut scratch, &curved, c, reach, tol)? else {
            return Ok(None);
        };
        spandrels.push(spandrel);
    }
    // Each crease's wall's far side, with its spandrel.
    let mut open: Vec<Option<(Shape, Shape)>> = vec![None; creases.len()];
    let mut tools = Vec::with_capacity(creases.len());
    for (i, c) in creases.iter().enumerate() {
        let (mut tool, mut named) = spandrels[i].clone();
        for (j, wall, joined) in neighbours[i].iter().flatten() {
            let (past, either) = match &open[*j] {
                Some(s) => s.clone(),
                None => {
                    let past = match &curved.faces[*wall].surface {
                        SurfaceGeometry::Plane(_) => {
                            let at = creases[*j].a;
                            let Some(n) = outward(&curved, *wall, at) else {
                                return Ok(None);
                            };
                            let mut planes = vec![(-n, -n.dot(at - Point::ORIGIN))];
                            boxed(&mut planes);
                            polytope(&mut scratch, &planes, tol)?
                        }
                        _ => {
                            let mut planes = Vec::new();
                            boxed(&mut planes);
                            match polytope(&mut scratch, &planes, tol)? {
                                Some(bounds) => past_wall(
                                    &mut scratch,
                                    &curved,
                                    &creases[*j],
                                    *wall,
                                    &bounds,
                                    tol,
                                )?,
                                None => None,
                            }
                        }
                    };
                    let Some(past) = past else {
                        return Ok(None);
                    };
                    let Ok(both) = ogeom_bool::fuse(&mut scratch, &past, &spandrels[*j].0, tol)
                    else {
                        return Ok(None);
                    };
                    open[*j] = Some((past.clone(), both.shape.clone()));
                    (past, both.shape)
                }
            };
            if !*joined {
                // The band ends where its ball touches the sharp edge, in the
                // plane through it square to its own crease.
                let v = if creases[*j].ends.0 == c.ends.0 || creases[*j].ends.1 == c.ends.0 {
                    c.ends.0
                } else {
                    c.ends.1
                };
                let at = curved.points[&v];
                let on = c.leaving(&curved, v, tol)?;
                let mut planes = vec![(-on, -on.dot(at - Point::ORIGIN))];
                boxed(&mut planes);
                let Some(kept_side) = polytope(&mut scratch, &planes, tol)? else {
                    return Ok(None);
                };
                // Within its own side of that plane, and past the other's
                // wall: a section standing off its own wall reaches across a
                // notch narrower than its standoff otherwise.
                for bound in [kept_side, past] {
                    let Ok(kept) = ogeom_bool::common(&mut scratch, &tool, &bound, tol) else {
                        return Ok(None);
                    };
                    named = named
                        .iter()
                        .flat_map(|f| kept.history.trace(f).to_vec())
                        .collect();
                    tool = kept.shape;
                }
                continue;
            }
            let Ok(kept) = ogeom_bool::common(&mut scratch, &tool, &either, tol) else {
                return Ok(None);
            };
            named = named
                .iter()
                .flat_map(|f| kept.history.trace(f).to_vec())
                .collect();
            tool = kept.shape;
        }
        for [plane, (toward, offset)] in ends_at[i].clone() {
            // The crease's side of the end plane, or the half away from it.
            let mut front = vec![plane];
            boxed(&mut front);
            let mut away = vec![(toward, offset)];
            boxed(&mut away);
            let (Some(front), Some(away)) = (
                polytope(&mut scratch, &front, tol)?,
                polytope(&mut scratch, &away, tol)?,
            ) else {
                return Ok(None);
            };
            let Ok(either) = ogeom_bool::fuse(&mut scratch, &front, &away, tol) else {
                return Ok(None);
            };
            let Ok(kept) = ogeom_bool::common(&mut scratch, &tool, &either.shape, tol) else {
                return Ok(None);
            };
            named = named
                .iter()
                .flat_map(|f| kept.history.trace(f).to_vec())
                .collect();
            tool = kept.shape;
        }
        if explore_unique(&scratch, &tool, ShapeType::Solid)?.is_empty() || named.is_empty() {
            return Ok(None);
        }
        tools.push(Tool {
            solid: tool,
            sense: Sense::Cut,
            faces: named.into_iter().map(|f| (f, c.shape.clone())).collect(),
        });
    }
    // At each corner whose rails meet only past the face, the ball rolling
    // round the sharp edge between the two bands' ends.
    for &(i, j, top, v) in &rolled {
        let Some((turned, band)) = rolled_corner(
            &mut scratch,
            &curved,
            [&creases[i], &creases[j]],
            top,
            v,
            tol,
        )?
        else {
            return Ok(None);
        };
        let mut faces = Vec::with_capacity(2 * band.len());
        for f in &band {
            faces.push((f.clone(), creases[i].shape.clone()));
            faces.push((f.clone(), creases[j].shape.clone()));
        }
        tools.push(Tool {
            solid: turned,
            sense: Sense::Cut,
            faces,
        });
    }
    let edges: Vec<Shape> = creases.iter().map(|c| c.shape.clone()).collect();
    let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    match super::guards::solid_fits(&scratch, &built.shape, tol) {
        Ok(()) => {
            *model = scratch;
            Ok(Some(built))
        }
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}

/// Whether `crease`'s band, carried on past its end `k` (0 its start),
/// leaves the solid: a point just past the vertex, a little inside the
/// wedge, stands outside it.
fn runs_out(
    model: &Model,
    solid: &Shape,
    crease: &Crease,
    k: usize,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let (at, away) = if k == 0 {
        (crease.a, crease.a - crease.b)
    } else {
        (crease.b, crease.b - crease.a)
    };
    let away = away / away.magnitude();
    let inward = -(crease.planes[0].0 + crease.planes[1].0);
    let inward = inward / inward.magnitude();
    let probe = at + away * (crease.radius * 0.5) + inward * (crease.radius * 0.1);
    Ok(
        ogeom_algo::classify_in_solid_exact(model, solid, probe, tol)?
            == ogeom_algo::Containment::Out,
    )
}

/// How many stations along a straight crease its contacts are read at on
/// its faces, ends apart.
const CONTACT_STATIONS: u32 = 16;

/// Whether a straight crease's ball, at its largest radius, sits on both
/// its faces along it: each contact, and the way to it from the crease
/// across its face, on the solid's boundary between the crease's ends (a
/// contact past a groove beside the crease lands on the boundary again,
/// the way to it over the groove does not).
///
/// # Errors
///
/// As classifying the contacts.
fn seated_along(
    model: &Model,
    solid: &Shape,
    crease: &Crease,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let [(na, _), (nb, _)] = crease.planes;
    let centre = -(na + nb) * (crease.radius / (1.0 + na.dot(nb)));
    for i in 1..CONTACT_STATIONS {
        let at = crease.a + (crease.b - crease.a) * (f64::from(i) / f64::from(CONTACT_STATIONS));
        for n in [na, nb] {
            let contact = at + centre + n * crease.radius;
            for j in 1..=ACROSS_STATIONS {
                let p = at + (contact - at) * (f64::from(j) / f64::from(ACROSS_STATIONS));
                if ogeom_algo::classify_in_solid_exact(model, solid, p, tol)?
                    != ogeom_algo::Containment::On
                {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

/// How many points across its face, out to its contact, a straight
/// crease's seat is read at.
const ACROSS_STATIONS: u32 = 8;

/// The rounds of `asked` on `solid` where every asked edge is a straight
/// convex crease between planes, no two meeting, each carried on past both
/// ends out of the solid (a chord's crease running out across a drum): the
/// band ends flush on the faces it meets past each end, as a free end at an
/// inside corner does (ADR 0001, amendment "M5e"). Each tool is the
/// crease's spandrel carried on, which must meet the solid in one piece:
/// the band alone, nothing it reaches farther along its line. `None`
/// otherwise, or where a tool or the boolean will not build.
///
/// # Errors
///
/// As reading the solid and building the tools, other than refusals,
/// which decline.
pub(crate) fn carried_out(
    model: &mut Model,
    solid: &Shape,
    asked: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if asked.is_empty() {
        return Ok(None);
    }
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let mut creases: Vec<Crease> = Vec::with_capacity(asked.len());
    for (edge, radius) in asked {
        let Some(crease) = crease_of(model, &curved, id(edge), edge, *radius, tol)? else {
            return Ok(None);
        };
        if !matches!(crease.kind, Kind::Straight) {
            return Ok(None);
        }
        // No two meeting or sharing a face: two bands on one face are
        // checked against each other where they are built together.
        let shares = |c: &Crease| {
            [c.ends.0, c.ends.1]
                .iter()
                .any(|v| *v == crease.ends.0 || *v == crease.ends.1)
                || c.faces.iter().any(|f| crease.faces.contains(f))
        };
        if creases.iter().any(shares) {
            return Ok(None);
        }
        for k in 0..2 {
            if !runs_out(model, solid, &crease, k, tol)? {
                return Ok(None);
            }
        }
        // Each end on the solid's own faces, each meeting the crease's at an
        // angle: one running on tangent from them (a band an earlier round
        // left) is that round's to meet, not a face the band runs out across.
        for v in [crease.ends.0, crease.ends.1] {
            let at = curved.points[&v];
            let [(na, _), (nb, _)] = crease.planes;
            for (f, face) in curved.faces.iter().enumerate() {
                if crease.faces.contains(&f)
                    || !face.loops.iter().flatten().any(|(x, _)| {
                        let g = &curved.edges[x];
                        g.ends.0 == v || g.ends.1 == v
                    })
                {
                    continue;
                }
                let Some(n) = outward(&curved, f, at) else {
                    return Ok(None);
                };
                if n.cross(na).magnitude() <= tol.angular()
                    || n.cross(nb).magnitude() <= tol.angular()
                {
                    return Ok(None);
                }
            }
        }
        if !seated_along(model, solid, &crease, tol)? {
            return Ok(None);
        }
        creases.push(crease);
    }
    let (mut lo, mut hi) = (
        Vector::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        Vector::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in curved.points.values() {
        let v = *p - Point::ORIGIN;
        lo = Vector::new(lo.x.min(v.x), lo.y.min(v.y), lo.z.min(v.z));
        hi = Vector::new(hi.x.max(v.x), hi.y.max(v.y), hi.z.max(v.z));
    }
    let reach = (hi - lo).magnitude() + 4.0 * creases.iter().map(|c| c.radius).fold(0.0, f64::max);
    let mut scratch = model.clone();
    let mut tools = Vec::with_capacity(creases.len());
    for c in &creases {
        let Some((tool, named)) = spandrel(&mut scratch, &curved, c, reach, tol)? else {
            return Ok(None);
        };
        let Ok(taken) = ogeom_bool::common(&mut scratch, &tool, solid, tol) else {
            return Ok(None);
        };
        if explore_unique(&scratch, &taken.shape, ShapeType::Solid)?.len() != 1 {
            return Ok(None);
        }
        tools.push(Tool {
            solid: tool,
            sense: Sense::Cut,
            faces: named.into_iter().map(|f| (f, c.shape.clone())).collect(),
        });
    }
    let edges: Vec<Shape> = creases.iter().map(|c| c.shape.clone()).collect();
    let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    match super::guards::solid_fits(&scratch, &built.shape, tol) {
        Ok(()) => {
            *model = scratch;
            Ok(Some(built))
        }
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}

/// Whether two creases meeting at the inside corner `v` take a junction
/// there: whether their contacts on the face they share, carried on past
/// the corner, cross within the solid. Where they cross only past its far
/// side (a crease turning back nearly on itself), the ball rolls round the
/// sharp edge instead ([`rolled_corner`]).
///
/// # Errors
///
/// As reading the creases and classifying the crossing.
#[allow(
    clippy::too_many_arguments,
    reason = "one corner's reading of the solid"
)]
fn joined_at(
    model: &Model,
    solid: &Shape,
    curved: &Curved,
    creases: [&Crease; 2],
    walls: [usize; 2],
    top: usize,
    v: Id,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let at = curved.points[&v];
    let Some(nt) = outward(curved, top, at) else {
        return Ok(true);
    };
    let mut lines = Vec::with_capacity(2);
    for (c, wall) in creases.iter().zip(walls) {
        let Some(nw) = outward(curved, wall, at) else {
            return Ok(true);
        };
        let centre = -(nt + nw) * (c.radius / (1.0 + nt.dot(nw)));
        lines.push((at + centre + nt * c.radius, -c.leaving(curved, v, tol)?));
    }
    let [(p, d), (q, e)] = [lines[0], lines[1]];
    let across = d.cross(e);
    if across.magnitude() <= tol.angular() {
        return Ok(false);
    }
    let s = (q - p).cross(e).dot(across) / across.dot(across);
    let u = (q - p).cross(d).dot(across) / across.dot(across);
    if s < 0.0 || u < 0.0 {
        return Ok(true);
    }
    let crossing = p + d * s - nt * (0.1 * creases[0].radius.min(creases[1].radius));
    Ok(
        ogeom_algo::classify_in_solid_exact(model, solid, crossing, tol)?
            == ogeom_algo::Containment::In,
    )
}

/// The sharp edge at `v` other than the two creases: a straight edge, its
/// direction out of the shared face `top`.
fn sharp_edge(curved: &Curved, creases: [&Crease; 2], v: Id) -> Option<Id> {
    let mine = [id(&creases[0].shape), id(&creases[1].shape)];
    let mut found = None;
    for (e, edge) in &curved.edges {
        if mine.contains(e) || !(edge.ends.0 == v || edge.ends.1 == v) {
            continue;
        }
        if found.is_some() || !matches!(edge.curve, Curve::Line(_)) {
            return None;
        }
        found = Some(*e);
    }
    found
}

/// Whether the ball rolls round the sharp edge at `v` between the two
/// creases' bands: constant rounds of one radius, the vertex's third edge
/// straight and square to the face they share.
fn rolls_round(curved: &Curved, creases: [&Crease; 2], top: usize, v: Id, tol: Tolerances) -> bool {
    if creases.iter().any(|c| c.law.is_some())
        || (creases[0].radius - creases[1].radius).abs() > tol.confusion()
    {
        return false;
    }
    let Some(e) = sharp_edge(curved, creases, v) else {
        return false;
    };
    let edge = &curved.edges[&e];
    let at = curved.points[&v];
    let other = if edge.ends.0 == v {
        edge.ends.1
    } else {
        edge.ends.0
    };
    let d = curved.points[&other] - at;
    outward(curved, top, at)
        .is_some_and(|n| n.cross(d).magnitude() <= tol.angular() * d.magnitude())
}

/// The ball rolling round the sharp edge at `v` between the two creases'
/// bands, touching it and the face both share: the bands' section, in the
/// plane through the edge where the first band ends (its ball touching the
/// edge), turned about the edge through the angle between the two bands'
/// end planes behind the corner. Its side on the shared face stands off it,
/// a quarter radius out, over the air above it. With the turned band's
/// faces. `None` where it will not build.
fn rolled_corner(
    model: &mut Model,
    curved: &Curved,
    creases: [&Crease; 2],
    top: usize,
    v: Id,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    use ogeom_geom::{CircleCurve, LineCurve, PlaneSurface};
    let built = (|| -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
        let at = curved.points[&v];
        let Some(up) = outward(curved, top, at) else {
            return Ok(None);
        };
        // Into each band's material, square to its wall at the corner.
        let mut into = Vec::with_capacity(2);
        for c in creases {
            let Some(wall) = c.faces.iter().copied().find(|f| *f != top) else {
                return Ok(None);
            };
            let Some(nw) = outward(curved, wall, at) else {
                return Ok(None);
            };
            into.push(-nw);
        }
        let r = creases[0].radius;
        let (from, to) = (into[0], into[1]);
        let turn = from.angle(to, tol)?;
        let axis = if from.cross(to).dot(up) >= 0.0 {
            up
        } else {
            -up
        };
        let off = 0.25 * r;
        let centre = at + from * r - up * r;
        let (contact, foot) = (at + from * r, at - up * r);
        let (lift_c, lift_a) = (contact + up * off, at + up * off);
        let [vc, vf, vla, vlc] =
            [contact, foot, lift_a, lift_c].map(|q| ogeom_algo::make_vertex(model, q).shape);
        let mut line = |a: Point, b: Point, va: &Shape, vb: &Shape| -> OgeomResult<Shape> {
            Ok(ogeom_algo::make_edge_between(
                model,
                LineCurve::segment(a, b, tol)?.into(),
                (0.0, a.distance(b)),
                va,
                vb,
                tol,
            )?
            .shape)
        };
        let up_axis = line(foot, lift_a, &vf, &vla)?;
        let over = line(lift_a, lift_c, &vla, &vlc)?;
        let down = line(lift_c, contact, &vlc, &vc)?;
        let x = contact - centre;
        let y = foot - centre;
        let frame = Frame::new(
            centre,
            Direction::new(x.cross(y), tol)?,
            Direction::new(x, tol)?,
            tol,
        )?;
        let arc_curve: Curve = CircleCurve::new(ogeom_math::Circle::new(frame, r, tol)?).into();
        let arc = ogeom_algo::make_edge_between(
            model,
            arc_curve,
            (0.0, core::f64::consts::FRAC_PI_2),
            &vc,
            &vf,
            tol,
        )?
        .shape;
        let wire = ogeom_algo::make_wire(model, &[arc.clone(), up_axis, over, down], tol)?.shape;
        let plane = PlaneSurface::new(ogeom_math::Plane::through(
            at,
            Direction::new(from.cross(up), tol)?,
        ));
        let face = ogeom_algo::make_face(model, plane.into(), &[wire], tol)?.shape;
        let turned = ogeom_algo::make_revolution(
            model,
            &face,
            Axis {
                location: at,
                direction: Direction::new(axis, tol)?,
            },
            turn,
            tol,
        )?;
        let band = turned.history.generated(&arc).to_vec();
        Ok((!band.is_empty()).then_some((turned.shape, band)))
    })();
    match built {
        Ok(found) => Ok(found),
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}

/// `edge` as a crease: straight, open, between two planar faces, convex;
/// `None` otherwise.
fn crease_of(
    model: &Model,
    curved: &Curved,
    edge: Id,
    shape: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Crease>> {
    let Some(c) = curved.edges.get(&edge) else {
        return Ok(None);
    };
    if c.ends.0 == c.ends.1 {
        return Ok(None);
    }
    if matches!(c.curve, Curve::Circle(_)) {
        return rim_of(model, curved, edge, shape, radius, tol);
    }
    if !matches!(c.curve, Curve::Line(_)) {
        return Ok(None);
    }
    // The crease's ends are its vertices brought onto its curve's line: a
    // vertex standing off its edge's line within its tolerance (a corner
    // the wall's plane does not quite reach) would carry the crease's seat
    // off its faces' planes.
    let (a, b) = {
        let Curve::Line(line) = &c.curve else {
            return Ok(None);
        };
        let axis = line.axis();
        let onto = |v: Point| {
            let along = axis.direction.vector();
            axis.location + along * along.dot(v - axis.location)
        };
        (
            onto(curved.points[&c.ends.0]),
            onto(curved.points[&c.ends.1]),
        )
    };
    let [fa, fb] = faces_of(curved, edge)[..] else {
        return Ok(None);
    };
    if [fa, fb]
        .iter()
        .any(|f| !matches!(curved.faces[*f].surface, SurfaceGeometry::Plane(_)))
    {
        return Ok(None);
    }
    let mid = a.midpoint(b);
    let (Some(na), Some(nb)) = (outward(curved, fa, mid), outward(curved, fb, mid)) else {
        return Ok(None);
    };
    // Convex: each face's points away from the edge stand behind the
    // other's plane (the material between them), read at its edges' ends
    // and middles: a D-shaped block's top has no vertex off its chord, its
    // arc's middle standing behind the chord's wall.
    let behind = |f: usize, n: Vector| {
        curved.faces[f].loops.iter().flatten().any(|(x, _)| {
            let g = &curved.edges[x];
            let middle = {
                use ogeom_geom::Curve3d as _;
                g.curve
                    .point_at(f64::midpoint(g.range.0, g.range.1), tol)
                    .ok()
            };
            [g.ends.0, g.ends.1]
                .iter()
                .map(|v| curved.points[v])
                .chain(middle)
                .any(|p| (p - mid).dot(n) < -tol.confusion())
        })
    };
    if !behind(fb, na) || !behind(fa, nb) || na.cross(nb).magnitude() <= tol.angular() {
        return Ok(None);
    }
    Ok(Some(Crease {
        shape: shape.clone(),
        ends: c.ends,
        a,
        b,
        faces: [fa, fb],
        planes: [
            (na, na.dot(a - Point::ORIGIN)),
            (nb, nb.dot(a - Point::ORIGIN)),
        ],
        radius,
        kind: Kind::Straight,
        law: None,
    }))
}

/// The spandrel of `crease` carried on along its line `reach` past both
/// ends: the wedge of its two faces' planes, capped square to the crease
/// and behind it, less that wedge with its edge rounded by the engine (an
/// exact cylinder); with the band's faces in it.
fn spandrel(
    model: &mut Model,
    curved: &Curved,
    crease: &Crease,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    if let Kind::Rim { range, top } = crease.kind {
        return rim_spandrel(model, curved, crease, range, top, reach, tol);
    }
    let eps = tol.confusion();
    let along = (crease.b - crease.a) / crease.a.distance(crease.b);
    let o = crease.a;
    let o_along = along.dot(o - Point::ORIGIN);
    let bisect = crease.planes[0].0 + crease.planes[1].0;
    let bisect = bisect / bisect.magnitude();
    let mut planes = crease.planes.to_vec();
    planes.push((along, o_along + crease.a.distance(crease.b) + reach));
    planes.push((-along, -(o_along - reach)));
    // Behind the crease, a few radii in: the wedge's back.
    planes.push((
        -bisect,
        -(bisect.dot(o - Point::ORIGIN) - 4.0 * crease.radius),
    ));
    let Some(wedge) = polytope(model, &planes, tol)? else {
        return Ok(None);
    };
    // The wedge's edge where its two faces' planes meet (an edge's own
    // vertices may stand off that line by a little).
    let line_point = |p: Point| {
        crease
            .planes
            .iter()
            .all(|(n, d)| (n.dot(p - Point::ORIGIN) - d).abs() <= 10.0 * eps)
    };
    let point = |model: &Model, v: &Shape| {
        model
            .node(v)
            .and_then(|n| n.data().as_vertex().map(|d| d.point))
    };
    let mut edge = None;
    for e in explore_unique(model, &wedge, ShapeType::Edge)? {
        if let Some((p, q)) = ogeom_algo::edge_vertices(model, &e)?
            && let (Some(p), Some(q)) = (point(model, &p), point(model, &q))
            && line_point(p)
            && line_point(q)
        {
            edge = Some(e);
        }
    }
    let Some(edge) = edge else {
        return Ok(None);
    };
    // A varying round's law over the crease, held level past its ends, laid
    // over the wedge's edge from the end before the crease's start.
    let section = match &crease.law {
        None => super::Section::Round {
            radius: crease.radius,
        },
        Some(profile) => {
            let length = crease.a.distance(crease.b);
            let before = crease.a - along * reach;
            let Some((p, q)) = ogeom_algo::edge_vertices(model, &edge)? else {
                return Ok(None);
            };
            let (Some(pp), Some(pq)) = (point(model, &p), point(model, &q)) else {
                return Ok(None);
            };
            let start = if pp.distance(before) <= pq.distance(before) {
                p
            } else {
                q
            };
            super::Section::Varying {
                law: crate::law::EdgeLaw::Laid(profile.stretched(reach / length, reach / length)),
                start,
            }
        }
    };
    let request = super::Request {
        blends: vec![super::request::Blend {
            edge: edge.clone(),
            section,
        }],
        mitre: false,
        splits: Vec::new(),
        corners: Vec::new(),
    };
    let Ok(Some(rounded)) = super::run(model, &wedge, &request, tol) else {
        return Ok(None);
    };
    let band: Vec<Shape> = rounded.history.generated(&edge).to_vec();
    let Ok(cut) = ogeom_bool::cut(model, &wedge, &rounded.shape, tol) else {
        return Ok(None);
    };
    // The band's faces as the spandrel holds them: traced through the cut,
    // or where it recorded none, the faces on the band's own surface.
    let kept = explore_unique(model, &cut.shape, ShapeType::Face)?;
    let mut named: Vec<Shape> = band
        .iter()
        .flat_map(|b| cut.history.trace(b).to_vec())
        .filter(|f| kept.iter().any(|k| k.is_same(f)))
        .collect();
    if named.is_empty() {
        let surface = |s: &Shape| {
            model
                .node(s)
                .and_then(|n| n.data().as_face().map(|d| d.surface))
        };
        named = kept
            .into_iter()
            .filter(|f| {
                band.iter()
                    .any(|b| surface(b).is_some() && surface(b) == surface(f))
            })
            .collect();
    }
    Ok(Some((cut.shape, named)))
}

/// `edge` as a rim crease: an arc of a circle between a plane square to
/// its axis and a coaxial cylinder or cone, convex; `None` otherwise.
fn rim_of(
    model: &Model,
    curved: &Curved,
    edge: Id,
    shape: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Crease>> {
    let ang = tol.angular();
    let c = &curved.edges[&edge];
    let Curve::Circle(circle) = &c.curve else {
        return Ok(None);
    };
    let axis = circle.circle().frame().z().vector();
    let [fa, fb] = faces_of(curved, edge)[..] else {
        return Ok(None);
    };
    let kind = |f: usize| match curved.faces[f].surface {
        SurfaceGeometry::Plane(_) => Some(true),
        SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_) => Some(false),
        _ => None,
    };
    let (top, wall) = match (kind(fa), kind(fb)) {
        (Some(true), Some(false)) => (fa, fb),
        (Some(false), Some(true)) => (fb, fa),
        _ => return Ok(None),
    };
    let (a, b) = (curved.points[&c.ends.0], curved.points[&c.ends.1]);
    let t_mid = 0.5 * (c.range.0 + c.range.1);
    let mid = c.curve.point_at(t_mid, tol)?;
    let (Some(nt), Some(nw)) = (outward(curved, top, mid), outward(curved, wall, mid)) else {
        return Ok(None);
    };
    if nt.cross(axis).magnitude() > ang {
        return Ok(None);
    }
    // Convex: the top runs from the rim behind the wall's tangent plane.
    let tangent = c.curve.d1_at(t_mid, tol)?;
    let length = c.curve.point_at(c.range.0, tol)?.distance(mid) * 2.0;
    let Some(into) = super::far_edge::into_face(
        model,
        curved,
        top,
        (mid, tangent / tangent.magnitude(), length),
        tol,
    )?
    else {
        return Ok(None);
    };
    if into.dot(nw) >= -ang {
        return Ok(None);
    }
    let [ta, wa] = [top, wall].map(|f| outward(curved, f, a).unwrap_or(nt));
    let planes = if top == fa {
        [
            (ta, ta.dot(a - Point::ORIGIN)),
            (wa, wa.dot(a - Point::ORIGIN)),
        ]
    } else {
        [
            (wa, wa.dot(a - Point::ORIGIN)),
            (ta, ta.dot(a - Point::ORIGIN)),
        ]
    };
    Ok(Some(Crease {
        shape: shape.clone(),
        ends: c.ends,
        a,
        b,
        faces: [fa, fb],
        planes,
        radius,
        kind: Kind::Rim {
            range: c.range,
            top,
        },
        law: None,
    }))
}

/// A rim's section at parameter `t` of its circle: the face in the plane
/// through the axis there, bounded by the rim's point, its two contacts
/// (a radius in from both faces) and the ball's arc between them; with
/// that arc's edge and the circle's unit tangent there.
fn rim_section(
    model: &mut Model,
    curved: &Curved,
    crease: &Crease,
    top: usize,
    t: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Shape, Vector)>> {
    use ogeom_geom::{CircleCurve, LineCurve, PlaneSurface};
    let r = crease.radius;
    let curve = &curved.edges[&id(&crease.shape)].curve;
    let wall = if crease.faces[0] == top {
        crease.faces[1]
    } else {
        crease.faces[0]
    };
    let p = curve.point_at(t, tol)?;
    let tangent = curve.d1_at(t, tol)?;
    let tangent = tangent / tangent.magnitude();
    let (Some(nt), Some(nw)) = (outward(curved, top, p), outward(curved, wall, p)) else {
        return Ok(None);
    };
    let g = nt.dot(nw);
    if (1.0 + g).abs() <= tol.angular() {
        return Ok(None);
    }
    let centre = p - (nt + nw) * (r / (1.0 + g));
    let (on_top, on_wall) = (centre + nt * r, centre + nw * r);
    let (x, y) = (on_top - centre, on_wall - centre);
    let normal = x.cross(y);
    if normal.magnitude() <= tol.confusion() * r {
        return Ok(None);
    }
    // The section's other sides stand off the wall's surface: from the top
    // contact along the top past the rim, out through the open side, and in
    // to the wall contact square to the wall. A tool's face on the wall's
    // own cone would meet the solid's along it, which the boolean does not
    // resolve; past the wall the tool takes only the open side.
    let off = r * 0.25;
    let past = p + (p - on_top) / p.distance(on_top) * off;
    let out = on_wall + nw * off;
    let [vt, vq, vo, vw] =
        [on_top, past, out, on_wall].map(|q| ogeom_algo::make_vertex(model, q).shape);
    let mut line = |from: Point, to: Point, a: &Shape, b: &Shape| -> OgeomResult<Shape> {
        Ok(ogeom_algo::make_edge_between(
            model,
            LineCurve::segment(from, to, tol)?.into(),
            (0.0, from.distance(to)),
            a,
            b,
            tol,
        )?
        .shape)
    };
    let back = line(on_wall, out, &vw, &vo)?;
    let across = line(out, past, &vo, &vq)?;
    let along_top = line(past, on_top, &vq, &vt)?;
    let frame = Frame::new(
        centre,
        Direction::new(normal, tol)?,
        Direction::new(x, tol)?,
        tol,
    )?;
    let arc_curve: Curve = CircleCurve::new(ogeom_math::Circle::new(frame, r, tol)?).into();
    let arc =
        ogeom_algo::make_edge_between(model, arc_curve, (0.0, x.angle(y, tol)?), &vt, &vw, tol)?
            .shape;
    let wire = ogeom_algo::make_wire(model, &[arc.clone(), back, across, along_top], tol)?.shape;
    let plane = PlaneSurface::new(ogeom_math::Plane::through(p, Direction::new(tangent, tol)?));
    let face = ogeom_algo::make_face(model, plane.into(), &[wire], tol)?.shape;
    Ok(Some((face, arc, tangent)))
}

/// A rim's spandrel carried on: its section turned about the axis over the
/// rim's arc (the torus band's spandrel), and past each end the end
/// section run straight along the tangent `reach` on; with the band's
/// faces in it.
fn rim_spandrel(
    model: &mut Model,
    curved: &Curved,
    crease: &Crease,
    range: (f64, f64),
    top: usize,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let curve = &curved.edges[&id(&crease.shape)].curve;
    let Curve::Circle(circle) = curve else {
        return Ok(None);
    };
    let circle = circle.circle();
    let built = (|| -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
        let (Some((f0, arc0, t0)), Some((f1, arc1, t1))) = (
            rim_section(model, curved, crease, top, range.0, tol)?,
            rim_section(model, curved, crease, top, range.1, tol)?,
        ) else {
            return Ok(None);
        };
        let turned = ogeom_algo::make_revolution(
            model,
            &f0,
            Axis {
                location: circle.centre(),
                direction: circle.frame().z(),
            },
            range.1 - range.0,
            tol,
        )?;
        let before = ogeom_algo::make_prism(model, &f0, -t0 * reach, tol)?;
        let after = ogeom_algo::make_prism(model, &f1, t1 * reach, tol)?;
        let mut band: Vec<Shape> = turned.history.generated(&arc0).to_vec();
        band.extend(before.history.generated(&arc0).iter().cloned());
        band.extend(after.history.generated(&arc1).iter().cloned());
        let one = ogeom_bool::fuse(model, &turned.shape, &before.shape, tol)?;
        band = band
            .iter()
            .flat_map(|f| one.history.trace(f).to_vec())
            .collect();
        let all = ogeom_bool::fuse(model, &one.shape, &after.shape, tol)?;
        band = band
            .iter()
            .flat_map(|f| all.history.trace(f).to_vec())
            .collect();
        let kept = explore_unique(model, &all.shape, ShapeType::Face)?;
        band.retain(|f| kept.iter().any(|k| k.is_same(f)));
        Ok((!band.is_empty()).then_some((all.shape, band)))
    })();
    match built {
        Ok(found) => Ok(found),
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}

/// The open side of a rim's curved wall within `bounds`: outside the wall's
/// cylinder or cone where the material stands inside it (a drum), inside it
/// where the material stands outside (a bore); the frustum taken a few
/// radii either way of the rim. `None` where the cone would close within
/// that, or a solid will not build.
fn past_wall(
    model: &mut Model,
    curved: &Curved,
    crease: &Crease,
    wall: usize,
    bounds: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let eps = tol.confusion();
    let Curve::Circle(circle) = &curved.edges[&id(&crease.shape)].curve else {
        return Ok(None);
    };
    let circle = circle.circle();
    let (o, z) = (circle.centre(), circle.frame().z().vector());
    let p = crease.a;
    let radial = (p - o) - z * (p - o).dot(z);
    let rho = radial.magnitude();
    let er = radial / rho;
    let Some(nw) = outward(curved, wall, p) else {
        return Ok(None);
    };
    let tangent = z.cross(er);
    let mut g = nw.cross(tangent);
    if g.dot(z) < 0.0 {
        g = -g;
    }
    if g.dot(z).abs() <= tol.angular() {
        return Ok(None);
    }
    let slope = g.dot(er) / g.dot(z);
    let h = 6.0 * crease.radius;
    let (lo, hi) = (rho - h * slope, rho + h * slope);
    if lo <= eps || hi <= eps {
        return Ok(None);
    }
    let frame = Frame::new(
        o - z * h,
        Direction::new(z, tol)?,
        Direction::new(er, tol)?,
        tol,
    )?;
    let built = (|| -> OgeomResult<Shape> {
        let frustum = match curved.faces[wall].surface {
            SurfaceGeometry::Cylinder(_) => {
                ogeom_algo::make_cylinder(model, frame, rho, 2.0 * h, tol)?
            }
            _ => ogeom_algo::make_cone(model, frame, lo, hi, 2.0 * h, tol)?,
        }
        .shape;
        if nw.dot(er) > 0.0 {
            Ok(ogeom_bool::cut(model, bounds, &frustum, tol)?.shape)
        } else {
            Ok(frustum)
        }
    })();
    match built {
        Ok(shape) => Ok(Some(shape)),
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}
