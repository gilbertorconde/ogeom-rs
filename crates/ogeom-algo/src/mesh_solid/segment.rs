//! Triangles gathered into regions: planar faces grown across shared
//! edges, the segmentation into recognized curved regions and planar faces
//! over the rest, and the chart of a canonical surface (a point's
//! coordinates, the point at some, which directions wrap, the branch round
//! a region's centre) that the later stages read curved regions through.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::Curve;
use ogeom_math::{Cone, Cylinder, Direction, Frame, Plane, Point, Sphere, Torus, Vector};

use super::frames::{align_axes, hole_frames, slit_bands};
use super::planner::EdgeSpec;
use super::regions::{
    merge_same_surface, patch_regions, recognized_regions, sphere_axes, split_disconnected,
    swept_regions,
};
use super::rounds::{faceted_rounds, plane_normals, tangent_blends, tangent_rounds};
use super::snap::plane_through;
use super::weld::{Adjacency, Half};
use super::{Carrier, Curved, Groups, MeshSolidOptions, Narrower, PatchRefusals};
use crate::recognize::Canonical;

/// The plane of a triangle, its normal by the winding, its `x` axis along
/// the first side. The triangle has area, so both are defined.
pub(super) fn plane_of(
    points: &[Point],
    [a, b, c]: [u32; 3],
    tol: Tolerances,
) -> OgeomResult<Plane> {
    let [a, b, c] = [a, b, c].map(|i| points[i as usize]);
    // Scaled to unit length before the kernel's own normalization, which
    // refuses a vector shorter than the confusion distance, as the cross
    // product of a small triangle's sides may be.
    let normal = (b - a).cross(c - a);
    let z = Direction::new(normal / normal.magnitude(), tol)?;
    let x = Direction::new((b - a) / (b - a).magnitude(), tol)?;
    Ok(Plane::new(Frame::new(a, z, x, tol)?))
}

pub(super) fn unit_normal(points: &[Point], [a, b, c]: [u32; 3]) -> Vector {
    let [a, b, c] = [a, b, c].map(|i| points[i as usize]);
    let n = (b - a).cross(c - a);
    n / n.magnitude()
}

/// Grow planar faces across shared edges from the largest triangles down,
/// among the triangles no face holds yet, taking a neighbour whose normal
/// is within the angle of the face's and whose corners all lie within the
/// distance of its plane. Measured against the face's own plane, not the
/// last triangle's, so a gently curved surface does not drift into one
/// face: the seed's plane at first, then the plane fitted to what the face
/// has gathered.
pub(super) fn coplanar_groups(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    options: &MeshSolidOptions,
    flat: f64,
    groups: &mut Groups,
    tol: Tolerances,
) -> OgeomResult<()> {
    let area = |t: usize| {
        let [a, b, c] = triangles[t].map(|i| points[i as usize]);
        (b - a).cross(c - a).magnitude()
    };
    let mut order: Vec<usize> = (0..triangles.len())
        .filter(|&t| groups.of[t] == usize::MAX)
        .collect();
    order.sort_by(|&x, &y| area(y).total_cmp(&area(x)));
    let cos = options.coplanar_angle.cos();
    let quantum = options.quantum.unwrap_or(0.0);
    let mut stack = Vec::new();
    for seed in order {
        if groups.of[seed] != usize::MAX {
            continue;
        }
        let g = groups.carriers.len();
        let mut plane = plane_of(points, triangles[seed], tol)?;
        groups.of[seed] = g;
        let mut members = vec![seed];
        let mut vertices: Vec<u32> = triangles[seed].to_vec();
        let mut rim = vec![seed];
        loop {
            let (origin, normal) = (plane.frame().origin(), plane.frame().z().vector());
            stack.clone_from(&rim);
            rim.clear();
            let grown = members.len();
            while let Some(t) = stack.pop() {
                for h in 3 * t..3 * t + 3 {
                    let Some(twin) = adjacency.twin[h] else {
                        continue;
                    };
                    let other = twin / 3;
                    if groups.of[other] != usize::MAX {
                        continue;
                    }
                    let [a, b, c] = triangles[other].map(|i| points[i as usize]);
                    let n = (b - a).cross(c - a);
                    // A small triangle's corners, each rounded by the
                    // mesh's encoding, tilt its normal by up to twice the
                    // rounding over its smallest altitude: on a mesh stored
                    // far from its origin, a tilt that is no turn of the
                    // surface. Its corners' distance decides it then.
                    let longest = (b - a)
                        .magnitude()
                        .max((c - b).magnitude())
                        .max((a - c).magnitude());
                    let altitude = n.magnitude() / longest.max(f64::MIN_POSITIVE);
                    let slop = (2.0 * quantum / altitude.max(f64::MIN_POSITIVE))
                        .min(1.0)
                        .asin();
                    // A triangle no higher than the distance (three corners
                    // all but on a line) has no plane of its own at that
                    // distance: built alone it is a face of almost no area,
                    // turned whichever way its rounding points.
                    let cos = if altitude <= flat {
                        -1.0
                    } else {
                        cos.min(slop.cos())
                    };
                    if n.dot(normal) < cos * n.magnitude()
                        || [a, b, c]
                            .iter()
                            .any(|p| (*p - origin).dot(normal).abs() > flat)
                    {
                        rim.push(t);
                        continue;
                    }
                    groups.of[other] = g;
                    members.push(other);
                    vertices.extend(triangles[other]);
                    stack.push(other);
                }
            }
            // A plane through one triangle leans with its corners' slop, and
            // across a large face that lean carries the far side past the
            // distance. Refitted to everything gathered, where the fit still
            // holds it all, the face grows again from its rim.
            if members.len() == grown || rim.is_empty() {
                break;
            }
            vertices.sort_unstable();
            vertices.dedup();
            let at: Vec<Point> = vertices.iter().map(|&v| points[v as usize]).collect();
            let Some((centre, fitted)) = plane_through(&at, tol) else {
                break;
            };
            let fitted = if fitted.vector().dot(normal) < 0.0 {
                fitted.reversed()
            } else {
                fitted
            };
            if at
                .iter()
                .any(|p| (*p - centre).dot(fitted.vector()).abs() > flat)
            {
                break;
            }
            let x = plane.frame().x().vector();
            let x = x - fitted.vector() * x.dot(fitted.vector());
            plane = Plane::new(Frame::new(centre, fitted, Direction::new(x, tol)?, tol)?);
            rim.sort_unstable();
            rim.dedup();
        }
        rim.clear();
        groups.carriers.push(Carrier::Plane(plane));
    }
    Ok(())
}

/// One face per triangle, for the triangles no face holds yet.
pub(super) fn one_each(
    points: &[Point],
    triangles: &[[u32; 3]],
    groups: &mut Groups,
    tol: Tolerances,
) -> OgeomResult<()> {
    for (t, triangle) in triangles.iter().enumerate() {
        if groups.of[t] == usize::MAX {
            groups.of[t] = groups.carriers.len();
            groups
                .carriers
                .push(Carrier::Plane(plane_of(points, *triangle, tol)?));
        }
    }
    Ok(())
}

/// The direction of a canonical surface's own normal at a point near it,
/// not normalized: the gradient of its distance.
/// How far a sample's normals must turn before the smallest sample stage
/// is fitted: twenty degrees.
pub(super) const SMALL_STAGE_TURN: f64 = 0.35;

/// The widest angle between any two of a sample's normals.
pub(super) fn turn_of(normals: &[Vector]) -> f64 {
    let mut widest: f64 = 0.0;
    for (i, a) in normals.iter().enumerate() {
        for b in &normals[i + 1..] {
            widest = widest.max(a.dot(*b).clamp(-1.0, 1.0).acos());
        }
    }
    widest
}

/// Whether a facet whose corners lie on `shape` leans only as the surface
/// turns under it: its normal within the spread of the surface's normals
/// at its corners (and [`FACET_LEAN`]) of their mean, the surface turning
/// no more than sixty degrees under it. A facet of a coarse mesh over a
/// tight bend passes however far it turns from its neighbours; a flat
/// cap's triangle with its corners on a cylinder's rim, square to every
/// one of those normals, does not.
pub(super) fn leans_as_the_surface(shape: &Canonical, corners: [Point; 3], normal: Vector) -> bool {
    let mut at = [Vector::ZERO; 3];
    for (n, p) in at.iter_mut().zip(corners) {
        let g = gradient(shape, p);
        let m = g.magnitude();
        if m == 0.0 {
            return false;
        }
        *n = if g.dot(normal) < 0.0 { -g / m } else { g / m };
    }
    let angle = |a: Vector, b: Vector| a.dot(b).clamp(-1.0, 1.0).acos();
    let spread = angle(at[0], at[1])
        .max(angle(at[1], at[2]))
        .max(angle(at[0], at[2]));
    // A facet of any mesh turns far less than this across itself; one that
    // spans more is a chord across the surface, not a piece of it.
    if spread > FACET_TURN {
        return false;
    }
    let mean = at[0] + at[1] + at[2];
    let m = mean.magnitude();
    m > 0.0 && angle(mean / m, normal) <= spread + FACET_LEAN
}

/// How far a facet's normal may lean past the spread of the surface's
/// normals at its corners: a skewed triangle across a cylinder, joining
/// points at different heights and angles, tilts along the axis out of the
/// plane its corners' normals span. About eleven degrees.
const FACET_LEAN: f64 = 0.2;

/// The most a surface may turn under one facet: sixty degrees.
const FACET_TURN: f64 = core::f64::consts::FRAC_PI_3;

pub(super) fn gradient(shape: &Canonical, p: Point) -> Vector {
    let radial = |o: Point, z: Vector| {
        let w = p - o;
        let r = w - z * w.dot(z);
        let m = r.magnitude();
        (if m > 0.0 { r / m } else { Vector::ZERO }, w.dot(z))
    };
    match shape {
        Canonical::Plane(plane) => plane.frame().z().vector(),
        Canonical::Cylinder(c) => radial(c.frame().origin(), c.frame().z().vector()).0,
        Canonical::Cone(c) => {
            let z = c.frame().z().vector();
            let (out, _) = radial(c.frame().origin(), z);
            out - z * c.half_angle().tan()
        }
        Canonical::Sphere(s) => p - s.centre(),
        Canonical::Torus(t) => {
            let z = t.frame().z().vector();
            let (out, _) = radial(t.frame().origin(), z);
            p - (t.frame().origin() + out * t.major_radius())
        }
        Canonical::Swept(s) => {
            use ogeom_geom::Surface as _;
            s.foot(p)
                .and_then(|f| {
                    s.surface
                        .normal_at(f.parameters.0, f.parameters.1, s.tol)
                        .ok()
                })
                .map_or(Vector::ZERO, |n| n.vector())
        }
    }
}

/// The axis of a surface of revolution: its frame.
pub(super) fn axis_frame(shape: &Canonical) -> Option<Frame> {
    match shape {
        Canonical::Cylinder(c) => Some(c.frame()),
        Canonical::Cone(c) => Some(c.frame()),
        Canonical::Torus(t) => Some(t.frame()),
        Canonical::Sphere(s) => Some(s.frame()),
        Canonical::Plane(_) | Canonical::Swept(_) => None,
    }
}

/// The cone about the same axis, turned about it so its seam (where its
/// angle is zero) runs through `at`; `None` where `at` is on the axis.
pub(super) fn cone_seamed_through(cone: &Cone, at: Point, tol: Tolerances) -> Option<Cone> {
    let old = cone.frame();
    let x = Direction::new(at - old.origin(), tol).ok()?;
    let mut frame = Frame::new(old.origin(), old.z(), x, tol).ok()?;
    if frame.handedness() != old.handedness() {
        frame = frame.mirrored();
    }
    Cone::new(frame, cone.reference_radius(), cone.half_angle(), tol).ok()
}

/// A ring's planned edges in walking order, repeats run together, as
/// `entry` names each half-edge's.
pub(super) fn walk_entries(
    ring: &[Half],
    entry: impl Fn(Half) -> (usize, bool),
) -> Vec<(usize, bool)> {
    let mut entries: Vec<(usize, bool)> = Vec::new();
    for &h in ring {
        let e = entry(h);
        if entries.last() != Some(&e) {
            entries.push(e);
        }
    }
    if entries.len() > 1 && entries.first() == entries.last() {
        entries.pop();
    }
    entries
}

/// Where a cone's rim starts, walked as `entries` run: the rim one full
/// circle (its own start), or arcs of one parallel of the cone, each a
/// circle about the cone's axis through the same height within `reach`,
/// together going once round (the first one's start as walked). `None`
/// for any other rim.
pub(super) fn cone_rim_start(
    cone: &Cone,
    specs: &[EdgeSpec],
    entries: &[(usize, bool)],
    reach: f64,
    tol: Tolerances,
) -> Option<Point> {
    use ogeom_geom::Curve3d as _;
    if entries.iter().any(|&(e, _)| e == usize::MAX) {
        return None;
    }
    if let [(edge, _)] = entries {
        let spec = &specs[*edge];
        return if spec.closed_circle {
            spec.curve.point_at(0.0, tol).ok()
        } else {
            None
        };
    }
    let frame = cone.frame();
    let (o, z) = (frame.origin(), frame.z().vector());
    let mut height: Option<f64> = None;
    let (mut sweep, mut slack) = (0.0, 0.0);
    for &(edge, _) in entries {
        let spec = &specs[edge];
        let Curve::Circle(c) = &spec.curve else {
            return None;
        };
        let circle = c.circle();
        let w = circle.frame().origin() - o;
        let along = w.dot(z);
        if circle.frame().z().vector().cross(z).magnitude() > tol.angular()
            || (w - z * along).magnitude() > reach
            || (cone.radius_at(along) - circle.radius()).abs() > reach
            || height.is_some_and(|h| (h - along).abs() > reach)
        {
            return None;
        }
        height = Some(along);
        sweep += (spec.range.1 - spec.range.0).abs();
        slack = reach / circle.radius();
    }
    if (sweep - core::f64::consts::TAU).abs() > slack {
        return None;
    }
    let (edge, forward) = entries[0];
    let spec = &specs[edge];
    spec.curve
        .point_at(if forward { spec.range.0 } else { spec.range.1 }, tol)
        .ok()
}

/// The same surface on another frame whose axis is the same line: its
/// radii kept, a cone's reference radius carried to the new origin.
pub(super) fn on_frame(shape: &Canonical, frame: Frame, tol: Tolerances) -> Option<Canonical> {
    Some(match shape {
        Canonical::Cylinder(c) => Canonical::Cylinder(Cylinder::new(frame, c.radius(), tol).ok()?),
        Canonical::Cone(c) => {
            // The new origin's height on the old axis, and which way the
            // new axis runs against the old.
            let old = c.frame();
            let shift = (frame.origin() - old.origin()).dot(old.z().vector());
            let same = frame.z().vector().dot(old.z().vector()) > 0.0;
            if !same {
                return None;
            }
            let r0 = c.radius_at(shift);
            Canonical::Cone(Cone::new(frame, r0.max(tol.confusion()), c.half_angle(), tol).ok()?)
        }
        Canonical::Torus(t) => {
            Canonical::Torus(Torus::new(frame, t.major_radius(), t.minor_radius(), tol).ok()?)
        }
        Canonical::Sphere(s) => Canonical::Sphere(Sphere::new(frame, s.radius(), tol).ok()?),
        Canonical::Plane(_) | Canonical::Swept(_) => return None,
    })
}

/// A point's raw chart coordinates on a canonical surface.
pub(super) fn chart(shape: &Canonical, p: Point, tol: Tolerances) -> Option<(f64, f64)> {
    use ogeom_math::elementary as e;
    match shape {
        Canonical::Plane(plane) => Some(e::plane_parameters(plane, p)),
        Canonical::Cylinder(c) => e::cylinder_parameters(c, p, tol).ok(),
        Canonical::Cone(c) => e::cone_parameters(c, p, tol).ok(),
        Canonical::Sphere(s) => e::sphere_parameters(s, p, tol).ok(),
        Canonical::Torus(t) => e::torus_parameters(t, p, tol).ok(),
        Canonical::Swept(s) => s.foot(p).map(|f| f.parameters),
    }
}

pub(super) fn evaluate(shape: &Canonical, (u, v): (f64, f64)) -> Point {
    use ogeom_math::elementary as e;
    match shape {
        Canonical::Plane(plane) => e::plane_at(plane, u, v).point,
        Canonical::Cylinder(c) => e::cylinder_at(c, u, v).point,
        Canonical::Cone(c) => e::cone_at(c, u, v).point,
        Canonical::Sphere(s) => e::sphere_at(s, u, v).point,
        Canonical::Torus(t) => e::torus_at(t, u, v).point,
        Canonical::Swept(s) => {
            use ogeom_geom::Surface as _;
            s.surface.point_at(u, v, s.tol).unwrap_or(Point::ORIGIN)
        }
    }
}

/// Which chart directions are angles, and so wrap.
pub(super) fn periodic(shape: &Canonical) -> (bool, bool) {
    match shape {
        Canonical::Plane(_) => (false, false),
        Canonical::Cylinder(_) | Canonical::Cone(_) | Canonical::Sphere(_) => (true, false),
        Canonical::Torus(_) => (true, true),
        Canonical::Swept(s) => {
            use ogeom_geom::Surface as _;
            (s.surface.is_periodic_u(), s.surface.is_periodic_v())
        }
    }
}

/// Chart coordinates on the branch around `centre`.
pub(super) fn unwrapped(curved: &Curved, p: Point, tol: Tolerances) -> Option<(f64, f64)> {
    let (u, v) = chart(&curved.shape, p, tol)?;
    let (pu, pv) = periodic(&curved.shape);
    let near = |x: f64, c: f64, wraps: bool| {
        if wraps {
            c + ogeom_math::elementary::wrap_signed_angle(x - c)
        } else {
            x
        }
    };
    Some((near(u, curved.centre.0, pu), near(v, curved.centre.1, pv)))
}

/// The circular mean of angles, and the widest gap between them.
pub(super) fn angular_spread(angles: &mut [f64]) -> (f64, f64) {
    let (s, c) = angles
        .iter()
        .fold((0.0, 0.0), |(s, c), a| (s + a.sin(), c + a.cos()));
    angles.sort_by(f64::total_cmp);
    let mut gap: f64 = 0.0;
    for w in angles.windows(2) {
        gap = gap.max(w[1] - w[0]);
    }
    if let (Some(first), Some(last)) = (angles.first(), angles.last()) {
        gap = gap.max(first + core::f64::consts::TAU - last);
    }
    (s.atan2(c), gap)
}

/// Grow regions of triangles that recognition says lie on one curved
/// canonical surface, then planar faces over the rest.
///
/// A region starts at a triangle with a curved edge (one across which the
/// surface turns by less than the crease angle but more than the coplanar
/// angle) and first grows across such edges only, which keeps a flat face
/// tangent to a fillet out of the fillet's first samples. Once it holds
/// enough vertices it is recognized, and from then on grows across any
/// smooth edge to a triangle whose corners lie on the surface and whose
/// normal agrees with it, the surface refitted to everything held as the
/// region doubles. A region whose samples are free-form, or also flat,
/// returns its triangles to the planar pass.
#[allow(clippy::too_many_arguments, reason = "the segmentation's inputs")]
pub(super) fn segment(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    options: &MeshSolidOptions,
    flat: f64,
    tol: Tolerances,
) -> OgeomResult<Groups> {
    let n = triangles.len();
    let mut groups = Groups {
        of: vec![usize::MAX; n],
        carriers: Vec::new(),
        refused: PatchRefusals::default(),
        fans: std::collections::BTreeSet::new(),
        narrower: ogeom_core::FastMap::default(),
    };
    if !options.merge_coplanar {
        one_each(points, triangles, &mut groups, tol)?;
        return Ok(groups);
    }
    if options.recognize {
        recognized_regions(
            points,
            triangles,
            adjacency,
            options,
            flat,
            &mut groups,
            tol,
        );
        if options.sweeps {
            swept_regions(
                points,
                triangles,
                adjacency,
                options,
                flat,
                &mut groups,
                tol,
            );
        }
        if options.patches {
            patch_regions(
                points,
                triangles,
                adjacency,
                options,
                flat,
                &mut groups,
                tol,
            )?;
        }
        split_disconnected(triangles, adjacency, &mut groups);
        merge_same_surface(points, triangles, adjacency, &mut groups, flat);
        sphere_axes(points, triangles, adjacency, &mut groups, flat, tol);
        let mut planes = groups.clone();
        coplanar_groups(
            points,
            triangles,
            adjacency,
            options,
            flat,
            &mut planes,
            tol,
        )?;
        if faceted_rounds(
            points,
            triangles,
            adjacency,
            &mut groups,
            &planes,
            options.crease.cos(),
            flat,
            tol,
        ) {
            planes = groups.clone();
            coplanar_groups(
                points,
                triangles,
                adjacency,
                options,
                flat,
                &mut planes,
                tol,
            )?;
        }
        tangent_rounds(
            points,
            triangles,
            adjacency,
            &mut groups,
            &planes,
            options.crease.cos(),
            flat,
            None,
            tol,
        );
        let normals = plane_normals(&planes);
        align_axes(points, &mut groups, &normals, flat, tol);
        // Pieces of one torus met as spheres are put on it apart, and
        // are one region.
        if tangent_blends(
            points,
            triangles,
            adjacency,
            &mut groups,
            &planes,
            options.crease.cos(),
            flat,
            None,
            tol,
        ) {
            merge_same_surface(points, triangles, adjacency, &mut groups, flat);
        }
        hole_frames(points, triangles, adjacency, &mut groups, tol);
        slit_bands(points, triangles, adjacency, &mut groups, tol);
    }
    coplanar_groups(
        points,
        triangles,
        adjacency,
        options,
        flat,
        &mut groups,
        tol,
    )?;
    Ok(groups)
}

/// Give each curved region of `groups`, found at the wider `distance`, the
/// curved regions of `narrow`, found at the default, that hold any of its
/// triangles (none where the default left them all planar), for the build
/// to put back in its place.
pub(super) fn found_within(groups: &mut Groups, narrow: &Groups, distance: f64) {
    let mut within: ogeom_core::FastMap<usize, Vec<usize>> = ogeom_core::FastMap::default();
    for (&w, &n) in groups.of.iter().zip(&narrow.of) {
        if !matches!(groups.carriers.get(w), Some(Carrier::Curved(_))) {
            continue;
        }
        let found = within.entry(w).or_default();
        if matches!(narrow.carriers.get(n), Some(Carrier::Curved(_))) && !found.contains(&n) {
            found.push(n);
        }
    }
    for (w, found) in within {
        let regions = found
            .into_iter()
            .map(|n| {
                let held = (0..narrow.of.len())
                    .filter(|&t| narrow.of[t] == n)
                    .collect();
                (narrow.carriers[n].clone(), held)
            })
            .collect();
        groups.narrower.insert(w, Narrower { distance, regions });
    }
}
