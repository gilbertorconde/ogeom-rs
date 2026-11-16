//! Closed tools for an isolated vertex: exact corner surfaces and their caps.

use ogeom_algo::Built;
use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{CircleCurve, PlaneSurface};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape};

use super::curved::{SheetEdge, SheetFace, build_faces};
use crate::round::section::{PlanarSection, RoundSpec};

/// The portion of a convex polytope outside its tangent ball. Each wall
/// either touches the ball or passes through its centre. On the latter,
/// the two radial sides meeting at the centre become a great-circle arc;
/// those arcs bound the same spherical polygon P3 uses for edge rounds.
pub(super) fn ball_polytope(
    model: &mut Model,
    walls: &[(Vector, Point)],
    centre: Point,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let block = super::vertex_geometry::convex_block(model, walls, tol)?;
    let Some(planar) = super::classify::read_planar(model, &block, tol)? else {
        ogeom_bail!(
            Construction,
            "a corner compartment is not a planar polytope"
        );
    };
    let mut keys: Vec<_> = planar.points.keys().copied().collect();
    keys.sort_unstable();
    let points: Vec<Point> = keys.iter().map(|k| planar.points[k]).collect();
    let indices: FastMap<_, _> = keys.iter().enumerate().map(|(i, k)| (*k, i)).collect();
    let Some(pole) = points
        .iter()
        .position(|p| p.distance(centre) <= tol.confusion() * 100.0)
    else {
        ogeom_bail!(
            Construction,
            "a corner compartment has no vertex at its ball's centre"
        );
    };
    let mut edges = Vec::<SheetEdge>::new();
    let mut segments = FastMap::<(usize, usize), usize>::default();
    let mut faces = Vec::new();
    let mut arcs = Vec::new();
    for face in &planar.faces {
        let mut ring: Vec<usize> = face.corners.iter().map(|v| indices[v]).collect();
        let cap = ring.iter().position(|v| *v == pole);
        if let Some(k) = cap {
            ring.rotate_left(k);
            ring.remove(0);
        }
        let mut boundary = Vec::new();
        for k in 0..ring.len() {
            let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
            if cap.is_some() && k + 1 == ring.len() {
                let (x, y) = (points[a] - centre, points[b] - centre);
                if (x.magnitude() - radius).abs() > tol.confusion() * 100.0
                    || (y.magnitude() - radius).abs() > tol.confusion() * 100.0
                {
                    ogeom_bail!(
                        Construction,
                        "a corner compartment's cap misses the ball's contacts"
                    );
                }
                let frame = Frame::new(
                    centre,
                    Direction::new(x.cross(y), tol)?,
                    Direction::new(x, tol)?,
                    tol,
                )?;
                let angle = x.cross(y).magnitude().atan2(x.dot(y));
                let curve = CircleCurve::new(Circle::new(frame, radius, tol)?);
                let edge = edges.len();
                edges.push(SheetEdge {
                    ends: (a, b),
                    curve: Some((curve.into(), (0.0, angle))),
                });
                boundary.push((edge, true));
                arcs.push((edge, false));
            } else {
                let key = (a.min(b), a.max(b));
                let edge = *segments.entry(key).or_insert_with(|| {
                    let index = edges.len();
                    edges.push(SheetEdge {
                        ends: key,
                        curve: None,
                    });
                    index
                });
                boundary.push((edge, a < b));
            }
        }
        let frame = Frame::new(
            points[ring[0]],
            Direction::new(face.normal, tol)?,
            Direction::new(points[ring[1]] - points[ring[0]], tol)?,
            tol,
        )?;
        faces.push(SheetFace {
            surface: Some(PlaneSurface::new(Plane::new(frame)).into()),
            ring: boundary,
            holes: Vec::new(),
        });
    }
    let mut ring = Vec::new();
    let Some(first) = arcs.pop() else {
        ogeom_bail!(
            Construction,
            "a corner compartment has no spherical opening"
        );
    };
    let ends = |(e, forward): (usize, bool)| {
        let (a, b) = edges[e].ends;
        if forward { (a, b) } else { (b, a) }
    };
    ring.push(first);
    while !arcs.is_empty() {
        let end = ends(ring[ring.len() - 1]).1;
        let Some(k) = arcs.iter().position(|side| ends(*side).0 == end) else {
            ogeom_bail!(Construction, "a corner's spherical arcs do not close");
        };
        ring.push(arcs.remove(k));
    }
    if ends(ring[ring.len() - 1]).1 != ends(first).0 {
        ogeom_bail!(Construction, "a corner's spherical arcs do not close");
    }
    let feet: Vec<Point> = ring.iter().map(|side| points[ends(*side).0]).collect();
    let surface = crate::round::sphere::corner_sphere(centre, radius, &feet, tol)?;
    ring.reverse();
    for side in &mut ring {
        side.1 = !side.1;
    }
    faces.push(SheetFace {
        surface: Some(surface),
        ring,
        holes: Vec::new(),
    });
    let mut built = build_faces(model, &points, &edges, &faces, tol)?;
    let Some(sphere) = built.pop() else {
        ogeom_bail!(Construction, "a corner compartment has no spherical face");
    };
    built.push(sphere.reversed());
    let shell = ogeom_algo::make_shell(model, &built)?.shape;
    let solid = ogeom_algo::make_solid(model, &[shell])?.shape;
    super::guards::solid_fits(model, &solid, tol)?;
    Ok(solid)
}

/// A circular strip's spandrel, closed by its two face flaps and its exact
/// end sections. The section and sheet builders are shared with P1/P4.
pub(super) fn strip(
    model: &mut Model,
    seat: &crate::support::Seat,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let section = PlanarSection::new(seat, &RoundSpec::Circular(radius), tol)?;
    let length = section.length();
    let [a, b] = section.order();
    let points = [
        seat.start,
        seat.end,
        section.contact(a, 0.0),
        section.contact(b, 0.0),
        section.contact(a, length),
        section.contact(b, length),
    ];
    let mut edges: Vec<SheetEdge> = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (3, 5)]
        .into_iter()
        .map(|ends| SheetEdge { ends, curve: None })
        .collect();
    edges.push(SheetEdge {
        ends: (2, 3),
        curve: Some(section.arc_at(0.0, tol)?),
    });
    edges.push(SheetEdge {
        ends: (4, 5),
        curve: Some(section.arc_at(length, tol)?),
    });
    let face = |ring| SheetFace {
        surface: None,
        ring,
        holes: Vec::new(),
    };
    let mut faces = vec![
        face(vec![(0, true), (3, true), (5, false), (1, false)]),
        face(vec![(0, false), (2, true), (6, true), (4, false)]),
        face(vec![(1, true), (7, true), (2, false)]),
        face(vec![(4, true), (8, false), (3, false)]),
    ];
    faces.push(SheetFace {
        surface: Some(section.surface((0.0, length), tol)?),
        ring: vec![(7, true), (6, true), (8, false), (5, false)],
        holes: Vec::new(),
    });
    let mut built = build_faces(model, &points, &edges, &faces, tol)?;
    let Some(band) = built.pop() else {
        ogeom_bail!(Construction, "a corner's strip has no band");
    };
    built.push(band.reversed());
    let shell = ogeom_algo::make_shell(model, &built)?.shape;
    let tool = ogeom_algo::make_solid(model, &[shell])?.shape;
    // A spandrel tapers to zero at each rail: the planar flap and the
    // inward cylinder are tangent there with opposite outward normals.
    // P5's no-fold rule applies to the retained solid, not this cutter.
    let checked = ogeom_algo::check(model, &tool, tol)?;
    if !checked.is_valid() {
        ogeom_bail!(
            Construction,
            "a corner strip's cutter fails its check: {:?}",
            checked.problems
        );
    }
    Ok(tool)
}

pub(super) fn cut_strip(
    model: &mut Model,
    solid: &Shape,
    seat: &crate::support::Seat,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let tool = strip(model, seat, radius, tol)?;
    ogeom_bool::cut(model, solid, &tool, tol)
}

/// A truncated radial compartment outside the corner sphere. Its side
/// planes contain the ball centre; its far cap lies entirely past the
/// sphere. The input solid supplies the curved host boundaries at the cut.
pub(super) fn cone(
    model: &mut Model,
    centre: Point,
    radius: f64,
    rays: &[Vector],
    tips: &[Point],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let n = rays.len();
    let mut feet: Vec<Point> = rays.iter().map(|r| centre + *r * radius).collect();
    let mut tips = tips.to_vec();
    let middle = tips.iter().fold(Vector::ZERO, |sum, p| sum + (*p - centre));
    let turn = (tips[1] - tips[0]).cross(tips[2] - tips[0]);
    if turn.dot(middle) < 0.0 {
        feet.reverse();
        tips.reverse();
    }
    let points: Vec<Point> = feet.iter().chain(&tips).copied().collect();
    let mut edges = Vec::new();
    for k in 0..n {
        edges.push(SheetEdge {
            ends: (k, n + k),
            curve: None,
        });
    }
    for k in 0..n {
        edges.push(SheetEdge {
            ends: (n + k, n + (k + 1) % n),
            curve: None,
        });
    }
    for k in 0..n {
        let (x, y) = (feet[k] - centre, feet[(k + 1) % n] - centre);
        let frame = Frame::new(
            centre,
            Direction::new(x.cross(y), tol)?,
            Direction::new(x, tol)?,
            tol,
        )?;
        let angle = x.cross(y).magnitude().atan2(x.dot(y));
        let curve = CircleCurve::new(Circle::new(frame, radius, tol)?);
        edges.push(SheetEdge {
            ends: (k, (k + 1) % n),
            curve: Some((curve.into(), (0.0, angle))),
        });
    }
    let mut faces = vec![SheetFace {
        surface: None,
        ring: (0..n).map(|k| (n + k, true)).collect(),
        holes: Vec::new(),
    }];
    for k in 0..n {
        faces.push(SheetFace {
            surface: None,
            ring: vec![
                (n + k, false),
                (k, false),
                (2 * n + k, true),
                ((k + 1) % n, true),
            ],
            holes: Vec::new(),
        });
    }
    let surface = crate::round::sphere::corner_sphere(centre, radius, &feet, tol)?;
    faces.push(SheetFace {
        surface: Some(surface),
        ring: (0..n).map(|k| (2 * n + k, true)).collect(),
        holes: Vec::new(),
    });
    let mut built = build_faces(model, &points, &edges, &faces, tol)?;
    let Some(sphere) = built.pop() else {
        ogeom_bail!(Construction, "a radial compartment has no spherical face");
    };
    built.push(sphere.reversed());
    let shell = ogeom_algo::make_shell(model, &built)?.shape;
    let solid = ogeom_algo::make_solid(model, &[shell])?.shape;
    super::guards::solid_fits(model, &solid, tol)?;
    Ok(solid)
}
