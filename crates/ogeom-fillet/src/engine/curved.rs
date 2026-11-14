//! P0 for solids with curved faces: the solid read once into faces on
//! their surfaces, each with every loop it walks, and edges on their
//! curves, each with its ends and the faces walking it; and the sheet
//! assembly that builds faces back from such rings over shared edges.
//!
//! Unlike the planar read, an edge may close on itself (a rim, its ends
//! one vertex) and a face may walk one edge twice, once each way (a
//! cylinder's seam), and a face may hold holes. Every face keeps its own
//! surface, so a face built back from the same rings is the same face.

use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Orientation, Shape, ShapeType, explore_unique};

use super::{Id, id};

/// A face of the solid: its surface where it stands, its loops in the order
/// the face stores them, each a cycle of edges with the way it walks each,
/// and whether it stands reversed on its surface in the solid.
#[derive(Debug, Clone)]
pub(crate) struct CurvedFace {
    pub shape: Shape,
    pub surface: SurfaceGeometry,
    pub loops: Vec<Vec<(Id, bool)>>,
    pub reversed: bool,
}

/// An edge of the solid: its curve where it stands and the range of it the
/// edge covers, from its first vertex to its last (the same vertex for a
/// closed edge).
#[derive(Debug, Clone)]
pub(crate) struct CurvedEdge {
    pub shape: Shape,
    pub curve: Curve,
    pub range: (f64, f64),
    pub ends: (Id, Id),
}

/// A solid read for the engine.
#[derive(Debug, Clone)]
pub(crate) struct Curved {
    pub faces: Vec<CurvedFace>,
    pub edges: FastMap<Id, CurvedEdge>,
    pub points: FastMap<Id, Point>,
    pub shapes: FastMap<Id, Shape>,
}

/// `solid` read into faces and edges, or `None` where it is not one solid of
/// one shell, a shape in it is placed (its node then stands in two places),
/// or an edge has no curve. A degenerate edge (a sphere's pole) is left out
/// of the faces' loops, and [`build`] puts one back wherever a face's chart
/// breaks at a pole.
///
/// # Errors
///
/// As reading the shapes' geometry.
pub(crate) fn read(model: &Model, solid: &Shape, tol: Tolerances) -> OgeomResult<Option<Curved>> {
    if model.kind_of(solid)? != ShapeType::Solid
        || explore_unique(model, solid, ShapeType::Shell)?.len() != 1
    {
        return Ok(None);
    }
    let placed = |s: &Shape| !s.location().is_identity();
    let mut points: FastMap<Id, Point> = FastMap::default();
    let mut shapes: FastMap<Id, Shape> = FastMap::default();
    for v in explore_unique(model, solid, ShapeType::Vertex)? {
        if placed(&v) {
            return Ok(None);
        }
        let Some(p) = model
            .node(&v)
            .and_then(|n| n.data().as_vertex().map(|d| d.point))
        else {
            return Ok(None);
        };
        points.insert(id(&v), p);
        shapes.insert(id(&v), v);
    }
    let mut edges: FastMap<Id, CurvedEdge> = FastMap::default();
    let mut poles: Vec<Id> = Vec::new();
    for e in explore_unique(model, solid, ShapeType::Edge)? {
        if placed(&e) {
            return Ok(None);
        }
        let e = e.oriented(Orientation::Forward);
        if model
            .node(&e)
            .and_then(|n| n.data().as_edge())
            .is_some_and(|d| d.degenerate)
        {
            poles.push(id(&e));
            continue;
        }
        let Some((a, b)) = ogeom_algo::edge_vertices(model, &e)? else {
            return Ok(None);
        };
        let Some((curve, range)) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
            return Ok(None);
        };
        edges.insert(
            id(&e),
            CurvedEdge {
                shape: e.clone(),
                curve,
                range,
                ends: (id(&a), id(&b)),
            },
        );
        shapes.insert(id(&e), e);
    }
    let mut faces = Vec::new();
    for shape in explore_unique(model, solid, ShapeType::Face)? {
        if placed(&shape) {
            return Ok(None);
        }
        let Some(data) = model.node(&shape).and_then(|n| n.data().as_face()) else {
            return Ok(None);
        };
        if !data.location.is_identity() {
            return Ok(None);
        }
        let Some(surface) = model.geometry().surface(data.surface).cloned() else {
            return Ok(None);
        };
        let mut loops = Vec::new();
        for wire in model.children_of(&shape)? {
            let walk: Vec<(Id, bool)> = model
                .children_of(&wire)?
                .iter()
                .filter(|e| !poles.contains(&id(e)))
                .map(|e| (id(e), e.orientation() == Orientation::Forward))
                .collect();
            loops.push(walk);
        }
        faces.push(CurvedFace {
            reversed: shape.orientation() == Orientation::Reversed,
            shape,
            surface,
            loops,
        });
    }
    Ok(Some(Curved {
        faces,
        edges,
        points,
        shapes,
    }))
}

/// One face to build: its surface, its rings over the sheet's edges (each
/// edge with the way it is walked, in the surface's own sense), and whether
/// it is turned over once built.
#[derive(Debug, Clone)]
pub(crate) struct SheetRings {
    pub surface: SurfaceGeometry,
    pub rings: Vec<Vec<(usize, bool)>>,
    pub reversed: bool,
}

/// An edge of the sheet: its ends, as indices of the sheet's points, and
/// its curve and range from the first to the second.
#[derive(Debug, Clone)]
pub(crate) struct SheetCurve {
    pub ends: (usize, usize),
    pub curve: Curve,
    pub range: (f64, f64),
}

/// What a sheet is built into: its vertices, edges and faces, each in the
/// sheet's order.
#[derive(Debug, Clone)]
pub(crate) struct Sheet {
    pub vertices: Vec<Shape>,
    pub edges: Vec<Shape>,
    pub faces: Vec<Shape>,
}

/// The faces of a sheet over shared vertices and edges, in the sheet's
/// order, each on its surface with the pcurves its edges need.
///
/// # Errors
///
/// As making the edges, wires and faces, and fitting a bare edge's pcurve.
pub(crate) fn build(
    model: &mut Model,
    points: &[Point],
    edges: &[SheetCurve],
    faces: &[SheetRings],
    tol: Tolerances,
) -> OgeomResult<Sheet> {
    let vertices: Vec<Shape> = points
        .iter()
        .map(|p| ogeom_algo::make_vertex(model, *p).shape)
        .collect();
    let mut made = Vec::with_capacity(edges.len());
    for e in edges {
        let (a, b) = e.ends;
        made.push(
            ogeom_algo::make_edge_between(
                model,
                e.curve.clone(),
                e.range,
                &vertices[a],
                &vertices[b],
                tol,
            )?
            .shape,
        );
    }
    let mut out = Vec::with_capacity(faces.len());
    for face in faces {
        let mut wires = Vec::with_capacity(face.rings.len());
        for ring in &face.rings {
            let walked: Vec<Shape> = ring
                .iter()
                .map(|(e, forward)| {
                    if *forward {
                        made[*e].clone()
                    } else {
                        made[*e].reversed()
                    }
                })
                .collect();
            wires.push(ogeom_algo::make_wire(model, &walked, tol)?.shape);
        }
        let mut built = ogeom_algo::make_face(model, face.surface.clone(), &wires, tol)?.shape;
        ogeom_algo::trimmed_or_fitted_where_bare(model, &built, tol)?;
        // A ring through a pole (a sphere's seam walked to it and back)
        // meets itself there in space but not in the chart: the face is
        // built again with the pole's degenerate edge closing it.
        let surface = model
            .node(&built)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface);
        if let Some(surface) = surface {
            let mut closed = Vec::with_capacity(face.rings.len());
            let mut grew = false;
            for wire in &wires {
                let uses = model.children_of(wire)?;
                let with_poles = ogeom_algo::closed_at_poles(model, surface, &uses, tol)?;
                grew |= with_poles.len() > uses.len();
                closed.push(with_poles);
            }
            if grew {
                let wires = closed
                    .iter()
                    .map(|uses| Ok(ogeom_algo::make_wire(model, uses, tol)?.shape))
                    .collect::<OgeomResult<Vec<_>>>()?;
                built = ogeom_algo::make_face_on(model, surface, &wires, tol)?.shape;
            }
        }
        out.push(if face.reversed {
            built.reversed()
        } else {
            built
        });
    }
    Ok(Sheet {
        vertices,
        edges: made,
        faces: out,
    })
}

/// One edge of a curved sheet: its two ends (indices into the sheet's
/// points) and its curve with the range between them, or `None` for the
/// straight segment.
#[derive(Debug, Clone)]
pub(crate) struct SheetEdge {
    pub ends: (usize, usize),
    pub curve: Option<(Curve, (f64, f64))>,
}

/// One face of a curved sheet: its surface (`None` for the plane its ring
/// spans), its ring of edges each with the way it is walked, wound
/// counter-clockwise seen from outside the solid.
#[derive(Debug, Clone)]
pub(crate) struct SheetFace {
    pub surface: Option<SurfaceGeometry>,
    pub ring: Vec<(usize, bool)>,
}

/// The faces of a sheet built over shared vertices and edges, in the
/// sheet's order, by [`build`]: each `None` curve the straight segment
/// between its ends, each `None` surface the plane its ring spans, each
/// face its one ring. The faces are left for the caller to assemble with
/// the rest.
///
/// # Errors
///
/// As [`build`], and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// a planar ring spans no plane.
pub(crate) fn build_faces(
    model: &mut Model,
    points: &[Point],
    edges: &[SheetEdge],
    faces: &[SheetFace],
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let curves = edges
        .iter()
        .map(|e| {
            let (a, b) = e.ends;
            Ok(match &e.curve {
                Some((curve, range)) => SheetCurve {
                    ends: e.ends,
                    curve: curve.clone(),
                    range: *range,
                },
                None => SheetCurve {
                    ends: e.ends,
                    curve: Curve::Line(ogeom_geom::LineCurve::segment(points[a], points[b], tol)?),
                    range: (0.0, points[a].distance(points[b])),
                },
            })
        })
        .collect::<OgeomResult<Vec<_>>>()?;
    let rings = faces
        .iter()
        .map(|face| {
            let surface = match &face.surface {
                Some(surface) => surface.clone(),
                None => {
                    let corners: Vec<Point> = face
                        .ring
                        .iter()
                        .map(|(e, forward)| {
                            let (a, b) = edges[*e].ends;
                            points[if *forward { a } else { b }]
                        })
                        .collect();
                    let mut n = Vector::ZERO;
                    for k in 1..corners.len().saturating_sub(1) {
                        n += (corners[k] - corners[0]).cross(corners[k + 1] - corners[0]);
                    }
                    let normal = ogeom_math::Direction::new(n, tol)?;
                    let x = ogeom_math::Direction::new(corners[1] - corners[0], tol)?;
                    let plane =
                        ogeom_math::Plane::new(ogeom_math::Frame::new(corners[0], normal, x, tol)?);
                    ogeom_geom::PlaneSurface::new(plane).into()
                }
            };
            Ok(SheetRings {
                surface,
                rings: vec![face.ring.clone()],
                reversed: false,
            })
        })
        .collect::<OgeomResult<Vec<_>>>()?;
    Ok(build(model, points, &curves, &rings, tol)?.faces)
}

/// `solid` read and built back face for face from what was read: the same
/// solid on new nodes, and the history from its faces, edges and vertices
/// to theirs.
///
/// # Errors
///
/// As [`read`] and [`build`], and as making the shell and solid.
#[cfg(test)]
pub(crate) fn rebuilt(
    model: &mut Model,
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_algo::Built>> {
    let Some(read) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let mut ids: Vec<Id> = read.points.keys().copied().collect();
    ids.sort_unstable();
    let index: FastMap<Id, usize> = ids.iter().enumerate().map(|(i, v)| (*v, i)).collect();
    let points: Vec<Point> = ids.iter().map(|v| read.points[v]).collect();
    let mut edge_ids: Vec<Id> = read.edges.keys().copied().collect();
    edge_ids.sort_unstable();
    let edge_index: FastMap<Id, usize> =
        edge_ids.iter().enumerate().map(|(i, e)| (*e, i)).collect();
    let edges: Vec<SheetCurve> = edge_ids
        .iter()
        .map(|e| {
            let edge = &read.edges[e];
            SheetCurve {
                ends: (index[&edge.ends.0], index[&edge.ends.1]),
                curve: edge.curve.clone(),
                range: edge.range,
            }
        })
        .collect();
    // A face reversed in the solid reads each edge walked the other way,
    // its loops in their stored order: on its surface's own sense each walk
    // turns back.
    let faces: Vec<SheetRings> = read
        .faces
        .iter()
        .map(|f| SheetRings {
            surface: f.surface.clone(),
            rings: f
                .loops
                .iter()
                .map(|l| {
                    l.iter()
                        .map(|(e, w)| (edge_index[e], *w != f.reversed))
                        .collect()
                })
                .collect(),
            reversed: f.reversed,
        })
        .collect();
    let sheet = build(model, &points, &edges, &faces, tol)?;
    let shell = ogeom_algo::make_shell(model, &sheet.faces)?.shape;
    let solid_out = ogeom_algo::make_solid(model, &[shell])?.shape;
    let mut history = ogeom_algo::History::new();
    for (f, b) in read.faces.iter().zip(&sheet.faces) {
        history.modify(&f.shape, b.clone());
    }
    for (e, b) in edge_ids.iter().zip(&sheet.edges) {
        history.modify(&read.edges[e].shape, b.clone());
    }
    for (v, b) in ids.iter().zip(&sheet.vertices) {
        history.modify(&read.shapes[v], b.clone());
    }
    Ok(Some(ogeom_algo::Built::new(solid_out, history)))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::{SheetEdge, SheetFace, build_faces};
    use ogeom_core::Tolerances;
    use ogeom_geom::{CircleCurve, Curve, CylinderSurface, LineCurve, SurfaceGeometry};
    use ogeom_math::{Circle, Cylinder, Direction, Frame, Point, Vector};
    use ogeom_topo::{Model, ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    /// A drum of radius 2 and height 3 as a sheet: two discs, and the
    /// wall's ring the bottom circle, the seam up, the top circle back and
    /// the seam down.
    #[test]
    fn a_sheet_with_a_seam_and_closed_edges_builds_a_drum() {
        let frame = |z: f64| {
            Frame::new(
                Point::new(0.0, 0.0, z),
                Direction::new(Vector::new(0.0, 0.0, 1.0), T).unwrap(),
                Direction::new(Vector::new(1.0, 0.0, 0.0), T).unwrap(),
                T,
            )
            .unwrap()
        };
        let circle =
            |z: f64| Curve::Circle(CircleCurve::new(Circle::new(frame(z), 2.0, T).unwrap()));
        let points = vec![Point::new(2.0, 0.0, 0.0), Point::new(2.0, 0.0, 3.0)];
        let tau = core::f64::consts::TAU;
        let edges = vec![
            SheetEdge {
                ends: (0, 0),
                curve: Some((circle(0.0), (0.0, tau))),
            },
            SheetEdge {
                ends: (1, 1),
                curve: Some((circle(3.0), (0.0, tau))),
            },
            SheetEdge {
                ends: (0, 1),
                curve: Some((
                    Curve::Line(LineCurve::segment(points[0], points[1], T).unwrap()),
                    (0.0, 3.0),
                )),
            },
        ];
        let wall: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(frame(0.0), 2.0, T).unwrap(), (0.0, 3.0))
                .unwrap()
                .into();
        // The bottom disc's plane faces into the drum: built on its natural
        // side, then turned over below.
        let disc = |z: f64| -> SurfaceGeometry {
            ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(frame(z))).into()
        };
        let faces = vec![
            SheetFace {
                surface: Some(disc(0.0)),
                ring: vec![(0, true)],
            },
            SheetFace {
                surface: Some(disc(3.0)),
                ring: vec![(1, true)],
            },
            SheetFace {
                surface: Some(wall),
                ring: vec![(0, true), (2, true), (1, false), (2, false)],
            },
        ];
        let mut m = Model::new();
        let mut built = build_faces(&mut m, &points, &edges, &faces, T).unwrap();
        built[0] = built[0].reversed();
        let shell = ogeom_algo::make_shell(&mut m, &built).unwrap().shape;
        let solid = ogeom_algo::make_solid(&mut m, &[shell]).unwrap().shape;
        let found = ogeom_algo::check(&m, &solid, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let v = ogeom_algo::volume_properties(
            &m,
            &solid,
            ogeom_mesh::Deflection {
                chord: 1e-5,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap()
        .mass;
        assert!((v - core::f64::consts::PI * 4.0 * 3.0).abs() < 1e-6, "{v}");
    }

    fn volume(m: &Model, s: &ogeom_topo::Shape) -> f64 {
        ogeom_algo::volume_properties(
            m,
            s,
            ogeom_mesh::Deflection {
                chord: 1e-5,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap()
        .mass
    }

    /// Read and built back face for face, a solid is the same solid: valid,
    /// its volume to the bit's neighbourhood, its faces on the same
    /// surfaces.
    fn round_trips(m: &mut Model, s: &ogeom_topo::Shape) {
        let before = volume(m, s);
        let built = super::rebuilt(m, s, T).unwrap().unwrap();
        let found = ogeom_algo::check(m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let after = volume(m, &built.shape);
        assert!(
            (after - before).abs() <= before * 1e-9,
            "{after} against {before}"
        );
        for kind in [ShapeType::Face, ShapeType::Edge, ShapeType::Vertex] {
            let made = explore_unique(m, &built.shape, kind).unwrap();
            for x in explore_unique(m, s, kind).unwrap() {
                let images = built.history.trace(&x);
                assert_eq!(images.len(), 1, "{kind:?}");
                assert!(made.iter().any(|y| y.is_same(&images[0])), "{kind:?}");
            }
        }
    }

    /// A cylinder (seam and closed rims), a cone, a box with a hole drilled
    /// through it (a face with a hole, a bore's wall) and a box with a boss
    /// on it read and build back whole.
    #[test]
    fn curved_solids_read_and_build_back_whole() {
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, 2.0, 3.0, T)
            .unwrap()
            .shape;
        round_trips(&mut m, &drum);
        let cone = ogeom_algo::make_cone(&mut m, Frame::WORLD, 3.0, 1.0, 4.0, T)
            .unwrap()
            .shape;
        round_trips(&mut m, &cone);
        let block = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 4.0), T)
            .unwrap()
            .shape;
        let at = |x: f64, y: f64, z: f64| {
            Frame::new(
                Point::new(x, y, z),
                Direction::new(Vector::new(0.0, 0.0, 1.0), T).unwrap(),
                Direction::new(Vector::new(1.0, 0.0, 0.0), T).unwrap(),
                T,
            )
            .unwrap()
        };
        let bore = ogeom_algo::make_cylinder(&mut m, at(5.0, 4.0, -1.0), 1.5, 6.0, T)
            .unwrap()
            .shape;
        let holed = ogeom_bool::cut(&mut m, &block, &bore, T).unwrap().shape;
        round_trips(&mut m, &holed);
        let boss = ogeom_algo::make_cylinder(&mut m, at(5.0, 4.0, 4.0), 1.5, 2.0, T)
            .unwrap()
            .shape;
        let bossed = ogeom_bool::fuse(&mut m, &block, &boss, T).unwrap().shape;
        round_trips(&mut m, &bossed);
    }

    /// A 10 cube with its top front edge rounded at 2, built from a sheet:
    /// six planar faces, the cylinder strip between the contact lines and
    /// the two arcs ending it. Valid, and the cube less the round's
    /// section (`(1 - pi/4) r^2`) along its 10.
    #[test]
    fn a_sheet_with_a_cylinder_strip_builds_the_rounded_cube() {
        let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
        let points = vec![
            p(0.0, 0.0, 0.0),
            p(10.0, 0.0, 0.0),
            p(10.0, 10.0, 0.0),
            p(0.0, 10.0, 0.0),
            p(0.0, 0.0, 8.0),
            p(10.0, 0.0, 8.0),
            p(10.0, 2.0, 10.0),
            p(0.0, 2.0, 10.0),
            p(10.0, 10.0, 10.0),
            p(0.0, 10.0, 10.0),
        ];
        let arc = |x: f64| -> Curve {
            let frame = Frame::new(
                p(x, 2.0, 8.0),
                Direction::new(Vector::new(-1.0, 0.0, 0.0), T).unwrap(),
                Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
                T,
            )
            .unwrap();
            Curve::Circle(ogeom_geom::CircleCurve::new(
                Circle::new(frame, 2.0, T).unwrap(),
            ))
        };
        let line = |a: usize, b: usize| SheetEdge {
            ends: (a, b),
            curve: None,
        };
        let edges = vec![
            line(0, 1), // 0
            line(1, 2), // 1
            line(2, 3), // 2
            line(3, 0), // 3
            line(0, 4), // 4
            line(1, 5), // 5
            line(2, 8), // 6
            line(3, 9), // 7
            line(4, 5), // 8: the front contact
            line(7, 6), // 9: the top contact
            line(6, 8), // 10
            line(8, 9), // 11
            line(9, 7), // 12
            SheetEdge {
                ends: (4, 7),
                curve: Some((arc(0.0), (0.0, std::f64::consts::FRAC_PI_2))),
            }, // 13
            SheetEdge {
                ends: (5, 6),
                curve: Some((arc(10.0), (0.0, std::f64::consts::FRAC_PI_2))),
            }, // 14
        ];
        let plane = |ring: Vec<(usize, bool)>| SheetFace {
            surface: None,
            ring,
        };
        let axis = Frame::new(
            p(0.0, 2.0, 8.0),
            Direction::new(Vector::new(1.0, 0.0, 0.0), T).unwrap(),
            Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
            T,
        )
        .unwrap();
        let cylinder =
            ogeom_geom::CylinderSurface::new(Cylinder::new(axis, 2.0, T).unwrap(), (-1.0, 11.0))
                .unwrap();
        let faces = vec![
            plane(vec![(3, false), (2, false), (1, false), (0, false)]), // bottom
            plane(vec![(9, true), (10, true), (11, true), (12, true)]),  // top
            plane(vec![(0, true), (5, true), (8, false), (4, false)]),   // front
            plane(vec![(7, true), (11, false), (6, false), (2, true)]),  // back
            plane(vec![
                (4, true),
                (13, true),
                (12, false),
                (7, false),
                (3, true),
            ]), // left
            plane(vec![
                (1, true),
                (6, true),
                (10, false),
                (14, false),
                (5, false),
            ]), // right
            SheetFace {
                surface: Some(SurfaceGeometry::Cylinder(cylinder)),
                ring: vec![(8, true), (14, true), (9, false), (13, false)],
            },
        ];
        let mut m = Model::new();
        let built = build_faces(&mut m, &points, &edges, &faces, T).unwrap();
        let shell = ogeom_algo::make_shell(&mut m, &built).unwrap().shape;
        let solid = ogeom_algo::make_solid(&mut m, &[shell]).unwrap().shape;
        let found = ogeom_algo::check(&m, &solid, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let volume = ogeom_algo::volume_properties(
            &m,
            &solid,
            ogeom_mesh::Deflection {
                chord: 1e-5,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap()
        .mass;
        let want = 1000.0 - (1.0 - std::f64::consts::FRAC_PI_4) * 4.0 * 10.0;
        assert!((volume - want).abs() < 1e-6, "{volume} against {want}");
    }
}
