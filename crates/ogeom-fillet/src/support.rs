//! Shared scaffolding for the subtractive blends.
//!
//! A chamfer and a constant-radius fillet on a straight edge between planar
//! faces differ only in the face that replaces the edge: a bevel plane for
//! one, a tangent cylinder for the other. Everything around that face (the
//! seat on the solid, the legs running along the adjacent faces, faces
//! assembled from explicit curves with exact pcurves) is one piece of
//! scaffolding, kept here so the two operations cannot drift apart.

use ogeom_algo::{make_edge, make_edge_between, make_revolution_band, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{CircleCurve, Curve, CylinderSurface, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Orientation, Shape, ShapeType, explore};

/// Where a blend sits on a solid: a straight edge and the two planar faces
/// meeting there, reduced to the numbers the wedge construction runs on.
pub(crate) struct Seat {
    /// The edge's start, at its lower parameter.
    pub start: Point,
    /// The edge's end.
    pub end: Point,
    /// Unit direction from start to end.
    pub along: Vector,
    /// Outward unit normals of the two faces, in discovery order.
    pub normals: [Vector; 2],
    /// The two faces themselves, in the same order, so a caller naming a
    /// face can find which leg it owns.
    pub faces: [Shape; 2],
    /// Whether the edge is convex: material inside the dihedral, so a blend
    /// subtracts. A concave edge's blend adds, with every sign mirrored.
    pub convex: bool,
}

impl Seat {
    /// On face `i`, the unit direction perpendicular to the edge that walks
    /// away from the *other* face, into the material the blend cuts back
    /// along.
    pub fn leg(&self, i: usize, tol: Tolerances) -> OgeomResult<Vector> {
        let own = self.normals[i];
        let other = self.normals[1 - i];
        let mut t = own.cross(self.along);
        if t.dot(other) > 0.0 {
            t = -t;
        }
        let m = t.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "a face is tangent to its own edge");
        }
        Ok(t / m)
    }
}

/// Whether `candidate` is the same edge occurrence as `edge`: the same node,
/// placed the same way in the world. Two references through different
/// location paths that compose to one placement are the same occurrence; the
/// same node re-instanced elsewhere (a prism's profile edge at both caps)
/// is not.
pub(crate) fn same_occurrence(
    model: &Model,
    candidate: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> bool {
    if candidate.node() != edge.node() {
        return false;
    }
    if candidate.is_same(edge) {
        return true;
    }
    let (Ok(a), Ok(b)) = (
        candidate.transform(model.datums()),
        edge.transform(model.datums()),
    ) else {
        return false;
    };
    use ogeom_math::Point;
    [
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
        Point::new(0.0, 0.0, 1.0),
    ]
    .into_iter()
    .all(|p| a.apply(p).distance(b.apply(p)) <= tol.confusion())
}

/// An edge's 3D curve and range, placed in the world.
pub(crate) fn edge_curve(
    model: &Model,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<(Curve, (f64, f64))> {
    let Some(node) = model.node(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data() else {
        ogeom_bail!(Construction, "expected an edge");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "the edge has no curve to blend along");
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    // Baked into the world: the edge lives wherever its placement puts it,
    // and every seat is measured there. A rigid placement preserves the
    // parameterization, so the stored range carries over unchanged.
    use ogeom_geom::Transformable as _;
    let placed = geometry
        .clone()
        .transformed(&edge.transform(model.datums())?, tol)?;
    Ok((placed, *range))
}

/// Whether every face of `solid` bounding `edge` is planar: the straight
/// seat's own forms speak only planes, and a straight edge on a curved
/// face (a ruling of a drum where it meets a wall) is the march's.
pub(crate) fn hosts_planar(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| same_occurrence(model, e, edge, tol));
        if !touches {
            continue;
        }
        let planar = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, SurfaceGeometry::Plane(_)));
        if !planar {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Find the seat of a blend: the straight edge's ends and direction, and the
/// outward normals of the exactly two planar faces of `solid` meeting there.
///
/// Refuses tangent edges, which have no corner. Convexity is read from the
/// solid and carried in the seat: a convex edge's wedge subtracts, a
/// concave edge's adds, with every sign mirrored.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the edge is
/// not straight, is not shared by exactly two planar faces of `solid`, or the
/// edge is tangent.
pub(crate) fn planar_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Seat> {
    let (curve, range) = edge_curve(model, edge, tol)?;
    let Curve::Line(_) = &curve else {
        ogeom_bail!(
            Construction,
            "blending a curved edge needs the marching blend machinery; this \
             is the straight-edge form"
        );
    };
    let start = curve.point_at(range.0, tol)?;
    let end = curve.point_at(range.1, tol)?;
    let along = (end - start) / start.distance(end);

    let mut normals: Vec<Vector> = Vec::new();
    let mut faces: Vec<Shape> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| same_occurrence(model, e, edge, tol));
        if !touches {
            continue;
        }
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
            ogeom_bail!(
                Construction,
                "blending an edge of a curved face needs the marching blend \
                 machinery; this is the planar form"
            );
        };
        let placement = face.transform(model.datums())?;
        let mut normal = placement.apply_vector(plane.plane().normal().vector());
        if face.orientation() == ogeom_topo::Orientation::Reversed {
            normal = -normal;
        }
        normals.push(normal);
        faces.push(face);
    }
    if normals.len() != 2 {
        ogeom_bail!(
            Construction,
            "a blend needs an edge shared by exactly two faces, found {}",
            normals.len()
        );
    }

    let mut seat = Seat {
        start,
        end,
        along,
        normals: [normals[0], normals[1]],
        faces: [faces[0].clone(), faces[1].clone()],
        convex: true,
    };
    // Convexity is read from the face itself, not derived from the normals:
    // the leg construction cannot answer it, because it *chooses* its side.
    // Which way the first face actually extends from the edge (sampled
    // against its own trim) leans behind the other face's plane on a convex
    // edge and in front of it on a concave one.
    let raw = {
        let t = normals[0].cross(along);
        let m = t.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "a face is tangent to its own edge");
        }
        t / m
    };
    let mid = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
    let mut face_side: Option<Vector> = None;
    'scales: for scale in [1e-3, 1e-2, 5e-2] {
        let eps = start.distance(end) * scale;
        let deflection = ogeom_mesh::Deflection {
            chord: eps * 0.1,
            ..ogeom_mesh::Deflection::default()
        };
        for dir in [raw, -raw] {
            if crate::support::on_face_side(
                model,
                &seat.faces[0],
                mid + dir * eps,
                deflection,
                tol,
            )? {
                face_side = Some(dir);
                break 'scales;
            }
        }
    }
    let Some(extends) = face_side else {
        ogeom_bail!(
            Construction,
            "cannot read which way the edge's face extends; the face is \
             thinner than the probe can resolve"
        );
    };
    let lean = extends.dot(seat.normals[1]);
    if lean.abs() <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the edge's faces are tangent; there is no corner to blend"
        );
    }
    seat.convex = lean < 0.0;
    Ok(seat)
}

/// Where a revolved blend sits: a circular rim shared by one perpendicular
/// planar cap and one coaxial cylindrical wall, reduced to the numbers the
/// revolved wedge runs on.
///
/// Four seats, one parameterization. With `sigma` the wall's outward radial
/// sign and `tau` telling whether the wall extends away from the cap's
/// outward side, `tau` alone decides whether the wedge subtracts (the
/// external rim and the hole's rim, both convex) or fuses (the boss base and
/// the blind hole's floor, both concave).
pub(crate) struct RevolvedSeat {
    /// The rim's centre.
    pub centre: Point,
    /// The rim's radius, also the wall's.
    pub radius: f64,
    /// The cap's outward unit normal.
    pub up: Vector,
    /// The rim frame's `x`, so every ring built on the seat shares a
    /// parameter origin.
    pub x_ref: Direction,
    /// The wall's outward radial sign: `+1` material inside, `-1` a bore.
    pub sigma: f64,
    /// `+1` when the wall extends away from the cap's outward side (the rim
    /// configurations) and `-1` alongside it, the concave seats.
    pub tau: f64,
    /// The planar cap.
    pub cap_face: Shape,
    /// The cylindrical wall.
    pub wall_face: Shape,
}

impl RevolvedSeat {
    /// Whether the wedge fuses rather than subtracts.
    pub const fn additive(&self) -> bool {
        self.tau < 0.0
    }

    /// A frame at `origin` sharing the seat's axis and parameter origin.
    pub fn frame_at(&self, origin: Point, tol: Tolerances) -> OgeomResult<Frame> {
        Frame::new(origin, Direction::new(self.up, tol)?, self.x_ref, tol)
    }

    /// A full ring on the seat's axis: radius `r` in the plane through
    /// `origin`.
    pub fn ring(
        &self,
        model: &mut Model,
        origin: Point,
        r: f64,
        tol: Tolerances,
    ) -> OgeomResult<Shape> {
        let circle = Circle::new(self.frame_at(origin, tol)?, r, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        let domain = ogeom_geom::Curve3d::domain(&curve);
        Ok(make_edge(model, curve, domain, tol)?.shape)
    }
}

/// Find a revolved blend's seat: the exactly one perpendicular planar cap and
/// one coaxial cylindrical wall meeting at a circular rim.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// rim's faces are not that pair, or the wall has no far ring to read its
/// side from.
pub(crate) fn revolved_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    rim: &CircleCurve,
    tol: Tolerances,
) -> OgeomResult<RevolvedSeat> {
    let rim_circle = rim.circle();
    let rim_centre = rim_circle.centre();
    let rim_radius = rim_circle.radius();
    let rim_axis = rim_circle.frame().z().vector();

    // The two faces at the rim: exactly one plane and one coaxial cylinder.
    let mut cap: Option<(Shape, Vector)> = None;
    let mut wall: Option<(Shape, f64)> = None;
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| same_occurrence(model, e, edge, tol));
        if !touches {
            continue;
        }
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let placement = face.transform(model.datums())?;
        let reversed = face.orientation() == Orientation::Reversed;
        match model.geometry().surface(data.surface) {
            Some(SurfaceGeometry::Plane(p)) => {
                let mut normal = placement.apply_vector(p.plane().normal().vector());
                if reversed {
                    normal = -normal;
                }
                if normal.cross(rim_axis).magnitude() > tol.angular()
                    || p.plane()
                        .distance_to(placement.inverse()?.apply(rim_centre))
                        > tol.confusion()
                {
                    ogeom_bail!(
                        Construction,
                        "the rim's cap is not the perpendicular plane through \
                         it; that seat needs the marching blend machinery"
                    );
                }
                cap = Some((face.clone(), normal));
            }
            Some(SurfaceGeometry::Cylinder(c)) => {
                let cyl = c.cylinder();
                let axis_point = placement.apply(cyl.frame().origin());
                let axis_z = placement.apply_vector(cyl.frame().z().vector());
                let off_axis = {
                    let to_rim = rim_centre - axis_point;
                    (to_rim - axis_z * to_rim.dot(axis_z)).magnitude()
                };
                if axis_z.cross(rim_axis).magnitude() > tol.angular()
                    || off_axis > tol.confusion()
                    || (cyl.radius() - rim_radius).abs() > tol.confusion()
                {
                    ogeom_bail!(
                        Construction,
                        "the rim's wall is not the coaxial cylinder through \
                         it; that seat needs the marching blend machinery"
                    );
                }
                let sigma = if reversed { -1.0 } else { 1.0 };
                wall = Some((face.clone(), sigma));
            }
            Some(_) => ogeom_bail!(
                Construction,
                "the rim meets a face that is neither plane nor cylinder; \
                 that seat needs the marching blend machinery"
            ),
            None => ogeom_bail!(Dangling, "face refers to a surface not in this model"),
        }
    }
    let (Some((cap_face, up_raw)), Some((wall_face, sigma))) = (cap, wall) else {
        ogeom_bail!(
            Construction,
            "a revolved blend needs the edge shared by one planar cap and one \
             cylindrical wall"
        );
    };
    let up = up_raw / up_raw.magnitude();

    // Which side of the cap the wall extends: read at the edge itself, just
    // above and just below the cap's plane, where the wall stands on one
    // side only; a wall that runs on past the cap elsewhere (a boss standing
    // over a block's edge, its wall going down beside the block) stands on
    // both sides of the plane somewhere, so its far rings alone do not say.
    // Where the probes do not decide, the wall band's other ring does.
    // `tau` positive means away from the cap's outward side (the rim
    // configurations) and negative means alongside it, the concave seats.
    let local_tau = {
        let (curve, range) = edge_curve(model, edge, tol)?;
        let at = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
        let eps = rim_radius * 1e-2;
        let deflection = ogeom_mesh::Deflection {
            chord: eps * 0.1,
            ..ogeom_mesh::Deflection::default()
        };
        let stands = |p: Point| -> OgeomResult<bool> {
            Ok(
                ogeom_algo::classify_on_face(model, &wall_face, p, deflection, tol)?
                    == ogeom_algo::Containment::In,
            )
        };
        match (stands(at + up * eps)?, stands(at - up * eps)?) {
            (true, false) => Some(-1.0),
            (false, true) => Some(1.0),
            _ => None,
        }
    };
    let tau = if let Some(tau) = local_tau {
        tau
    } else {
        let mut side = None;
        for e in explore(model, &wall_face, Filter::OfType(ShapeType::Edge))? {
            if same_occurrence(model, &e, edge, tol) {
                continue;
            }
            // Any circular edge of the wall at another height says which
            // side the wall extends: whether the ring survives as one
            // closed edge or as the arcs a boolean rebuilt it into.
            let Ok((curve, _)) = edge_curve(model, &e, tol) else {
                continue;
            };
            let Curve::Circle(c) = curve else {
                continue;
            };
            let lean = (c.circle().centre() - rim_centre).dot(up);
            if lean.abs() > tol.confusion() * 10.0 {
                side = Some(lean);
                break;
            }
        }
        let Some(lean) = side else {
            ogeom_bail!(
                Construction,
                "the rim's wall has no far ring to read its side from; the \
                 partial seat needs the marching blend machinery"
            );
        };
        if lean < 0.0 { 1.0 } else { -1.0 }
    };

    Ok(RevolvedSeat {
        centre: rim_centre,
        radius: rim_radius,
        up,
        x_ref: rim_circle.frame().x(),
        sigma,
        tau,
        cap_face,
        wall_face,
    })
}

/// Whether a blend along an open arc of a rim must run the whole turn.
///
/// An open arc is two different seats: a piece a boolean split off a full
/// rim, whose blend must run the whole turn, or a rim that genuinely stops,
/// a stadium's rounded end or a rounded corner. The two faces answer: the
/// rim goes on round the arc's complement only where the wall stands just
/// below it and the cap lies just beside it, on the cap's own side of the
/// rim. A wall that carries on past where the cap ends (a boss standing
/// over the edge of a block) holds no rim there.
pub(crate) fn runs_whole_turn(
    model: &Model,
    seat: &RevolvedSeat,
    arc: &CircleCurve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Curve3d as _;
    let frame = seat.frame_at(seat.centre, tol)?;
    let angle_of = |t: f64| -> OgeomResult<f64> {
        let local = frame.to_local(arc.point_at(t, tol)?);
        Ok(local.y.atan2(local.x))
    };
    // The arc's own middle: which way round it runs about the seat's axis
    // depends on whether the cap faces along the circle's normal or against
    // it, and the parameter's middle is the arc's middle either way.
    let mid = angle_of(f64::midpoint(range.0, range.1))?;
    let complement = mid + core::f64::consts::PI;
    let dir = frame.x().vector() * complement.cos() + frame.y().vector() * complement.sin();
    let eps = seat.radius * 1e-2;
    let probe = seat.centre + dir * seat.radius - seat.up * (seat.tau * eps);
    let deflection = ogeom_mesh::Deflection {
        chord: eps * 0.1,
        ..ogeom_mesh::Deflection::default()
    };
    if ogeom_algo::classify_on_face(model, &seat.wall_face, probe, deflection, tol)?
        != ogeom_algo::Containment::In
    {
        return Ok(false);
    }
    // The cap lies inside the rim where the material is inside the wall and
    // the wall runs away from the cap's outward side, or where both are
    // turned; outside it otherwise.
    let outward = -seat.sigma * seat.tau;
    let beside = seat.centre + dir * (seat.radius + outward * eps);
    Ok(
        ogeom_algo::classify_on_face(model, &seat.cap_face, beside, deflection, tol)?
            == ogeom_algo::Containment::In,
    )
}

/// The flanks every revolved wedge shares: the band of the wall down to the
/// tangency ring, the annulus of the cap out to its own, and the three rings
/// bounding them. Only the face between the two tangency rings differs:
/// a quarter-tube for the fillet, a cone for the chamfer.
pub(crate) struct RevolvedFlanks {
    /// The band of the wall between the rim and `wall_ring`.
    pub wall_band: Shape,
    /// The annulus of the cap between the rim and `cap_ring`.
    pub annulus: Shape,
    /// The tangency ring on the wall, `wall_depth` from the rim.
    pub wall_ring: Shape,
    /// The tangency ring on the cap, at radius `cap_rho`.
    pub cap_ring: Shape,
}

/// Build the flanks: legs `wall_depth` down the wall and in to `cap_rho`
/// along the cap, each coincident with the solid's own face, aligned when
/// subtracting and opposed when fusing, which is what the melt needs.
pub(crate) fn revolved_flanks(
    model: &mut Model,
    seat: &RevolvedSeat,
    wall_depth: f64,
    cap_rho: f64,
    tol: Tolerances,
) -> OgeomResult<RevolvedFlanks> {
    let wall_level = seat.centre - seat.up * (seat.tau * wall_depth);
    let apex_ring = seat.ring(model, seat.centre, seat.radius, tol)?;
    let wall_ring = seat.ring(model, wall_level, seat.radius, tol)?;
    let cap_ring = seat.ring(model, seat.centre, cap_rho, tol)?;

    let wall_band = {
        let origin = seat.centre - seat.up * (wall_depth + 1.0);
        let surface: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(seat.frame_at(origin, tol)?, seat.radius, tol)?,
            (0.0, 2.0 * wall_depth + 2.0),
        )?
        .into();
        let band = make_revolution_band(model, &surface, &wall_ring, &apex_ring, tol)?;
        if seat.sigma * seat.tau < 0.0 {
            band.reversed()
        } else {
            band
        }
    };

    let annulus = {
        let plane = Plane::through(seat.centre, Direction::new(seat.up * seat.tau, tol)?);
        let reach = (seat.radius + wall_depth + cap_rho) * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        // Both rings turn about `up`. The wider is the outer boundary and
        // turns positively about the face's normal, the narrower the hole
        // and turns against it, so the face keeps itself on their left.
        let (wide, narrow) = if cap_rho > seat.radius {
            (&cap_ring, &apex_ring)
        } else {
            (&apex_ring, &cap_ring)
        };
        let with_normal = seat.tau > 0.0;
        let turned = |ring: &Shape, positive: bool| {
            if positive == with_normal {
                ring.clone()
            } else {
                ring.reversed()
            }
        };
        let outer = ogeom_algo::make_wire(model, &[turned(wide, true)], tol)?.shape;
        let inner = ogeom_algo::make_wire(model, &[turned(narrow, false)], tol)?.shape;
        let face = ogeom_algo::make_face(model, surface.clone(), &[outer, inner], tol)?.shape;
        let surface_id = {
            let Some(node) = model.node(&face) else {
                ogeom_bail!(Dangling, "the face just built is not in this model");
            };
            let NodeData::Face(data) = node.data() else {
                ogeom_bail!(Construction, "the face holds no face data");
            };
            data.surface
        };
        for pedge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            let (curve, prange) = edge_curve(model, &pedge, tol)?;
            let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&curve, &surface, tol) else {
                ogeom_bail!(
                    Construction,
                    "an annulus edge has no closed-form pcurve on its plane"
                );
            };
            ogeom_algo::attach_pcurve(
                model,
                &pedge,
                pcurve,
                surface_id,
                ogeom_topo::Location::identity(),
                prange,
            )?;
        }
        face
    };

    Ok(RevolvedFlanks {
        wall_band,
        annulus,
        wall_ring,
        cap_ring,
    })
}

/// Whether `probe` sits strictly inside the face's trim.
pub(crate) fn on_face_side(
    model: &Model,
    face: &Shape,
    probe: ogeom_math::Point,
    deflection: ogeom_mesh::Deflection,
    tol: Tolerances,
) -> OgeomResult<bool> {
    Ok(
        ogeom_algo::classify_on_face(model, face, probe, deflection, tol)?
            == ogeom_algo::Containment::In,
    )
}

/// An edge along the segment from `from` to `to`, parameterized by arc
/// length, joining two existing vertices.
///
/// A wire chains through shared vertex *objects*, not through coincident
/// coordinates, which is why this takes the vertices and not just the
/// points.
pub(crate) fn segment_between(
    model: &mut Model,
    from: (&Shape, Point),
    to: (&Shape, Point),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let line = LineCurve::segment(from.1, to.1, tol)?;
    let curve = Curve::Line(line);
    let domain = curve.domain();
    Ok(make_edge_between(model, curve, domain, from.0, to.0, tol)?.shape)
}

/// A face on `surface` bounded by `edges` in traversal order, with an exact
/// same-parameter pcurve attached to every edge.
///
/// [`ogeom_algo::make_face_with_pcurves`] with one wire: the blend keeps this
/// thin name because every wedge face is a single loop.
pub(crate) fn face_from_edges(
    model: &mut Model,
    surface: SurfaceGeometry,
    edges: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    Ok(ogeom_algo::make_face_with_pcurves(model, surface, &[edges.to_vec()], tol)?.shape)
}

/// Record every face of `built`'s result that no face of `solid` reaches
/// through its history as generated from `edge`: a blend or bevel face
/// comes from the tool the boolean took it from, and the edge it replaces
/// is where a caller naming faces through the operation finds it.
pub(crate) fn credit_new_faces(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    built: &mut ogeom_algo::Built,
) -> OgeomResult<()> {
    // A wedge set aside for later leaves the solid as it was: nothing new.
    if built.shape.is_same(solid) {
        return Ok(());
    }
    let inputs = ogeom_topo::explore_unique(model, solid, ShapeType::Face)?;
    for face in ogeom_topo::explore_unique(model, &built.shape, ShapeType::Face)? {
        let reached = inputs.iter().any(|g| {
            built
                .history
                .modified(g)
                .iter()
                .chain(built.history.generated(g))
                .any(|x| x.is_same(&face))
        });
        if !reached {
            built.history.generate(edge, face);
        }
    }
    Ok(())
}

/// A blend's wedge set aside to be applied with others in one boolean.
pub(crate) struct Wedge {
    /// The wedge as a solid.
    pub(crate) solid: Shape,
    /// Fused on (a concave edge) rather than cut away.
    pub(crate) additive: bool,
    /// The edge it rounds.
    pub(crate) edge: Option<Shape>,
    /// What the blend is refused as when the boolean cannot apply the
    /// wedge, where its builder knows a contact the melt does not resolve.
    pub(crate) melt_refusal: Option<&'static str>,
}

/// Name the refusal the wedge last set aside is refused as, should the
/// boolean fail to apply it. Nothing when no caller is collecting.
pub(crate) fn note_melt_refusal(refusal: &'static str) {
    WEDGES.with(|held| {
        if let Some(last) = held.borrow_mut().as_mut().and_then(|w| w.last_mut()) {
            last.melt_refusal = Some(refusal);
        }
    });
}

thread_local! {
    /// Where [`apply_wedge`] sets its wedge aside, while a caller collects.
    static WEDGES: std::cell::RefCell<Option<Vec<Wedge>>> = const { std::cell::RefCell::new(None) };
}

/// Run `f` with every wedge [`apply_wedge`] builds set aside rather than
/// applied, and hand the wedges back with its result.
pub(crate) fn collecting_wedges<T>(f: impl FnOnce() -> T) -> (T, Vec<Wedge>) {
    struct Restore(Option<Vec<Wedge>>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let outer = self.0.take();
            WEDGES.with(|held| *held.borrow_mut() = outer);
        }
    }
    let outer = WEDGES.with(|held| held.borrow_mut().replace(Vec::new()));
    let restore = Restore(outer);
    let out = f();
    let wedges = WEDGES
        .with(|held| held.borrow_mut().take())
        .unwrap_or_default();
    drop(restore);
    (out, wedges)
}

/// Sew the wedge's faces, demand a closed shell, and apply it to the solid:
/// subtracted on a convex edge, fused on a concave one. Either way the
/// history reads the same truth: the edge the blend replaces is gone.
pub(crate) fn apply_wedge(
    model: &mut Model,
    solid: &Shape,
    edge: Option<&Shape>,
    faces: &[Shape],
    additive: bool,
    tol: Tolerances,
) -> OgeomResult<ogeom_algo::Built> {
    let sewn = ogeom_algo::sew(model, faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the blend wedge did not close");
    }
    let wedge = ogeom_algo::make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;
    let set_aside = WEDGES.with(|held| {
        held.borrow_mut().as_mut().map(|wedges| {
            wedges.push(Wedge {
                solid: wedge.shape.clone(),
                additive,
                edge: edge.cloned(),
                melt_refusal: None,
            });
        })
    });
    if set_aside.is_some() {
        let mut history = ogeom_algo::History::new();
        if let Some(edge) = edge {
            history.delete(edge);
        }
        return Ok(ogeom_algo::Built::new(solid.clone(), history));
    }
    let mut result = if additive {
        ogeom_bool::fuse(model, solid, &wedge.shape, tol)?
    } else {
        ogeom_bool::cut(model, solid, &wedge.shape, tol)?
    };
    if let Some(edge) = edge {
        result.history.delete(edge);
    }
    Ok(result)
}

/// A planar face over explicit corners, with `outward` as its plane normal.
///
/// The corners must be coplanar and `outward` perpendicular to them: the
/// callers know both by construction, which is why this takes the normal
/// instead of rediscovering it.
pub(crate) fn planar_face(
    model: &mut Model,
    corners: &[Point],
    outward: Vector,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let normal = Direction::new(outward, tol)?;
    let plane = Plane::through(corners[0], normal);
    let mut reach = 1.0_f64;
    for p in corners {
        reach = reach.max(p.distance(corners[0]) * 2.0);
    }
    let surface = PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?;
    // The ring keeps the face on its left about `outward`: corners given
    // turning the other way about it are walked in reverse.
    let turn = (0..corners.len()).fold(Vector::ZERO, |sum, i| {
        let (p, q) = (corners[i], corners[(i + 1) % corners.len()]);
        sum + p.to_vector().cross(q.to_vector())
    });
    let mut corners = corners.to_vec();
    if turn.dot(outward) < 0.0 {
        corners.reverse();
    }
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|p| make_vertex(model, *p).shape)
        .collect();
    let mut edges = Vec::with_capacity(corners.len());
    for i in 0..corners.len() {
        let j = (i + 1) % corners.len();
        edges.push(segment_between(
            model,
            (&vertices[i], corners[i]),
            (&vertices[j], corners[j]),
            tol,
        )?);
    }
    face_from_edges(model, surface.into(), &edges, tol)
}

/// How far `face` reaches from the line through `from` along `leg`: the
/// greatest distance any of its edges stands back in that direction.
pub(crate) fn face_reach(
    model: &Model,
    face: &Shape,
    from: ogeom_math::Point,
    leg: ogeom_math::Vector,
    tol: Tolerances,
) -> OgeomResult<f64> {
    use ogeom_geom::Curve3d as _;
    let mut reach = 0.0_f64;
    for edge in ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Edge)? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            continue;
        };
        let placement = edge.transform(model.datums())?;
        for k in 0..=16 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 16.0;
            let p = placement.apply(geometry.point_at(t, tol)?);
            reach = reach.max((p - from).dot(leg));
        }
    }
    Ok(reach)
}

/// Whether a rolling ball of `radius` sits on both faces along `edge`: at
/// stations along it, each face's contact is set back into the face by
/// `radius · tan(φ/2)` (φ the turn between the faces' outward normals,
/// taken on the surfaces where the edge runs), found on the face's surface
/// and classified against the face. A contact outside a face is a band
/// that runs past it; the refusal names the edge and the distance. Edges
/// that are not shared by exactly two faces, and faces meeting tangent,
/// are left to the blend itself. The contacts found come back, each with
/// its face and its setback, for a chain to check its bands against each
/// other.
///
/// `anywhere` asks only that the ball sit somewhere along the edge: a face
/// narrowing toward one end (a pyramid's side under a rounded apex) holds
/// the ball along most of the edge and lets it run off near that end,
/// where the blend stops flush. Only a ball that fits at no station is
/// refused then.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where a contact falls outside its face (at every station, for
/// `anywhere`).
pub(crate) fn ball_fits(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    anywhere: bool,
    tol: Tolerances,
) -> OgeomResult<Vec<Contact>> {
    use ogeom_geom::Surface as _;
    let mut contacts = Vec::new();
    let Ok((curve, range)) = edge_curve(model, edge, tol) else {
        return Ok(contacts);
    };
    // Each face holding the edge, with the edge's sense in that face's walk.
    let mut sides: Vec<(Shape, bool)> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        for e in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            if same_occurrence(model, &e, edge, tol) {
                // The walk as the face presents it: exploring under the
                // face composes its orientation into the edge's already.
                let reversed = e.orientation() == Orientation::Reversed;
                sides.push((face.clone(), reversed));
                break;
            }
        }
    }
    let [(face0, rev0), (face1, rev1)] = sides.as_slice() else {
        return Ok(contacts);
    };
    let normal_at = |face: &Shape, p: Point| -> OgeomResult<Option<(SurfaceGeometry, Vector)>> {
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
            return Ok(None);
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            return Ok(None);
        };
        use ogeom_geom::Transformable as _;
        let placed = surface
            .clone()
            .transformed(&face.transform(model.datums())?, tol)?;
        let found = ogeom_algo::project_on_surface(&placed, p, 16, tol)?;
        let (u, v) = found.parameters;
        let Ok(n) = placed.normal_at(u, v, tol) else {
            return Ok(None);
        };
        let n = if face.orientation() == Orientation::Reversed {
            -n.vector()
        } else {
            n.vector()
        };
        Ok(Some((placed, n)))
    };
    // The edge's curve comes placed. Stations stay clear of its ends,
    // where a neighbouring face in a chain carries the band on.
    let mut missed: Option<String> = None;
    let mut held = false;
    'stations: for fraction in [0.25, 0.5, 0.75] {
        let t = range.0 + (range.1 - range.0) * fraction;
        let p = curve.point_at(t, tol)?;
        let d = curve.d1_at(t, tol)?;
        let m = d.magnitude();
        if m <= tol.confusion() {
            continue;
        }
        let tangent = d / m;
        let (Some((s0, n0)), Some((s1, n1))) = (normal_at(face0, p)?, normal_at(face1, p)?) else {
            continue;
        };
        let turn = n0.dot(n1).clamp(-1.0, 1.0).acos();
        if turn <= tol.angular() * 1e3 {
            continue;
        }
        let setback = radius * (turn * 0.5).tan();
        // Into each face: across the edge, square to the normal, on the
        // side a hair's step from the edge lands inside the face. The walk
        // the orientation flags imply is the first guess; the face itself
        // settles it.
        let into = |face: &Shape,
                    surface: &SurfaceGeometry,
                    n: Vector,
                    reversed: bool|
         -> OgeomResult<Option<Vector>> {
            let walk = if reversed { -tangent } else { tangent };
            let d = n.cross(walk);
            let l = d.magnitude();
            if l <= 0.0 {
                return Ok(None);
            }
            let d = d / l;
            // A hair's step, and a hair whatever the radius: a step scaled
            // by a huge setback lands off the face on both sides.
            let hair = (setback * 1e-3).clamp(tol.confusion() * 10.0, tol.confusion() * 1e5);
            let lands = |step: Vector| -> OgeomResult<bool> {
                let (u, v) = ogeom_algo::project_on_surface(surface, p + step, 16, tol)?.parameters;
                chart_holds(model, face, ogeom_math::Point2::new(u, v), tol)
            };
            let near = [lands(d * hair)?, lands(-d * hair)?];
            // A hair's step lands in the face on both sides where the
            // face's trim turns back close to the edge, inside the chord
            // its chart is drawn with. The whole setback, landing in the
            // face on one side only, settles it there.
            if near == [true, true] {
                let far = [lands(d * setback)?, lands(-d * setback)?];
                if far == [false, true] {
                    return Ok(Some(-d));
                }
            }
            Ok(match near {
                [true, _] => Some(d),
                [false, true] => Some(-d),
                [false, false] => None,
            })
        };
        let (Some(d0), Some(d1)) = (into(face0, &s0, n0, *rev0)?, into(face1, &s1, n1, *rev1)?)
        else {
            continue;
        };
        for (face, surface, step) in [(face0, &s0, d0), (face1, &s1, d1)] {
            let guess = p + step * setback;
            let found = ogeom_algo::project_on_surface(surface, guess, 16, tol)?;
            let (u, v) = found.parameters;
            // On its own face, or on any face of the solid it has run on
            // to (a coplanar piece the face was split from, the next face
            // of a chain): anywhere on the boundary the ball can sit.
            if !chart_holds(model, face, ogeom_math::Point2::new(u, v), tol)?
                && !on_boundary(model, solid, found.point, tol)?
            {
                let refusal = format!(
                    "a blend of radius {radius} on the edge through {p:?} sets back {setback} \
                     across a face that does not reach so far"
                );
                if !anywhere {
                    ogeom_bail!(Construction, "{refusal}");
                }
                missed = Some(refusal);
                continue 'stations;
            }
            contacts.push(Contact {
                face: face.node(),
                at: found.point,
                setback,
            });
        }
        held = true;
    }
    if let Some(refusal) = missed
        && !held
    {
        ogeom_bail!(Construction, "{refusal}");
    }
    Ok(contacts)
}

/// Where a rolling ball touches one of an edge's faces.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Contact {
    /// The face touched.
    pub face: ogeom_topo::TShapeId,
    /// The point of contact.
    pub at: Point,
    /// How far back from the edge it stands.
    pub setback: f64,
}

/// Refuse a chain whose bands overlap across a face they share.
///
/// Each edge's band covers the strip of a face between the edge and its
/// contact line. Two edges of one face that do not meet at a vertex (the
/// two long edges of a narrow strip, a wall's two rims) each need their own
/// strip; where one's contact stands inside the other's strip, the two
/// balls overlap and no pair of blends is the rounding. Edges meeting at a
/// vertex are corner mates, whose bands meet there by design.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where two bands overlap.
pub(crate) fn bands_clear(
    model: &Model,
    edges: &[Shape],
    contacts: &[Vec<Contact>],
    tol: Tolerances,
) -> OgeomResult<()> {
    let vertices = |edge: &Shape| -> OgeomResult<Vec<ogeom_topo::TShapeId>> {
        Ok(ogeom_algo::edge_vertices(model, edge)?
            .map(|(a, b)| vec![a.node(), b.node()])
            .unwrap_or_default())
    };
    // Each edge's box: a contact further from it than a strip is further
    // from the edge too, and the scan along the edge is not needed.
    let mut boxes = Vec::with_capacity(edges.len());
    for edge in edges {
        boxes.push(ogeom_algo::shape_bounds(model, edge, tol)?);
    }
    for (i, first) in edges.iter().enumerate() {
        for (j, second) in edges.iter().enumerate() {
            if i == j {
                continue;
            }
            let (ends, others) = (vertices(first)?, vertices(second)?);
            if ends.iter().any(|v| others.contains(v)) {
                continue;
            }
            let Ok((curve, range)) = edge_curve(model, second, tol) else {
                continue;
            };
            for contact in &contacts[i] {
                // The second edge's own setback on this face, if it has one.
                let Some(strip) = contacts[j]
                    .iter()
                    .filter(|c| c.face == contact.face)
                    .map(|c| c.setback)
                    .reduce(f64::max)
                else {
                    continue;
                };
                if boxes[j].distance_to(contact.at) >= strip {
                    continue;
                }
                let away = distance_to_edge(&curve, range, contact.at, tol)?;
                if away < strip - tol.confusion() * 1e3 {
                    ogeom_bail!(
                        Construction,
                        "two blends of the chain overlap across a face they share: a \
                         contact stands {away} from the other edge, inside its band of {strip}"
                    );
                }
            }
        }
    }
    Ok(())
}

/// The distance from `p` to a curve over `range`: the nearest of a scan,
/// refined in the interval around it.
fn distance_to_edge(
    curve: &Curve,
    range: (f64, f64),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<f64> {
    use ogeom_geom::Curve3d as _;
    const SCAN: u32 = 64;
    let at = |t: f64| -> OgeomResult<f64> { Ok(curve.point_at(t, tol)?.distance(p)) };
    let step = (range.1 - range.0) / f64::from(SCAN);
    let mut best = (range.0, at(range.0)?);
    for k in 1..=SCAN {
        let t = range.0 + step * f64::from(k);
        let d = at(t)?;
        if d < best.1 {
            best = (t, d);
        }
    }
    let (mut lo, mut hi) = ((best.0 - step).max(range.0), (best.0 + step).min(range.1));
    for _ in 0..60 {
        let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if at(a)? < at(b)? {
            hi = b;
        } else {
            lo = a;
        }
    }
    Ok(at(f64::midpoint(lo, hi))?.min(best.1))
}

thread_local! {
    /// Face meshes kept while a caller checks many edges against one
    /// unchanging solid: a face with a hundred holes is otherwise meshed
    /// once per edge on it.
    static FACE_MESHES: std::cell::RefCell<Option<std::collections::HashMap<ogeom_topo::SameKey, std::rc::Rc<ogeom_topo::Triangulation>>>> =
        const { std::cell::RefCell::new(None) };
}

/// Run `f` with face meshes kept from one call of [`chart_holds`] to the
/// next. The solid must not change while it runs.
pub(crate) fn keeping_face_meshes<T>(f: impl FnOnce() -> T) -> T {
    struct Restore(
        Option<
            std::collections::HashMap<ogeom_topo::SameKey, std::rc::Rc<ogeom_topo::Triangulation>>,
        >,
    );
    impl Drop for Restore {
        fn drop(&mut self) {
            let outer = self.0.take();
            FACE_MESHES.with(|held| *held.borrow_mut() = outer);
        }
    }
    let outer =
        FACE_MESHES.with(|held| held.borrow_mut().replace(std::collections::HashMap::new()));
    let _restore = Restore(outer);
    f()
}

/// A face's mesh at the default deflection, kept where a caller asked.
fn face_mesh(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<std::rc::Rc<ogeom_topo::Triangulation>> {
    let key = ogeom_topo::SameKey(face.clone());
    if let Some(mesh) =
        FACE_MESHES.with(|held| held.borrow().as_ref().and_then(|m| m.get(&key).cloned()))
    {
        return Ok(mesh);
    }
    let mesh = std::rc::Rc::new(ogeom_mesh::triangulate_face(
        model,
        face,
        ogeom_mesh::Deflection::default(),
        tol,
    )?);
    FACE_MESHES.with(|held| {
        if let Some(meshes) = held.borrow_mut().as_mut() {
            meshes.insert(key, std::rc::Rc::clone(&mesh));
        }
    });
    Ok(mesh)
}

/// Whether a chart point lies in a face's region: inside one of the
/// triangles its own mesh lays over its chart, or on one's side.
fn chart_holds(
    model: &Model,
    face: &Shape,
    at: ogeom_math::Point2,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let mesh = face_mesh(model, face, tol)?;
    let chart = |i: u32| {
        let (u, v) = mesh.parameters[i as usize];
        ogeom_math::Point2::new(u, v)
    };
    let span = mesh
        .parameters
        .iter()
        .fold(0.0_f64, |m, (u, v)| m.max(u.abs()).max(v.abs()))
        .max(1.0);
    let eps = span * 1e-9;
    Ok(mesh.triangles.iter().any(|t| {
        let [a, b, c] = t.map(chart);
        let cross = |p: ogeom_math::Point2, q: ogeom_math::Point2| {
            (q.x - p.x) * (at.y - p.y) - (q.y - p.y) * (at.x - p.x)
        };
        let (d1, d2, d3) = (cross(a, b), cross(b, c), cross(c, a));
        let neg = d1 < -eps || d2 < -eps || d3 < -eps;
        let pos = d1 > eps || d2 > eps || d3 > eps;
        !(neg && pos)
    }))
}

/// Whether a point lies on some face of `solid`: on its surface, within
/// the model's reach, and inside its chart region.
fn on_boundary(model: &Model, solid: &Shape, point: Point, tol: Tolerances) -> OgeomResult<bool> {
    use ogeom_geom::Transformable as _;
    let reach = tol.confusion() * 100.0;
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            continue;
        };
        let placed = surface
            .clone()
            .transformed(&face.transform(model.datums())?, tol)?;
        let found = ogeom_algo::project_on_surface(&placed, point, 16, tol)?;
        if found.distance > reach {
            continue;
        }
        let (u, v) = found.parameters;
        if chart_holds(model, &face, ogeom_math::Point2::new(u, v), tol)? {
            return Ok(true);
        }
    }
    Ok(false)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use std::collections::HashMap;

    use ogeom_core::Tolerances;
    use ogeom_math::{Frame, Point};
    use ogeom_topo::{Filter, Model, Orientation, Shape, ShapeType, explore, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    /// How many edges two faces of `shape` walk the same way.
    fn walked_one_way(model: &Model, shape: &Shape) -> usize {
        let mut walks: HashMap<_, Vec<Orientation>> = HashMap::new();
        for face in explore(model, shape, Filter::OfType(ShapeType::Face)).unwrap() {
            for wire in model.children_of(&face).unwrap() {
                for edge in model.children_of(&wire).unwrap() {
                    if !model
                        .node(&edge)
                        .unwrap()
                        .data()
                        .as_edge()
                        .unwrap()
                        .degenerate
                    {
                        walks
                            .entry(edge.node())
                            .or_default()
                            .push(edge.orientation());
                    }
                }
            }
        }
        walks
            .values()
            .filter(|w| w.len() == 2 && w[0] == w[1])
            .count()
    }

    /// The edge of `shape` whose middle stands nearest `at`.
    fn edge_at(model: &Model, shape: &Shape, at: Point) -> Shape {
        let off = |e: &Shape| {
            let b = ogeom_algo::shape_bounds(model, e, T).unwrap();
            b.low().unwrap().midpoint(b.high().unwrap()).distance(at)
        };
        explore_unique(model, shape, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .min_by(|a, b| off(a).total_cmp(&off(b)))
            .unwrap()
    }

    /// The wedges a straight edge's fillet and chamfer and a drum rim's
    /// fillet take away walk each edge between two of their faces once each
    /// way, and hold the volume the blend removes: a quarter round's
    /// section, or a half square's, along the edge or round the rim at the
    /// section's centroid.
    #[test]
    fn blend_wedges_walk_each_edge_once_each_way() {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T)
            .unwrap()
            .shape;
        let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 4.0, T)
            .unwrap()
            .shape;
        let ridge = edge_at(&model, &block, Point::new(5.0, 0.0, 4.0));
        let rim = edge_at(&model, &drum, Point::new(0.0, 0.0, 4.0));
        let pi = core::f64::consts::PI;
        let section = 1.0 - pi / 4.0;
        let centroid = (5.0 / 6.0 - pi / 4.0) / section;
        let cases = [
            (&block, &ridge, false, section * 10.0),
            (&block, &ridge, true, 0.5 * 10.0),
            (&drum, &rim, false, 2.0 * pi * (3.0 - centroid) * section),
        ];
        for (solid, edge, bevel, removed) in cases {
            let (made, wedges) = super::collecting_wedges(|| {
                if bevel {
                    crate::chamfer_edge(&mut model, solid, edge, 1.0, T)
                } else {
                    crate::fillet_edge(&mut model, solid, edge, 1.0, T)
                }
            });
            made.unwrap();
            let [wedge] = &wedges[..] else {
                panic!("one wedge, not {}", wedges.len());
            };
            assert_eq!(walked_one_way(&model, &wedge.solid), 0, "bevel {bevel}");
            let volume = ogeom_algo::volume_properties(
                &model,
                &wedge.solid,
                ogeom_mesh::Deflection::default(),
                T,
            )
            .unwrap()
            .mass;
            assert!(
                (volume - removed).abs() < removed * 1e-6,
                "bevel {bevel}: {volume} against {removed}"
            );
        }
    }

    /// A stadium prism from `z0` up `height`: two runs of `length` along
    /// x joined by semicircles of radius `r`, centred on (0, 0) and
    /// (`length`, 0).
    fn stadium(model: &mut Model, length: f64, r: f64, z0: f64, height: f64) -> Shape {
        use ogeom_math::{Circle, Direction, Plane, Vector};
        let corner = |x: f64, y: f64| Point::new(x, y, z0);
        let [tl, tr, br, bl] = [(0.0, r), (length, r), (length, -r), (0.0, -r)]
            .map(|(x, y)| ogeom_algo::make_vertex(model, corner(x, y)).shape);
        let arc = |model: &mut Model, cx: f64, x: Direction, from: &Shape, to: &Shape| {
            let frame = Frame::new(corner(cx, 0.0), Direction::Z, x, T).unwrap();
            let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(
                Circle::new(frame, r, T).unwrap(),
            ));
            ogeom_algo::make_edge_between(model, curve, (0.0, core::f64::consts::PI), from, to, T)
                .unwrap()
                .shape
        };
        let top = super::segment_between(model, (&tl, corner(0.0, r)), (&tr, corner(length, r)), T)
            .unwrap();
        let right = arc(model, length, -Direction::Y, &br, &tr);
        let bottom =
            super::segment_between(model, (&br, corner(length, -r)), (&bl, corner(0.0, -r)), T)
                .unwrap();
        let left = arc(model, 0.0, Direction::Y, &tl, &bl);
        let plane = ogeom_geom::PlaneSurface::over(
            Plane::through(corner(0.0, 0.0), Direction::Z),
            (-100.0, 100.0),
            (-100.0, 100.0),
        )
        .unwrap();
        let face = ogeom_algo::make_face_with_pcurves(
            model,
            plane.into(),
            &[vec![left, bottom.reversed(), right, top.reversed()]],
            T,
        )
        .unwrap()
        .shape;
        ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, height), T)
            .unwrap()
            .shape
    }

    /// The wedges an arc of a rim takes away or fills walk each edge between
    /// two of their faces once each way, whichever side of the rim the cap
    /// and the material stand on, and hold the volume the blend removes or
    /// adds: a quarter round's section swept round the half turn at its
    /// centroid's radius, inside the rim where the wall faces out and
    /// outside it where the wall faces the axis or the rim is concave.
    #[test]
    fn rim_arc_wedges_walk_each_edge_once_each_way() {
        let mut model = Model::new();
        let (length, r) = (10.0, 5.0);
        let post = stadium(&mut model, length, r, 0.0, 4.0);
        let block = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(-10.0, -10.0, 0.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            (30.0, 20.0, 4.0),
            T,
        )
        .unwrap()
        .shape;
        let bore = stadium(&mut model, length, r, -1.0, 6.0);
        let slot = ogeom_bool::cut(&mut model, &block, &bore, T).unwrap().shape;
        let plate = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(-10.0, -10.0, -4.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            (30.0, 20.0, 4.0),
            T,
        )
        .unwrap()
        .shape;
        let boss = ogeom_bool::fuse(&mut model, &plate, &post, T)
            .unwrap()
            .shape;
        let pi = core::f64::consts::PI;
        let section = 1.0 - pi / 4.0;
        let centroid = (10.0 - 3.0 * pi) / (12.0 - 3.0 * pi);
        let inside = pi * (r - centroid) * section;
        let outside = pi * (r + centroid) * section;
        let at = |z: f64| Point::new(length + r / 2.0, 0.0, z);
        let cases = [
            ("post top", &post, at(4.0), inside),
            ("post foot", &post, at(0.0), inside),
            ("slot top", &slot, at(4.0), outside),
            ("slot foot", &slot, at(0.0), outside),
            ("boss foot", &boss, at(0.0), outside),
        ];
        for (name, solid, near, expected) in cases {
            let rim = edge_at(&model, solid, near);
            let (made, wedges) =
                super::collecting_wedges(|| crate::fillet_edge(&mut model, solid, &rim, 1.0, T));
            made.unwrap();
            let [wedge] = &wedges[..] else {
                panic!("{name}: one wedge, not {}", wedges.len());
            };
            assert_eq!(walked_one_way(&model, &wedge.solid), 0, "{name}");
            let volume = ogeom_algo::volume_properties(
                &model,
                &wedge.solid,
                ogeom_mesh::Deflection::default(),
                T,
            )
            .unwrap()
            .mass;
            assert!(
                (volume - expected).abs() < expected * 1e-6,
                "{name}: {volume} against {expected}"
            );
        }
    }

    /// The edge of `shape` whose curve is the first `kind` accepts.
    fn edge_on(model: &Model, shape: &Shape, kind: impl Fn(&ogeom_geom::Curve) -> bool) -> Shape {
        explore_unique(model, shape, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| super::edge_curve(model, e, T).is_ok_and(|(curve, _)| kind(&curve)))
            .unwrap()
    }

    /// The volume of a solid, through the mesh where a face is fitted.
    fn volume(model: &Model, solid: &Shape) -> f64 {
        ogeom_algo::volume_properties(model, solid, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass
    }

    /// The wedges a marched band takes away or fills walk each edge between
    /// two of their faces once each way and pass the checker, as do the
    /// blended solids, and each holds
    /// what the blend adds to or removes from the solid: a closed band round
    /// a branch's foot, a band running out through a block's walls at both
    /// ends, and a band pinched where two equal drums turn tangent.
    #[test]
    #[ignore = "heavy"]
    fn marched_wedges_walk_each_edge_once_each_way() {
        use ogeom_geom::Curve;
        use ogeom_math::{Direction, Vector};
        let mut model = Model::new();
        let across = |x: f64, y: f64, z: f64| {
            Frame::new(Point::new(x, y, z), Direction::X, Direction::Y, T).unwrap()
        };
        let main = ogeom_algo::make_cylinder(&mut model, across(-20.0, 0.0, 0.0), 10.0, 40.0, T)
            .unwrap()
            .shape;
        let branch = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
            .unwrap()
            .shape;
        let branched = ogeom_bool::fuse(&mut model, &main, &branch, T)
            .unwrap()
            .shape;
        let foot = edge_on(&model, &branched, |c| {
            use ogeom_geom::Curve3d as _;
            let (lo, hi) = c.domain();
            matches!(c, Curve::BSpline(_))
                && c.point_at(lo, T)
                    .and_then(|p| c.point_at(hi, T).map(|q| p.distance(q)))
                    .is_ok_and(|d| d <= 1e-6)
        });

        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let tilt = 0.35_f64;
        let axis = Direction::new(Vector::new(0.0, tilt.cos(), -tilt.sin()), T).unwrap();
        let drum = ogeom_algo::make_cylinder(
            &mut model,
            Frame::new(Point::new(10.0, -5.0, 11.5), axis, Direction::X, T).unwrap(),
            4.0,
            30.0,
            T,
        )
        .unwrap()
        .shape;
        let grooved = ogeom_bool::cut(&mut model, &block, &drum, T).unwrap().shape;
        let groove = edge_at(&model, &grooved, Point::new(10.0, 14.84, 0.0));

        let upright = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
            .unwrap()
            .shape;
        let lying = ogeom_algo::make_cylinder(&mut model, across(-20.0, 0.0, 10.0), 5.0, 40.0, T)
            .unwrap()
            .shape;
        let crossed = ogeom_bool::fuse(&mut model, &upright, &lying, T)
            .unwrap()
            .shape;
        let pinch = edge_on(&model, &crossed, |c| matches!(c, Curve::Ellipse(_)));

        let cases = [
            ("closed band", &branched, &foot, 1.5),
            ("run-out band", &grooved, &groove, 1.0),
            ("pinched band", &crossed, &pinch, 1.0),
        ];
        for (name, solid, edge, radius) in cases {
            let (made, wedges) =
                super::collecting_wedges(|| crate::fillet_edge(&mut model, solid, edge, radius, T));
            made.unwrap();
            let [wedge] = &wedges[..] else {
                panic!("{name}: one wedge, not {}", wedges.len());
            };
            assert_eq!(walked_one_way(&model, &wedge.solid), 0, "{name}");
            // The checker asks the walk too, and that every vertex is held
            // as loosely as the edges meeting it.
            let found = ogeom_algo::check(&model, &wedge.solid, T).unwrap();
            assert!(found.is_usable(), "{name}: {found}");
            // What the blend adds is the wedge; what it removes is the part
            // of the wedge inside the solid, the band running on past the
            // walls where it runs out.
            let blended = crate::fillet_edge(&mut model, solid, edge, radius, T)
                .unwrap()
                .shape;
            let found = ogeom_algo::check(&model, &blended, T).unwrap();
            assert!(found.is_usable(), "{name}: the blended solid: {found}");
            let moved = (volume(&model, &blended) - volume(&model, solid)).abs();
            let held = if wedge.additive {
                volume(&model, &wedge.solid)
            } else {
                let inside = ogeom_bool::common(&mut model, &wedge.solid, solid, T)
                    .unwrap()
                    .shape;
                volume(&model, &inside)
            };
            assert!(
                (held - moved).abs() < moved * 1e-3,
                "{name}: the wedge holds {held}, the blend moves {moved}"
            );
        }
    }
}
