//! Shared scaffolding for the subtractive blends.
//!
//! A chamfer and a constant-radius fillet on a straight edge between planar
//! faces differ only in the face that replaces the edge — a bevel plane for
//! one, a tangent cylinder for the other. Everything around that face — the
//! seat on the solid, the legs running along the adjacent faces, faces
//! assembled from explicit curves with exact pcurves — is one piece of
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
    /// The two faces themselves, in the same order — so a caller naming a
    /// face can find which leg it owns.
    pub faces: [Shape; 2],
    /// Whether the edge is convex: material inside the dihedral, so a blend
    /// subtracts. A concave edge's blend adds, with every sign mirrored.
    pub convex: bool,
}

impl Seat {
    /// On face `i`, the unit direction perpendicular to the edge that walks
    /// away from the *other* face — into the material the blend cuts back
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

/// An edge's 3D curve and range, cloned out of the model.
pub(crate) fn edge_curve(model: &Model, edge: &Shape) -> OgeomResult<(Curve, (f64, f64))> {
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
    Ok((geometry.clone(), *range))
}

/// Find the seat of a blend: the straight edge's ends and direction, and the
/// outward normals of the exactly two planar faces of `solid` meeting there.
///
/// Refuses concave and tangent edges: the wedge these blends subtract lies in
/// the material only when the edge is convex. A concave blend *adds* material
/// and is a different construction — docs/PARITY.md, fillet.edge-blends.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the edge is
/// not straight, is not shared by exactly two planar faces of `solid`, or the
/// edge is concave or tangent.
pub(crate) fn planar_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Seat> {
    let (curve, range) = edge_curve(model, edge)?;
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
            .any(|e| e.node() == edge.node());
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
    // Which way the first face actually extends from the edge — sampled
    // against its own trim — leans behind the other face's plane on a convex
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
    /// The rim's radius — also the wall's.
    pub radius: f64,
    /// The cap's outward unit normal.
    pub up: Vector,
    /// The rim frame's `x`, so every ring built on the seat shares a
    /// parameter origin.
    pub x_ref: Direction,
    /// The wall's outward radial sign: `+1` material inside, `-1` a bore.
    pub sigma: f64,
    /// `+1` when the wall extends away from the cap's outward side — the rim
    /// configurations — and `-1` alongside it, the concave seats.
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
            .any(|e| e.node() == edge.node());
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

    // Which side of the cap the wall extends: read from the wall band's other
    // ring. `tau` positive means away from the cap's outward side — the rim
    // configurations — and negative means alongside it, the concave seats.
    let tau = {
        let mut side = None;
        for e in explore(model, &wall_face, Filter::OfType(ShapeType::Edge))? {
            if e.node() == edge.node() {
                continue;
            }
            // Any circular edge of the wall at another height says which
            // side the wall extends — whether the ring survives as one
            // closed edge or as the arcs a boolean rebuilt it into.
            let Ok((curve, _)) = edge_curve(model, &e) else {
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

/// The flanks every revolved wedge shares: the band of the wall down to the
/// tangency ring, the annulus of the cap out to its own, and the three rings
/// bounding them. Only the face between the two tangency rings differs —
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
/// along the cap, each coincident with the solid's own face — aligned when
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
        let outer = ogeom_algo::make_wire(model, std::slice::from_ref(&apex_ring), tol)?.shape;
        let inner = ogeom_algo::make_wire(model, std::slice::from_ref(&cap_ring), tol)?.shape;
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
            let (curve, prange) = edge_curve(model, &pedge)?;
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
/// coordinates — which is why this takes the vertices and not just the
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

/// Sew the wedge's faces, demand a closed shell, and apply it to the solid —
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
/// The corners must be coplanar and `outward` perpendicular to them — the
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
