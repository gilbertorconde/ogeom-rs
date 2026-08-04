//! Where a point sits relative to a shape.
//!
//! Three answers, never two: inside, outside, or *on* the boundary within
//! tolerance. The third is not a hedge. Geometry that meets is the normal case
//! in a kernel — a boolean's whole job is finding it — and a classifier that
//! forces every point to one side has to pick, silently, for exactly the points
//! where the choice matters most.
//!
//! # Accuracy
//!
//! [`classify_in_solid`] and [`classify_on_face`] work from the tessellation,
//! so a point nearer the boundary than the deflection cannot be told from one
//! on it. That is reported as [`Containment::On`] rather than guessed: the
//! band the answer is uncertain within is the deflection, and saying so is the
//! difference between an approximate answer and a wrong one.
//!
//! Tightening the deflection narrows the band. It never removes it — the
//! exact question needs ray/surface intersection, which is
//! [`classify_in_solid_exact`]: rays cast against the faces' true surfaces,
//! where the uncertain band shrinks from the deflection to the tolerance.

use og_core::{OgResult, Tolerances, og_bail};
use og_math::{Direction, Point, Point2, Vector};
use og_mesh::{Deflection, face_boundary, inside_boundary, triangulate};
use og_topo::{Model, NodeData, Shape, ShapeType};

use crate::measure::project_on_surface;

/// Where a point sits relative to a shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Containment {
    /// Strictly inside.
    In,
    /// On the boundary, within tolerance of it.
    On,
    /// Strictly outside.
    Out,
}

impl Containment {
    /// Whether the point is inside or on the boundary.
    #[must_use]
    pub const fn is_inside_or_on(self) -> bool {
        matches!(self, Self::In | Self::On)
    }

    /// The classification of the same point against the complement.
    ///
    /// `In` and `Out` swap; `On` is its own opposite, since a boundary is
    /// shared by both sides.
    #[must_use]
    pub const fn inverted(self) -> Self {
        match self {
            Self::In => Self::Out,
            Self::On => Self::On,
            Self::Out => Self::In,
        }
    }
}

/// Where a point sits relative to a face.
///
/// A point off the face's surface is [`Containment::Out`]: a face is a patch of
/// surface, so "inside" can only mean inside its trimming, and a point in space
/// that does not lie on the surface at all is not inside anything.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `face` is not a
/// face, or the deflection settings are unusable;
/// [`OgError::Dangling`](og_core::OgError::Dangling) if a handle fails to
/// resolve.
pub fn classify_on_face(
    model: &Model,
    face: &Shape,
    point: Point,
    deflection: Deflection,
    tol: Tolerances,
) -> OgResult<Containment> {
    deflection.validate()?;
    if model.kind_of(face)? != ShapeType::Face {
        og_bail!(Construction, "expected a face");
    }
    let Some(node) = model.node(face) else {
        og_bail!(Dangling, "face is not in this model");
    };
    let NodeData::Face(data) = node.data() else {
        og_bail!(Construction, "face node holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        og_bail!(Dangling, "face refers to a surface not in this model");
    };

    // Into the surface's own frame first: the trimming lives in parameter
    // space, and the face may be placed anywhere.
    let placement = face.transform(model.datums())?;
    let local = placement.inverse()?.apply(point);

    // A grid dense enough to bracket a foot point on a surface that folds:
    // too coarse and Newton starts in the wrong basin and converges on a far
    // side of a cylinder.
    let projection = project_on_surface(surface, local, 32, tol)?;
    let reach = tol.confusion().max(data.tolerance.get());
    if projection.distance > reach {
        return Ok(Containment::Out);
    }

    let rings = face_boundary(model, face, deflection, tol)?;
    let (u, v) = projection.parameters;
    let at = Point2::new(u, v);

    // The uncertain band, converted from a distance in space into one in
    // parameter units through the surface's own scale. A fixed parameter
    // tolerance would be metres wide at a sphere's equator and nothing at its
    // pole.
    let band = parametric_band(surface, (u, v), reach, tol);
    if distance_to_rings(&rings, at) <= band {
        return Ok(Containment::On);
    }
    Ok(if inside_boundary(&rings, at) {
        Containment::In
    } else {
        Containment::Out
    })
}

/// Where a point sits relative to a closed shell or solid.
///
/// Ray casting against the tessellation: a ray from the point crosses the
/// boundary an odd number of times if and only if it started inside.
///
/// # Errors
///
/// As [`triangulate()`], plus
/// [`OgError::Construction`](og_core::OgError::Construction) if the boundary is
/// not closed — an open shell has no inside — and
/// [`OgError::NotDone`](og_core::OgError::NotDone) if every ray tried hit an
/// edge or a vertex, where the crossing count is ambiguous.
pub fn classify_in_solid(
    model: &Model,
    solid: &Shape,
    point: Point,
    deflection: Deflection,
    tol: Tolerances,
) -> OgResult<Containment> {
    deflection.validate()?;
    let mesh = triangulate(model, solid, deflection, tol)?;
    if mesh.is_empty() || !mesh.is_closed() {
        og_bail!(
            Construction,
            "the boundary is not closed, so there is no inside to be in"
        );
    }

    let triangles: Vec<[Point; 3]> = mesh
        .triangles
        .iter()
        .map(|t| t.map(|i| mesh.positions[i as usize]))
        .collect();

    // On the boundary beats either side, and is decided in space rather than
    // along a ray: a point sitting on a face is on the boundary from every
    // direction, and no crossing count says so.
    let reach = tol.confusion() + deflection.chord;
    for t in &triangles {
        if distance_to_triangle(point, *t) <= reach {
            return Ok(Containment::On);
        }
    }

    // A ray that grazes an edge or passes through a vertex is counted once by
    // one triangle and twice by its neighbour, or not at all. Rather than
    // patch the count, notice the near-miss and cast again somewhere else.
    for direction in RAY_DIRECTIONS {
        let ray = Direction::new(Vector::new(direction[0], direction[1], direction[2]), tol)?;
        if let Some(crossings) = count_crossings(&triangles, point, ray, tol) {
            return Ok(if crossings % 2 == 1 {
                Containment::In
            } else {
                Containment::Out
            });
        }
    }
    og_bail!(
        NotDone,
        "every ray tried met an edge or a vertex, where the crossing count is \
         ambiguous"
    )
}

/// Where a point sits relative to a closed shell or solid, decided against
/// the true surfaces.
///
/// This is the promise the tessellated classifier's documentation makes on
/// behalf of "a later layer", kept: rays are cast against each face's actual
/// geometry through the curve/surface intersector, so the band where the
/// answer is *On* rather than a side is the tolerance, not the deflection. A
/// point a micron off a sphere's wall classifies as the side it is on;
/// [`classify_in_solid`] at any practical deflection could only say *On*.
///
/// The crossing-parity argument is the same as the tessellated one, and so is
/// the discipline about degenerate hits. A ray that grazes a surface
/// tangentially, meets a face too near its boundary to be sure which side of
/// the trim it crossed, lies *in* a face's surface, or passes through a
/// pole or an apex, is not patched into a count — the ray is abandoned and
/// the next direction tried. Six directions, deterministic, none axis-aligned.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the shape is
/// not a solid or a shell, or its boundary is not closed;
/// [`OgError::NotDone`](og_core::OgError::NotDone) if every ray met a
/// degeneracy — which a handful of deliberately skew directions makes an
/// engineered case rather than an encountered one.
pub fn classify_in_solid_exact(
    model: &Model,
    solid: &Shape,
    point: Point,
    tol: Tolerances,
) -> OgResult<Containment> {
    let kind = model.kind_of(solid)?;
    if !matches!(kind, ShapeType::Solid | ShapeType::Shell) {
        og_bail!(Construction, "expected a solid or a shell, got {kind:?}");
    }
    let shells = if kind == ShapeType::Shell {
        vec![solid.clone()]
    } else {
        og_topo::explore_unique(model, solid, ShapeType::Shell)?
    };
    if shells.is_empty() {
        og_bail!(Construction, "the shape has no shell, so no boundary");
    }
    for shell in &shells {
        if !crate::build::is_shell_closed(model, shell)? {
            og_bail!(
                Construction,
                "the boundary is not closed, so there is no inside to be in"
            );
        }
    }

    // Anything outside the shape's bound is outside the shape, and the bound
    // also sets how long a ray must be to have left everything behind.
    let bound = crate::measure::shape_bounds(model, solid, tol)?;
    let reach = tol.confusion();
    let (Some(centre), diagonal) = (bound.centre(), bound.diagonal()) else {
        og_bail!(Construction, "the boundary bounds nothing");
    };
    if !bound.expanded(reach).contains(point) {
        return Ok(Containment::Out);
    }
    let length = point.distance(centre) + diagonal + 1.0;

    // The rings' own polylining error, spatially: they are only used to
    // decide which side of a face's trim a crossing landed, and a crossing
    // nearer the ring than this is ambiguous rather than decided.
    let ring_chord = tol.confusion() * 1e4;
    let ring_deflection = Deflection {
        chord: ring_chord,
        angular: 0.05,
        ..Deflection::default()
    };

    let faces = og_topo::explore_unique(model, solid, ShapeType::Face)?;
    let mut prepared = Vec::with_capacity(faces.len());
    for face in faces {
        // On the boundary beats either side, and each face answers exactly:
        // projection distance against the true surface, trimming in parameter
        // space.
        if classify_on_face(model, &face, point, ring_deflection, tol)? != Containment::Out {
            return Ok(Containment::On);
        }
        let Some(node) = model.node(&face) else {
            og_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            og_bail!(Dangling, "face refers to a surface not in this model");
        };
        let inverse = face.transform(model.datums())?.inverse()?;
        let rings = face_boundary(model, &face, ring_deflection, tol)?;
        prepared.push((surface.clone(), inverse, rings));
    }

    'directions: for direction in RAY_DIRECTIONS {
        let along = Vector::new(direction[0], direction[1], direction[2]);
        let far = point + along * length;
        let mut crossings = 0_usize;

        for (surface, inverse, rings) in &prepared {
            // Into the face's frame, as two points rather than a direction, so
            // a placement that scales still carries the ray faithfully.
            let from = inverse.apply(point);
            let to = inverse.apply(far);
            let ray: og_geom::Curve = og_geom::LineCurve::segment(from, to, tol)?.into();
            let found = og_intersect::intersect_curve_surface(
                &ray,
                surface,
                og_intersect::CurveSurfaceOptions::default(),
                tol,
            )?;
            if !found.lying.is_empty() {
                // The ray runs in this face's surface: it crosses nothing and
                // touches everything, which no parity expresses.
                continue 'directions;
            }
            for hit in &found.crossings {
                if hit.on_curve <= tol.confusion() {
                    // At the very start. The boundary test above said the
                    // point is off every face, but a crossing this close is
                    // still too near the start to trust its side.
                    continue 'directions;
                }
                let (u, v) = hit.on_surface;
                use og_geom::Surface as _;
                let Ok((du, dv)) = surface.d1_at(u, v, tol) else {
                    continue 'directions;
                };
                let normal = du.cross(dv);
                if normal.magnitude() <= tol.confusion() {
                    // A pole or an apex: no normal, no transversality.
                    continue 'directions;
                }
                let ray_direction = (to - from) / (to - from).magnitude();
                if normal.dot(ray_direction).abs() <= GRAZING * normal.magnitude() {
                    // Tangential. A grazing contact counts once where parity
                    // needs zero or two; abandon the ray rather than guess.
                    continue 'directions;
                }
                let at = Point2::new(u, v);
                let band = parametric_band(surface, (u, v), reach + ring_chord, tol);
                if distance_to_rings(rings, at) <= band {
                    // Too near the face's boundary to know which side of the
                    // trim it crossed — and a shared edge would be counted by
                    // both faces or neither.
                    continue 'directions;
                }
                if inside_boundary(rings, at) {
                    crossings += 1;
                }
            }
        }
        return Ok(if crossings % 2 == 1 {
            Containment::In
        } else {
            Containment::Out
        });
    }
    og_bail!(
        NotDone,
        "every ray tried met a tangency, a boundary, or a degenerate point, \
         where the crossing count is ambiguous"
    )
}

/// The sine of the shallowest crossing angle a counted ray/surface crossing
/// may make.
///
/// Below this the hit is treated as a graze: a tangential contact is one
/// crossing where parity arithmetic needs zero or two, and the seed of the
/// threshold is the same as the marching intersector's `SHALLOWEST` — beneath
/// a microradian, rounding in the evaluated normal can no longer tell a
/// crossing from a touch.
const GRAZING: f64 = 1e-6;

/// Directions to cast rays along, tried in order.
///
/// Deterministic, not random: a classifier that gives different answers on
/// different runs is worse than one that fails, because the failure can be
/// handled and the inconsistency cannot. They are deliberately not axis-aligned
/// and share no common plane, so a mesh built on a regular grid — where an
/// axis-aligned ray runs along a whole row of edges — does not defeat all of
/// them at once.
const RAY_DIRECTIONS: [[f64; 3]; 6] = [
    [0.577_35, 0.577_35, 0.577_35],
    [-0.301_5, 0.904_5, 0.301_5],
    [0.727_6, -0.485_1, 0.485_1],
    [0.259_5, 0.259_5, -0.930_0],
    [-0.816_5, -0.408_2, 0.408_2],
    [0.132_5, -0.662_3, -0.737_5],
];

/// Count how many triangles a ray from `from` crosses, or `None` if any hit was
/// too close to an edge or vertex to count reliably.
fn count_crossings(
    triangles: &[[Point; 3]],
    from: Point,
    along: Direction,
    tol: Tolerances,
) -> Option<usize> {
    let mut crossings = 0;
    for t in triangles {
        match ray_hits_triangle(from, along, *t, tol) {
            Hit::Crosses => crossings += 1,
            Hit::Misses => {}
            Hit::Ambiguous => return None,
        }
    }
    Some(crossings)
}

/// What a ray did to a triangle.
enum Hit {
    /// Passed through its interior, ahead of the start.
    Crosses,
    /// Did not meet it.
    Misses,
    /// Met an edge, a vertex, or the plane edge-on, where counting it once is
    /// as defensible as counting it twice or not at all.
    Ambiguous,
}

/// Möller–Trumbore, with the degenerate cases separated out rather than
/// rounded away.
fn ray_hits_triangle(from: Point, along: Direction, t: [Point; 3], tol: Tolerances) -> Hit {
    let direction = along.vector();
    let (e1, e2) = (t[1] - t[0], t[2] - t[0]);
    let h = direction.cross(e2);
    let determinant = e1.dot(h);

    // Scale the comparison by the triangle: a determinant is a volume, so a
    // fixed threshold rejects small triangles and accepts edge-on hits on
    // large ones.
    let scale = e1.magnitude() * e2.magnitude();
    let flat = tol.confusion() * scale;
    if determinant.abs() <= flat {
        // Edge-on. It cannot cross cleanly, but it may lie in the plane and
        // touch, which no crossing count expresses.
        let normal = e1.cross(e2);
        let reach = tol.confusion() * scale;
        return if normal.dot(from - t[0]).abs() <= reach {
            Hit::Ambiguous
        } else {
            Hit::Misses
        };
    }

    let inverse = 1.0 / determinant;
    let s = from - t[0];
    let u = inverse * s.dot(h);
    let q = s.cross(e1);
    let v = inverse * direction.dot(q);
    let w = 1.0 - u - v;

    // Barycentric coordinates near zero mean the ray passed along an edge, and
    // near one that it went through a vertex.
    let edge = tol.confusion();
    if [u, v, w].iter().any(|c| c.abs() <= edge) {
        // Only ambiguous if the ray would otherwise have hit: a grazing miss
        // well outside the triangle is a miss.
        return if u >= -edge && v >= -edge && w >= -edge {
            Hit::Ambiguous
        } else {
            Hit::Misses
        };
    }
    if u < 0.0 || v < 0.0 || w < 0.0 {
        return Hit::Misses;
    }

    let distance = inverse * e2.dot(q);
    if distance <= tol.confusion() {
        // Behind the start, or right at it — and right at it was already ruled
        // out by the on-boundary test before any ray was cast.
        return Hit::Misses;
    }
    Hit::Crosses
}

/// The distance from a point to a triangle.
fn distance_to_triangle(p: Point, t: [Point; 3]) -> f64 {
    // Clamp the projection onto the triangle's plane into the triangle, by
    // checking the three edge regions and the interior. Solving the 2×2 normal
    // equations directly and clamping is shorter than a region case analysis
    // and gives the same closest point.
    let (e1, e2) = (t[1] - t[0], t[2] - t[0]);
    let d = t[0] - p;
    let (a, b, c) = (e1.dot(e1), e1.dot(e2), e2.dot(e2));
    let (dd, e) = (e1.dot(d), e2.dot(d));
    let determinant = b.mul_add(-b, a * c);

    if determinant.abs() <= f64::MIN_POSITIVE {
        // A degenerate triangle is a segment; its edges still answer.
        return edge_distance(p, t);
    }
    let mut s = b.mul_add(e, -(c * dd)) / determinant;
    let mut u = b.mul_add(dd, -(a * e)) / determinant;

    if s >= 0.0 && u >= 0.0 && s + u <= 1.0 {
        let closest = t[0] + e1 * s + e2 * u;
        return p.distance(closest);
    }
    // Outside: the closest point is on an edge.
    s = s.clamp(0.0, 1.0);
    u = u.clamp(0.0, 1.0);
    let _ = (s, u);
    edge_distance(p, t)
}

/// The distance from a point to the nearest of a triangle's three edges.
fn edge_distance(p: Point, t: [Point; 3]) -> f64 {
    let mut best = f64::INFINITY;
    for i in 0..3 {
        best = best.min(segment_distance(p, t[i], t[(i + 1) % 3]));
    }
    best
}

/// The distance from a point to a segment.
fn segment_distance(p: Point, a: Point, b: Point) -> f64 {
    let d = b - a;
    let squared = d.dot(d);
    if squared <= f64::MIN_POSITIVE {
        return p.distance(a);
    }
    let t = ((p - a).dot(d) / squared).clamp(0.0, 1.0);
    p.distance(a + d * t)
}

/// The distance in parameter space from a point to the nearest ring.
fn distance_to_rings(rings: &[Vec<Point2>], p: Point2) -> f64 {
    let mut best = f64::INFINITY;
    for ring in rings {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            best = best.min(segment_distance_2d(p, a, b));
        }
    }
    best
}

/// The distance from a 2D point to a 2D segment.
fn segment_distance_2d(p: Point2, a: Point2, b: Point2) -> f64 {
    let d = b - a;
    let squared = d.dot(d);
    if squared <= f64::MIN_POSITIVE {
        return p.distance(a);
    }
    let t = ((p - a).dot(d) / squared).clamp(0.0, 1.0);
    p.distance(a + d * t)
}

/// How wide, in parameter units, a distance of `reach` in space is at `(u, v)`.
///
/// The surface's tangents give the conversion. Where a tangent vanishes — a
/// sphere's pole, a cone's apex — no parameter distance corresponds to a
/// spatial one, and the band opens to cover the whole neighbourhood rather than
/// closing to nothing.
fn parametric_band(
    surface: &og_geom::SurfaceGeometry,
    at: (f64, f64),
    reach: f64,
    tol: Tolerances,
) -> f64 {
    use og_geom::Surface;
    let Ok((du, dv)) = surface.d1_at(at.0, at.1, tol) else {
        return reach;
    };
    let scale = du.magnitude().min(dv.magnitude());
    if scale <= tol.confusion() {
        return f64::INFINITY;
    }
    reach / scale
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::make_box;
    use og_math::Frame;
    use og_topo::{ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-3,
            angular: 0.05,
            ..Deflection::default()
        }
    }

    #[test]
    fn a_point_inside_a_box_is_inside_it() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

        for p in [
            Point::new(1.0, 1.0, 1.0),
            Point::new(0.1, 0.1, 0.1),
            Point::new(1.9, 1.9, 1.9),
        ] {
            assert_eq!(
                classify_in_solid(&model, &built.shape, p, fine(), T).unwrap(),
                Containment::In,
                "{p:?} should be inside"
            );
        }
    }

    #[test]
    fn a_point_outside_a_box_is_outside_it() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

        for p in [
            Point::new(3.0, 1.0, 1.0),
            Point::new(-1.0, 1.0, 1.0),
            Point::new(1.0, 1.0, -0.5),
            Point::new(-5.0, -5.0, -5.0),
        ] {
            assert_eq!(
                classify_in_solid(&model, &built.shape, p, fine(), T).unwrap(),
                Containment::Out,
                "{p:?} should be outside"
            );
        }
    }

    #[test]
    fn a_point_on_a_boxs_face_is_on_it_rather_than_forced_to_a_side() {
        // The case a two-valued classifier has to guess at, and the case a
        // boolean spends all its time in.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

        for p in [
            Point::new(1.0, 1.0, 0.0), // face centre
            Point::new(0.0, 1.0, 1.0), // another face
            Point::new(2.0, 2.0, 1.0), // an edge
            Point::ORIGIN,             // a vertex
            Point::new(2.0, 2.0, 2.0), // the far vertex
        ] {
            assert_eq!(
                classify_in_solid(&model, &built.shape, p, fine(), T).unwrap(),
                Containment::On,
                "{p:?} should be on the boundary"
            );
        }
    }

    #[test]
    fn a_ray_along_a_grid_of_edges_does_not_defeat_the_classifier() {
        // A box tessellated into two triangles per face has a diagonal across
        // every face and an edge along every side. An axis-aligned ray from the
        // centre runs straight into a face centre; a diagonal one can run along
        // a triangle edge. The retry is what makes either survivable.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

        // The centre of a cube: every axis-aligned ray hits a face centre, and
        // the main diagonal goes through a vertex.
        assert_eq!(
            classify_in_solid(&model, &built.shape, Point::new(1.0, 1.0, 1.0), fine(), T).unwrap(),
            Containment::In
        );
    }

    #[test]
    fn an_open_shell_has_no_inside() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let face = explore_unique(&model, &built.shape, ShapeType::Face).unwrap()[0].clone();
        assert!(classify_in_solid(&model, &face, Point::ORIGIN, fine(), T).is_err());
    }

    #[test]
    fn a_point_on_a_face_is_inside_its_trimming_or_not() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();
        // The face at z = 0, spanning x in [0, 2] and y in [0, 3]. Found by
        // the role make_box gave it, which is what provenance is for.
        let bottom = explore_unique(&model, &built.shape, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                model.provenance_of(f).and_then(og_core::Provenance::role)
                    == Some(crate::primitive::roles::FACE_MIN_Z)
            })
            .expect("the box has a face at z = 0");

        assert_eq!(
            classify_on_face(&model, &bottom, Point::new(1.0, 1.5, 0.0), fine(), T).unwrap(),
            Containment::In
        );
        assert_eq!(
            classify_on_face(&model, &bottom, Point::new(5.0, 1.5, 0.0), fine(), T).unwrap(),
            Containment::Out,
            "on the surface's plane but outside the trimming"
        );
        assert_eq!(
            classify_on_face(&model, &bottom, Point::new(1.0, 1.5, 1.0), fine(), T).unwrap(),
            Containment::Out,
            "off the surface entirely"
        );
        assert_eq!(
            classify_on_face(&model, &bottom, Point::new(0.0, 1.5, 0.0), fine(), T).unwrap(),
            Containment::On,
            "on the trimming boundary"
        );
    }

    #[test]
    fn the_answers_invert_the_way_a_complement_does() {
        assert_eq!(Containment::In.inverted(), Containment::Out);
        assert_eq!(Containment::Out.inverted(), Containment::In);
        // A boundary belongs to both sides, so complementing leaves it alone.
        assert_eq!(Containment::On.inverted(), Containment::On);

        assert!(Containment::In.is_inside_or_on());
        assert!(Containment::On.is_inside_or_on());
        assert!(!Containment::Out.is_inside_or_on());
    }

    #[test]
    fn a_translated_box_classifies_the_same_way_translated_points() {
        let offset = Vector::new(10.0, -20.0, 30.0);
        let mut model = Model::new();
        let frame = Frame::new(Point::ORIGIN + offset, Direction::Z, Direction::X, T).unwrap();
        let built = make_box(&mut model, frame, (2.0, 2.0, 2.0), T).unwrap();

        assert_eq!(
            classify_in_solid(
                &model,
                &built.shape,
                Point::new(1.0, 1.0, 1.0) + offset,
                fine(),
                T
            )
            .unwrap(),
            Containment::In
        );
        assert_eq!(
            classify_in_solid(&model, &built.shape, Point::new(1.0, 1.0, 1.0), fine(), T).unwrap(),
            Containment::Out,
            "the untranslated point is nowhere near the translated box"
        );
    }

    #[test]
    fn the_exact_classifier_agrees_with_the_tessellated_one_on_a_box() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

        for (p, want) in [
            (Point::new(1.0, 1.0, 1.0), Containment::In),
            (Point::new(0.1, 0.1, 0.1), Containment::In),
            (Point::new(3.0, 1.0, 1.0), Containment::Out),
            (Point::new(1.0, 1.0, -0.5), Containment::Out),
            (Point::new(-50.0, -50.0, -50.0), Containment::Out),
            (Point::new(1.0, 1.0, 0.0), Containment::On),
            (Point::new(2.0, 2.0, 1.0), Containment::On),
            (Point::ORIGIN, Containment::On),
        ] {
            assert_eq!(
                classify_in_solid_exact(&model, &built.shape, p, T).unwrap(),
                want,
                "{p:?}"
            );
        }
    }

    #[test]
    fn the_exact_classifier_resolves_what_the_deflection_band_cannot() {
        // The reason this function exists. A point a micron off a sphere's
        // wall is far inside any practical deflection band — the tessellated
        // classifier must say On, because against a mesh it genuinely cannot
        // tell. Against the true sphere the side is knowable, and known.
        let mut model = Model::new();
        let built = crate::make_sphere(&mut model, Frame::WORLD, 2.0, T).unwrap();

        let barely_in = Point::new(0.0, 0.0, 2.0 - 1e-5);
        let barely_out = Point::new(0.0, 0.0, 2.0 + 1e-5);

        assert_eq!(
            classify_in_solid(&model, &built.shape, barely_in, fine(), T).unwrap(),
            Containment::On,
            "the mesh cannot tell a micron from the wall"
        );
        assert_eq!(
            classify_in_solid_exact(&model, &built.shape, barely_in, T).unwrap(),
            Containment::In
        );
        assert_eq!(
            classify_in_solid_exact(&model, &built.shape, barely_out, T).unwrap(),
            Containment::Out
        );
        // Exactly on the wall: On, decided by projection, not by a ray.
        assert_eq!(
            classify_in_solid_exact(&model, &built.shape, Point::new(0.0, 0.0, 2.0), T).unwrap(),
            Containment::On
        );
    }

    #[test]
    fn the_exact_classifier_handles_a_cylinder_wall_and_caps() {
        let mut model = Model::new();
        let built = crate::make_cylinder(&mut model, Frame::WORLD, 1.5, 4.0, T).unwrap();

        for (p, want) in [
            (Point::new(0.0, 0.0, 2.0), Containment::In),
            (Point::new(1.5 - 1e-5, 0.0, 2.0), Containment::In),
            (Point::new(1.5 + 1e-5, 0.0, 2.0), Containment::Out),
            (Point::new(0.3, 0.4, 4.0 - 1e-5), Containment::In),
            (Point::new(0.3, 0.4, 4.0 + 1e-5), Containment::Out),
            (Point::new(1.5, 0.0, 2.0), Containment::On),
            (Point::new(0.3, 0.4, 0.0), Containment::On),
        ] {
            assert_eq!(
                classify_in_solid_exact(&model, &built.shape, p, T).unwrap(),
                want,
                "{p:?}"
            );
        }
    }

    #[test]
    fn the_exact_classifier_walks_the_general_path_through_a_torus() {
        // No analytic ray/torus case exists, so every crossing here came from
        // the seeded Newton path — and a torus also puts the hole in the
        // middle, where a ray to the outside crosses the tube wall twice.
        let mut model = Model::new();
        let built = crate::make_torus(&mut model, Frame::WORLD, 3.0, 1.0, T).unwrap();

        for (p, want) in [
            (Point::new(3.0, 0.0, 0.0), Containment::In),
            (Point::new(3.0, 0.0, 0.9), Containment::In),
            (Point::ORIGIN, Containment::Out),
            (Point::new(3.0, 0.0, 1.5), Containment::Out),
            (Point::new(5.0, 5.0, 0.0), Containment::Out),
            (Point::new(3.0, 0.0, 1.0), Containment::On),
        ] {
            assert_eq!(
                classify_in_solid_exact(&model, &built.shape, p, T).unwrap(),
                want,
                "{p:?}"
            );
        }
    }

    #[test]
    fn the_exact_classifier_respects_a_placed_solid() {
        let offset = Vector::new(10.0, -20.0, 30.0);
        let mut model = Model::new();
        let frame = Frame::new(Point::ORIGIN + offset, Direction::Z, Direction::X, T).unwrap();
        let built = make_box(&mut model, frame, (2.0, 2.0, 2.0), T).unwrap();

        assert_eq!(
            classify_in_solid_exact(&model, &built.shape, Point::new(1.0, 1.0, 1.0) + offset, T)
                .unwrap(),
            Containment::In
        );
        assert_eq!(
            classify_in_solid_exact(&model, &built.shape, Point::new(1.0, 1.0, 1.0), T).unwrap(),
            Containment::Out
        );
    }

    #[test]
    fn the_exact_classifier_refuses_an_open_boundary() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let face = explore_unique(&model, &built.shape, ShapeType::Face).unwrap()[0].clone();
        assert!(classify_in_solid_exact(&model, &face, Point::ORIGIN, T).is_err());
    }

    #[test]
    fn distance_to_a_triangle_is_measured_from_the_nearest_part_of_it() {
        let t = [
            Point::ORIGIN,
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        // Above the interior: the plane distance.
        assert!((distance_to_triangle(Point::new(0.25, 0.25, 2.0), t) - 2.0).abs() < 1e-12);
        // Beyond a vertex: the distance to that vertex.
        assert!((distance_to_triangle(Point::new(-3.0, 0.0, 0.0), t) - 3.0).abs() < 1e-12);
        // On it: nothing.
        assert!(distance_to_triangle(Point::new(0.25, 0.25, 0.0), t) < 1e-12);
    }

    #[test]
    fn an_unusable_deflection_is_refused() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let bad = Deflection {
            chord: f64::NAN,
            ..Deflection::default()
        };
        assert!(classify_in_solid(&model, &built.shape, Point::ORIGIN, bad, T).is_err());
    }
}
