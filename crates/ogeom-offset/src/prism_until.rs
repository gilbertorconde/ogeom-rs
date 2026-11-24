//! A prism that stops on faces: each line of the sweep ends where it first
//! meets one of the limits, at a different length for each.

use ogeom_algo::sweep::roles::{SWEEP_BOTTOM, SWEEP_TOP};
use ogeom_algo::{Built, Containment, CurveHit, Hit, curve_hit};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, LineCurve, SurfaceGeometry};
use ogeom_math::{Aabb, Axis, Direction, Plane, Point, Transform, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Filter, Location, Model, NodeData, Shape, ShapeType, explore, explore_unique};

/// Extrude `profile` (a planar face) along `direction`, each line of the
/// sweep ending where it first meets one of `limits`, moved `offset` along
/// `direction` first: the prism "up to face" when a limit is a face, and
/// "up to next" when it is a part, whose faces then all count.
///
/// The rules it keeps:
///
/// - A limit is met as trimmed, never as its whole surface: each line
///   stops on the first limit face it crosses inside that face's trim, as
///   [`curve_hit`] finds it. A line that meets no limit is refused by name;
///   to run on to a face's whole surface, bound the prism by that
///   surface's half space ([`make_half_space`](ogeom_algo::make_half_space))
///   instead.
/// - Only what lies ahead of the profile counts: a limit, or the part of
///   one, behind the profile's plane stops nothing, so a limit wholly
///   behind leaves its lines with nothing to meet. The direction is the
///   caller's and is never turned round.
/// - Lines stop at different lengths: two limits over different parts of
///   the profile give a stepped end, each part ending on the nearer limit
///   above it, with a wall where the step is.
/// - A limit that bends back over itself, a cylinder's side across the
///   travel, stops each line where it first meets it: nothing beyond that
///   meeting survives, the far side included.
/// - A face the direction runs along (a plane holding the direction, a
///   cylinder about it) stops no line: its lines graze it and cross
///   nothing. A curved limit whose silhouette, where it turns edge-on to
///   the direction, falls over the profile is refused: the lines there
///   only touch it, and whether a touch stops a line is not decided by
///   guessing.
/// - A negative `offset` stops each line short of its limit by that much
///   along the direction, a positive one carries it past.
///
/// History: each profile edge generates its pieces of wall, the profile
/// survives as the near end, and each far-end face is derived from the
/// limit face it lies on. The near end carries
/// [`SWEEP_BOTTOM`](ogeom_algo::sweep::roles::SWEEP_BOTTOM) derived from
/// the profile, each far-end face
/// [`SWEEP_TOP`](ogeom_algo::sweep::roles::SWEEP_TOP) derived from its
/// limit face, which also generates it in the history.
///
/// Built as a prism long enough to pass every limit, less the region each
/// limit face casts ahead of itself along the direction: the face cut to
/// the column of the profile's lines ahead of its plane, and swept on past
/// the prism's end. The result is then measured line by
/// line against the limits before it is returned: a prism whose length
/// along any measured line differs from where [`curve_hit`] finds that
/// line's limit is refused rather than returned.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the profile is not a planar face, the direction lies in its plane, the
/// offset is not finite, the limits hold no face, every limit lies behind
/// the profile, some line meets no limit, or a curved limit face turns
/// edge-on to the direction over the profile;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the built
/// prism does not end where its lines meet the limits; as the prism and
/// the booleans otherwise.
pub fn make_prism_until(
    model: &mut Model,
    profile: &Shape,
    direction: Direction,
    limits: &[Shape],
    offset: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if model.kind_of(profile)? != ShapeType::Face {
        ogeom_bail!(Construction, "a prism up to faces extrudes a planar face");
    }
    let Some(plane) = ogeom_algo::find_plane(model, profile, tol)? else {
        ogeom_bail!(Construction, "a prism up to faces extrudes a planar face");
    };
    if !offset.is_finite() {
        ogeom_bail!(Construction, "an offset of {offset} is not a distance");
    }
    let d = direction.vector();
    let rise = d.dot(plane.normal().vector());
    if rise.abs() <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the direction lies in the profile's plane; a prism leaves its plane"
        );
    }
    let ahead = Ahead {
        origin: plane.origin(),
        normal: plane.normal().vector(),
        rise,
    };

    // The limits' faces, moved by the offset: each kept beside the face
    // the caller named, which the far end is derived from.
    let shift =
        (offset != 0.0).then(|| Location::of(model.add_datum(Transform::translation(d * offset))));
    let mut faces: Vec<(Shape, Shape)> = Vec::new();
    for limit in limits {
        for face in explore(model, limit, Filter::OfType(ShapeType::Face))? {
            let moved = shift
                .as_ref()
                .map_or_else(|| face.clone(), |by| face.moved(by));
            faces.push((face, moved));
        }
    }
    if faces.is_empty() {
        ogeom_bail!(Construction, "the limits hold no face to stop on");
    }

    // A prism reaching past every limit.
    let profile_box = ogeom_algo::shape_bounds(model, profile, tol)?;
    let mut far = f64::NEG_INFINITY;
    for (_, moved) in &faces {
        far = far.max(ahead.span(&ogeom_algo::shape_bounds(model, moved, tol)?).1);
    }
    if far <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "every limit lies behind the profile; a prism up to faces runs \
             along its direction only"
        );
    }
    let reach = far + 0.25 * (far + profile_box.diagonal()) + 1.0;
    let prism = ogeom_algo::make_prism(model, profile, d * reach, tol)?;
    let mut history = prism.history;
    let mut solid = prism.shape;

    // Less the region each limit face casts ahead of itself, each face
    // first cut to the column of lines the prism sweeps: what of it lies
    // outside that column, or behind the profile, stops none of them.
    let column = ogeom_algo::make_prism(model, profile, d * (2.0 * reach), tol)?.shape;
    for (_, moved) in &faces {
        for cast in casting(model, moved, &column, &ahead, direction, tol)? {
            let shadow = ogeom_algo::make_prism(model, &cast, d * (2.0 * reach), tol)?.shape;
            let cut = ogeom_bool::cut(model, &solid, &shadow, tol)?;
            history = history.then(&cut.history);
            solid = cut.shape;
            if explore(model, &solid, Filter::OfType(ShapeType::Face))?.is_empty() {
                ogeom_bail!(
                    Construction,
                    "the limits meet the profile itself; nothing is left to extrude"
                );
            }
        }
    }

    // A line that met nothing still reaches the long prism's end.
    let end = reach - (tol.confusion() * 100.0).max(reach * 1e-9);
    for vertex in explore_unique(model, &solid, ShapeType::Vertex)? {
        let at = point_of(model, &vertex)?;
        if ahead.along(at) >= end {
            ogeom_bail!(
                Construction,
                "the line through {:?} meets no limit; a prism up to faces stops \
                 every line on one",
                at - d * ahead.along(at)
            );
        }
    }

    measured(
        model, profile, &plane, &solid, &faces, direction, reach, tol,
    )?;
    rolled(model, profile, &ahead, &solid, &faces, &mut history, tol)?;
    history.generate(profile, solid.clone());
    Ok(Built::new(solid, history))
}

/// How far along the direction a point lies ahead of the profile's plane.
struct Ahead {
    origin: Point,
    normal: Vector,
    rise: f64,
}

impl Ahead {
    /// The travel along the direction from the profile's plane to `p`.
    fn along(&self, p: Point) -> f64 {
        (p - self.origin).dot(self.normal) / self.rise
    }

    /// The least and greatest travel over a box's corners.
    fn span(&self, bounds: &Aabb) -> (f64, f64) {
        bounds
            .corners()
            .into_iter()
            .map(|c| self.along(c))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), t| {
                (lo.min(t), hi.max(t))
            })
    }
}

/// A vertex's point, placed.
fn point_of(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_bail!(Construction, "a vertex node holds no vertex data");
    };
    Ok(vertex.transform(model.datums())?.apply(data.point))
}

/// A face's surface, placed.
fn placed_surface(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
    use ogeom_geom::Transformable as _;
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Construction, "a face node holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "a face refers to a surface not in this model");
    };
    surface.transformed(&face.transform(model.datums())?, tol)
}

/// What of `face` casts a region onto the prism's lines, as faces each
/// of which sweeps along `direction` into a solid: nothing for a face the
/// direction runs along, which stops no line, nor for one wholly behind
/// the profile; a plane wholly ahead as it is; any other face cut to
/// `column`, the profile's lines from its plane on past every limit.
///
/// A piece that still turns edge-on to the direction inside the column
/// (a cylinder whose silhouette falls over the profile) is refused: its
/// sweep would fold over itself.
fn casting(
    model: &mut Model,
    face: &Shape,
    column: &Shape,
    ahead: &Ahead,
    direction: Direction,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let d = direction.vector();
    let (lo, hi) = ahead.span(&ogeom_algo::shape_bounds(model, face, tol)?);
    if hi <= tol.confusion() {
        return Ok(Vec::new());
    }
    let planar = match placed_surface(model, face, tol)? {
        SurfaceGeometry::Plane(p) => {
            if p.plane().normal().vector().dot(d).abs() <= tol.angular() {
                return Ok(Vec::new());
            }
            true
        }
        SurfaceGeometry::Cylinder(c) => {
            let axis = c.cylinder().axis().direction.vector();
            if axis.cross(d).magnitude() <= tol.angular() {
                return Ok(Vec::new());
            }
            false
        }
        _ => false,
    };
    if planar && lo >= tol.confusion() {
        return Ok(vec![face.clone()]);
    }
    let kept = ogeom_bool::common(model, face, column, tol)?.shape;
    let pieces = explore(model, &kept, Filter::OfType(ShapeType::Face))?;
    if !planar {
        for piece in &pieces {
            let (lo, hi) = normal_spread(model, piece, d, tol)?;
            if lo <= 0.0 && hi >= 0.0 {
                ogeom_bail!(
                    Construction,
                    "a limit face turns edge-on to the direction over the profile; \
                     its silhouette there leaves lines that graze it, and a prism \
                     up to faces is not built across one"
                );
            }
        }
    }
    Ok(pieces)
}

/// The least and greatest of the unit normal's component along `d` over a
/// face's mesh.
fn normal_spread(
    model: &Model,
    face: &Shape,
    d: Vector,
    tol: Tolerances,
) -> OgeomResult<(f64, f64)> {
    let size = ogeom_algo::shape_bounds(model, face, tol)?.diagonal();
    let mesh = ogeom_mesh::triangulate(
        model,
        face,
        Deflection::with_chord((size * 1e-3).max(tol.confusion()))?,
        tol,
    )?;
    let mut spread = (f64::INFINITY, f64::NEG_INFINITY);
    for n in &mesh.normals {
        let m = n.magnitude();
        if m > 0.0 {
            let c = n.dot(d) / m;
            spread = (spread.0.min(c), spread.1.max(c));
        }
    }
    Ok(spread)
}

/// Refuse a prism whose length along a lattice of the profile's lines is
/// not where [`curve_hit`] finds each line's limit.
#[allow(clippy::too_many_arguments)]
fn measured(
    model: &Model,
    profile: &Shape,
    plane: &Plane,
    solid: &Shape,
    faces: &[(Shape, Shape)],
    direction: Direction,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<()> {
    const N: u32 = 9;
    let frame = plane.frame();
    let bounds = ogeom_algo::shape_bounds(model, profile, tol)?;
    let local: Vec<Point> = bounds
        .corners()
        .iter()
        .map(|c| frame.to_local(*c))
        .collect();
    let (mut ulo, mut uhi, mut vlo, mut vhi) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for p in &local {
        (ulo, uhi, vlo, vhi) = (ulo.min(p.x), uhi.max(p.x), vlo.min(p.y), vhi.max(p.y));
    }
    let deflection = Deflection::with_chord((bounds.diagonal() * 1e-4).max(tol.confusion()))?;
    let bound = (tol.confusion() * 100.0).max(reach * 1e-9);
    let start = (tol.confusion() * 10.0).max(reach * 1e-9);
    // Off the lattice a profile is usually drawn on, so lines rarely fall
    // on the edge of a step.
    let shift = 0.137_f64;
    for i in 0..N {
        for j in 0..N {
            let u = ulo + (uhi - ulo) * (f64::from(i) + 0.5 + shift) / f64::from(N);
            let v = vlo + (vhi - vlo) * (f64::from(j) + 0.5 - shift) / f64::from(N);
            let p = frame.to_world(Point::new(u, v, 0.0));
            if ogeom_algo::classify_on_face(model, profile, p, deflection, tol)? != Containment::In
            {
                continue;
            }
            let line: Curve = LineCurve::over(
                Axis {
                    location: p,
                    direction,
                },
                0.0,
                reach,
            )?
            .into();
            let mut met: Option<CurveHit> = None;
            for (_, moved) in faces {
                if let Some(hit) = curve_hit(model, moved, &line, (start, reach), Hit::First, tol)?
                    && met
                        .as_ref()
                        .is_none_or(|held| hit.parameter < held.parameter)
                {
                    met = Some(hit);
                }
            }
            let Some(met) = met else {
                ogeom_bail!(
                    Construction,
                    "the line through {p:?} meets no limit; a prism up to faces \
                     stops every line on one"
                );
            };
            let ends = curve_hit(model, solid, &line, (start, reach), Hit::First, tol)?;
            match ends {
                Some(end) if (end.parameter - met.parameter).abs() <= bound => {}
                _ => ogeom_bail!(
                    NotDone,
                    "the prism along the line through {p:?} ends at {:?}, where the \
                     line meets its limit at {}; a prism up to faces is not built here",
                    ends.map(|e| e.parameter),
                    met.parameter
                ),
            }
        }
    }
    Ok(())
}

/// Give the near end and the far-end faces their roles and history.
fn rolled(
    model: &mut Model,
    profile: &Shape,
    ahead: &Ahead,
    solid: &Shape,
    faces: &[(Shape, Shape)],
    history: &mut ogeom_algo::History,
    tol: Tolerances,
) -> OgeomResult<()> {
    let size = ogeom_algo::shape_bounds(model, solid, tol)?.diagonal();
    let deflection = Deflection::with_chord((size * 1e-4).max(tol.confusion()))?;
    for face in explore_unique(model, solid, ShapeType::Face)? {
        let vertices = explore_unique(model, &face, ShapeType::Vertex)?;
        let mut near = !vertices.is_empty();
        for vertex in &vertices {
            near &= ahead.along(point_of(model, vertex)?).abs() <= tol.confusion() * 10.0;
        }
        if near {
            model.set_derived(&face, std::slice::from_ref(profile), SWEEP_BOTTOM)?;
            continue;
        }
        let Some(inner) = inner_point(model, &face, tol)? else {
            continue;
        };
        for (named, moved) in faces {
            if ogeom_algo::classify_on_face(model, moved, inner, deflection, tol)?
                != Containment::Out
            {
                model.set_derived(&face, std::slice::from_ref(named), SWEEP_TOP)?;
                history.generate(named, face.clone());
                break;
            }
        }
    }
    Ok(())
}

/// A point of a face away from its boundary: the middle of its first mesh
/// triangle, taken onto its surface.
fn inner_point(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<Point>> {
    let size = ogeom_algo::shape_bounds(model, face, tol)?.diagonal();
    let mesh = ogeom_mesh::triangulate(
        model,
        face,
        Deflection::with_chord((size * 1e-2).max(tol.confusion()))?,
        tol,
    )?;
    let Some(t) = mesh.triangles.first() else {
        return Ok(None);
    };
    let [a, b, c] = t.map(|k| mesh.positions[k as usize]);
    let middle = Point::from_vector((a.to_vector() + b.to_vector() + c.to_vector()) / 3.0);
    let surface = placed_surface(model, face, tol)?;
    Ok(Some(
        ogeom_algo::project_on_surface(&surface, middle, 16, tol)?.point,
    ))
}
