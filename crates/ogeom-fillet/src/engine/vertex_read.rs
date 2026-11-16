//! Geometry and admissibility of an isolated convex corner.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Direction, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType};

/// The kind of corner a vertex of a solid is to a ball of some radius.
pub(crate) enum Corner {
    /// Fewer than three planes pass through the corner: a curved face
    /// meets there.
    Curved(Point),
    /// No one ball touches every plane through the corner, whose inward
    /// normals these are.
    Setback(Point, Vec<Vector>),
    /// One ball touches every plane through a convex corner.
    OneBall(OneBall),
}

/// A convex planar corner one ball touches every face of.
pub(crate) struct OneBall {
    /// The vertex where it stands.
    pub(crate) corner: Point,
    /// The ball's centre.
    pub(crate) far: Point,
    /// The inward normal of each plane through the corner.
    pub(crate) m: Vec<Vector>,
    /// The planes in their order round the vertex.
    pub(crate) ring: Vec<usize>,
    /// Edge `k`, from the vertex into the solid, where plane `ring[k]` and
    /// plane `ring[k + 1]` meet.
    pub(crate) d: Vec<Direction>,
    /// The points whose containment in the solid told the corner convex.
    pub(crate) probes: Vec<Point>,
}

impl OneBall {
    /// The box round the corner's block: the vertex, the foot of each edge
    /// on its plane through the ball's centre, the point where the ball
    /// touches each plane, and the centre.
    pub(crate) fn block_bounds(&self, radius: f64) -> ogeom_math::Aabb {
        let mut points = vec![self.corner, self.far];
        for (dk, mk) in self.d.iter().zip(&self.m) {
            points.push(self.foot(dk.vector()));
            points.push(self.far - *mk * radius);
        }
        ogeom_math::Aabb::of_points(&points)
    }

    /// Where the plane through the ball's centre square to `along` meets
    /// the edge leaving the vertex along it.
    pub(crate) fn foot(&self, along: Vector) -> Point {
        self.corner + along * (self.far - self.corner).dot(along)
    }
}

/// The corner `vertex` of `solid` is to a ball of `radius`, and every check
/// that the tool speaks it.
pub(crate) fn corner_at(
    model: &Model,
    solid: &Shape,
    vertex: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Corner> {
    if model.kind_of(vertex)? != ShapeType::Vertex {
        ogeom_bail!(Construction, "round_vertex rounds a vertex");
    }
    if !(radius.is_finite() && radius > tol.confusion()) {
        ogeom_bail!(Construction, "a blend radius must be a positive distance");
    }
    let Some(raw) = model
        .node(vertex)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
    else {
        ogeom_bail!(Construction, "the vertex holds no point");
    };
    // Where the vertex stands, not where its node was built: a prism's far
    // end is its near end moved, and read unplaced the far corner is the
    // near one, so the tool would round the wrong corner of the solid.
    let corner = vertex.transform(model.datums())?.apply(raw);

    // The corner's frame comes from the planes that pass through the
    // vertex's point, not from the vertex's own adjacency, which the very
    // sequence this tool serves destroys: after three fillets the tip is
    // consumed, but the three shrunk planes still contain the corner, and
    // still say exactly which corner it was. The vertex argument may
    // therefore come from an earlier state of the solid (the sharp box's
    // corner captured before the fillets) and anchors the history. Each
    // plane's material side is read off its face, which the fillets shrink
    // but never turn over.
    let mut m: Vec<Vector> = Vec::new();
    for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            continue;
        };
        if !matches!(surface, ogeom_geom::SurfaceGeometry::Plane(_)) {
            continue;
        }
        let (origin, outward) = ogeom_algo::face_normal(model, &face, tol)?;
        if (corner - origin).dot(outward).abs() > tol.confusion() * 100.0 {
            continue;
        }
        // One vote per plane: coplanar trims share it.
        if m.iter()
            .any(|n| n.cross(outward).magnitude() < tol.angular() * 10.0)
        {
            continue;
        }
        m.push(-outward);
    }
    let n = m.len();
    if n < 3 {
        // A curved face meets here: the one-ball corner on whatever
        // surfaces pass through the vertex.
        return Ok(Corner::Curved(corner));
    }
    // The ball's centre: the point a radius in from every plane. Three
    // planes that span always hold one; more only when they share a
    // tangent ball, which the least-squares fit's residual tells: a
    // square pyramid's apex does, a general N-edged vertex does not, and
    // that vertex takes the offset-region envelope instead.
    let centre_for = |m: &[Vector]| -> Option<(Point, f64)> {
        // m_k · (c − corner) = radius over all k: solved directly for
        // three planes, through the normal equations for more. Three are
        // solved directly for the last bit, not for speed: the boolean's
        // paving at the touch points is sensitive to an ulp of the centre,
        // and an oblique corner closes on the direct solve's value.
        let mut a = [[0.0_f64; 3]; 3];
        let mut b = [0.0_f64; 3];
        if m.len() == 3 {
            for (i, mk) in m.iter().enumerate() {
                a[i] = [mk.x, mk.y, mk.z];
                b[i] = radius;
            }
        } else {
            for mk in m {
                let v = [mk.x, mk.y, mk.z];
                for i in 0..3 {
                    b[i] += radius * v[i];
                    for j in 0..3 {
                        a[i][j] += v[i] * v[j];
                    }
                }
            }
        }
        let det3 = |a: &[[f64; 3]; 3]| -> f64 {
            a[0][0].mul_add(
                a[1][1].mul_add(a[2][2], -(a[1][2] * a[2][1])),
                -a[0][1].mul_add(
                    a[1][0].mul_add(a[2][2], -(a[1][2] * a[2][0])),
                    -(a[0][2] * a[1][0].mul_add(a[2][1], -(a[1][1] * a[2][0]))),
                ),
            )
        };
        let det = det3(&a);
        if !det.is_finite() || det.abs() <= 1e-12 {
            return None;
        }
        let mut x = [0.0_f64; 3];
        for (k, xk) in x.iter_mut().enumerate() {
            let mut ak = a;
            for (row, bi) in ak.iter_mut().zip(b) {
                row[k] = bi;
            }
            *xk = det3(&ak) / det;
        }
        let centre = corner + Vector::new(x[0], x[1], x[2]);
        let residual = m
            .iter()
            .map(|mk| (mk.dot(centre - corner) - radius).abs())
            .fold(0.0, f64::max);
        Some((centre, residual))
    };
    let Some((far, residual)) = centre_for(&m) else {
        ogeom_bail!(
            Construction,
            "the {n} planes through this vertex do not span a corner"
        );
    };
    if residual > tol.confusion() * 10.0 {
        // No one ball touches every face: the rounded corner is the
        // envelope of every ball a radius in from all of them, which is
        // more than one sphere. Its pieces are read off the region the
        // ball's centre may occupy.
        return Ok(Corner::Setback(corner, m));
    }
    // A concave vertex puts the ball's centre outside the material: a ball
    // there adds material instead of shedding it, and a tool built from a
    // cut cannot say so.
    let boundary = ogeom_algo::SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)?;
    let mut probes = vec![far];
    if boundary.holds(model, far, tol)? != ogeom_algo::Containment::In {
        ogeom_bail!(
            Construction,
            "no material a radius in from every face at this vertex; a \
             concave vertex gains a ball instead of shedding one, and this \
             tool cannot round it"
        );
    }
    // The corner's edges: where two of the planes meet along a ray that
    // lies inside all the others, directed from the vertex into the solid.
    // A simple convex corner has one edge per plane, and the edges chain the
    // planes into a ring round the vertex.
    let mut edges: Vec<(usize, usize, Direction)> = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let cross = m[i].cross(m[j]);
            let length = cross.magnitude();
            if length <= tol.angular() * 10.0 {
                continue;
            }
            for sign in [1.0, -1.0] {
                let d = cross * (sign / length);
                let inside = (0..n)
                    .filter(|&k| k != i && k != j)
                    .all(|k| m[k].dot(d) >= -tol.angular() * 10.0);
                if inside {
                    edges.push((i, j, Direction::new(d, tol)?));
                    break;
                }
            }
        }
    }
    if edges.len() != n {
        ogeom_bail!(
            Construction,
            "the {n} planes through this vertex meet along {} edges, not {n}; \
             the corner is not the simple convex one this tool speaks",
            edges.len()
        );
    }
    // The planes taken whole make any corner look convex, but a re-entrant
    // one (an L-bracket's step, where a concave edge meets two convex ones)
    // has material outside a plane beside one of its edges. At a convex
    // corner the solid lies within every plane near the vertex, however
    // much of it earlier blends have cut away, so a point just off one of
    // an edge's planes, and in from the other, is outside it.
    {
        let boundary = ogeom_algo::SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)?;
        let step = (radius * 1e-2).max(tol.confusion() * 1e3);
        for &(i, j, d) in &edges {
            for (k, other) in [(i, j), (j, i)] {
                let off = corner + (d.vector() - m[k] + m[other]) * step;
                probes.push(off);
                if boundary.holds(model, off, tol)? == ogeom_algo::Containment::In {
                    ogeom_bail!(
                        Construction,
                        "material stands outside a plane through this vertex beside \
                         one of its edges; the corner is re-entrant, not the simple \
                         convex one this tool speaks"
                    );
                }
            }
        }
    }
    let mut ring: Vec<usize> = vec![0];
    let mut last_edge: Option<usize> = None;
    loop {
        let here = *ring.last().unwrap_or(&0);
        let Some((e, other)) = edges.iter().enumerate().find_map(|(e, &(i, j, _))| {
            if Some(e) == last_edge {
                None
            } else if i == here {
                Some((e, j))
            } else if j == here {
                Some((e, i))
            } else {
                None
            }
        }) else {
            ogeom_bail!(
                Construction,
                "the planes through this vertex do not chain into a ring round it"
            );
        };
        last_edge = Some(e);
        if other == 0 {
            break;
        }
        if ring.contains(&other) || ring.len() >= n {
            ogeom_bail!(
                Construction,
                "the planes through this vertex do not chain into a ring round it"
            );
        }
        ring.push(other);
    }
    if ring.len() != n {
        ogeom_bail!(
            Construction,
            "the planes through this vertex do not chain into a ring round it"
        );
    }
    // Plane `ring[k]` and plane `ring[k+1]` meet along edge k.
    let d: Vec<Direction> = (0..n)
        .map(|k| {
            let (p, q) = (ring[k], ring[(k + 1) % n]);
            edges
                .iter()
                .find(|&&(i, j, _)| (i == p && j == q) || (i == q && j == p))
                .map(|&(_, _, dir)| dir)
        })
        .collect::<Option<_>>()
        .ok_or_else(|| {
            ogeom_core::ogeom_err!(Construction, "the corner's edges lost their ring")
        })?;
    Ok(Corner::OneBall(OneBall {
        corner,
        far,
        m,
        ring,
        d,
        probes,
    }))
}

/// A convex planar corner one ball touches every face of, read as
/// [`crate::round_vertex`] reads it: the vertex's point, the ball's centre, each
/// plane's inward normal, the planes in their order round the vertex, and
/// edge `k`'s direction into the solid, where plane `ring[k]` meets plane
/// `ring[k + 1]`.
pub(crate) struct CornerFrame {
    pub corner: Point,
    pub far: Point,
    pub m: Vec<Vector>,
    pub ring: Vec<usize>,
    pub d: Vec<Direction>,
}

/// The frame of `vertex` of `solid` as a corner one ball of `radius`
/// touches every plane of, refused by name where it is another kind.
pub(crate) fn corner_frame(
    model: &Model,
    solid: &Shape,
    vertex: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<CornerFrame> {
    match corner_at(model, solid, vertex, radius, tol)? {
        Corner::OneBall(one) => Ok(CornerFrame {
            corner: one.corner,
            far: one.far,
            m: one.m,
            ring: one.ring,
            d: one.d,
        }),
        Corner::Curved(_) => ogeom_bail!(
            Construction,
            "a setback corner where a curved face meets the vertex is still \
             owed; see docs/CAPABILITIES.md, fillet.setback-corners"
        ),
        Corner::Setback(..) => ogeom_bail!(
            Construction,
            "a setback corner whose planes share no tangent ball is still \
             owed; see docs/CAPABILITIES.md, fillet.setback-corners"
        ),
    }
}
