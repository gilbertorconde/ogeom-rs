//! Tangent propagation: the edges a blend picked on one edge runs on to.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, NodeData, Shape, ShapeType, explore_unique};

/// The widest turn, in radians, at which two edges meeting at a vertex
/// still run on as one tangent chain by default: one degree (ADR 0001,
/// amendment "tangent chains stop at a kink").
pub const TANGENT_CHAIN_ANGLE: f64 = core::f64::consts::PI / 180.0;

/// Which vertices a tangent chain runs on through. [`tangent_chain`] is
/// [`tangent_chain_with`] under [`ChainRules::new`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChainRules {
    /// The widest turn, in radians, still taken as tangent: between the
    /// way one edge leaves a vertex and the way the next arrives, and
    /// between two faces' normals for [`through_fading`](Self::through_fading).
    pub angle: f64,
    /// Also run on where both edges' creases fade at the shared vertex:
    /// each edge's two faces have normals agreeing within `angle` there,
    /// whatever the angle between the edges. Edges with no crease anywhere
    /// (a seam, a smooth edge) are never taken this way.
    pub through_fading: bool,
    /// Also run on where two edges meet in a cusp: both leave the vertex
    /// the same way, within `angle`.
    pub through_cusps: bool,
}

impl ChainRules {
    /// Edges running on into each other within `angle` only: no fading
    /// creases, no cusps.
    #[must_use]
    pub const fn new(angle: f64) -> Self {
        Self {
            angle,
            through_fading: false,
            through_cusps: false,
        }
    }
}

impl Default for ChainRules {
    fn default() -> Self {
        Self::new(TANGENT_CHAIN_ANGLE)
    }
}

/// `edges` and every edge of `solid` reached from them through vertices
/// where two edges run on into each other, turning by at most `angle`
/// radians there ([`TANGENT_CHAIN_ANGLE`] by default).
///
/// This is the tangent propagation of a fillet or chamfer request: a
/// picked side of a rounded outline brings the arcs and the sides past
/// them, and the chain stops where the outline turns a corner. A kink
/// wider than `angle` stops it, however slight it looks. The edges come
/// back with `edges` first, in their order, then each edge reached, in the
/// order the walk reaches it, so the same request gives the same list.
///
/// Two edges run on into each other at a vertex they share when the way
/// one leaves the vertex is the way the other arrives, within `angle`.
/// Edges meeting back to back (a cusp) do not, nor do two edges meeting at
/// an angle where their creases fade; [`tangent_chain_with`] takes those
/// on request. The angle chooses edges only: the blend engine runs bands
/// on into each other with no corner only where the edges are tangent, and
/// closes any other vertex of the chain as a corner.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `angle` is not in `[0, π/2)`, or an edge of `edges` is not one of
/// `solid`'s.
pub fn tangent_chain(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    angle: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    tangent_chain_with(model, solid, edges, ChainRules::new(angle), tol)
}

/// [`tangent_chain`] under `rules`: also through vertices where both
/// edges' creases fade ([`ChainRules::through_fading`]), or where two
/// edges meet in a cusp ([`ChainRules::through_cusps`]).
///
/// A crease fades at a vertex where its two faces turn tangent to each
/// other: two equal drums crossing at right angles meet along two
/// ellipses, and the four arcs meeting at each point where the drums
/// touch all fade there. Picking one arc takes its ellipse alone by
/// default, and the whole seam through fading creases (ADR 0001,
/// amendment "tangent chains through fading creases").
///
/// # Errors
///
/// As [`tangent_chain`], for `rules.angle`.
pub fn tangent_chain_with(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    rules: ChainRules,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let angle = rules.angle;
    if !(0.0..core::f64::consts::FRAC_PI_2).contains(&angle) {
        ogeom_bail!(
            Construction,
            "a tangent chain's angle of {angle} rad is not in [0, π/2)"
        );
    }
    let all = explore_unique(model, solid, ShapeType::Edge)?;
    let mut ends: Vec<Vec<End>> = all
        .iter()
        .map(|e| edge_ends(model, e, tol))
        .collect::<OgeomResult<_>>()?;
    let cos = angle.cos();
    if rules.through_fading {
        for (edge, ends) in all.iter().zip(&mut ends) {
            mark_fading(model, solid, edge, ends, cos, tol)?;
        }
    }
    let mut taken = vec![false; all.len()];
    let mut chain: Vec<Shape> = Vec::with_capacity(edges.len());
    let mut queue: Vec<usize> = Vec::new();
    for edge in edges {
        let Some(i) = all.iter().position(|e| e.is_same(edge)) else {
            ogeom_bail!(
                Construction,
                "an edge of the tangent chain is not one of the solid's edges"
            );
        };
        if !taken[i] {
            taken[i] = true;
            chain.push(edge.clone());
            queue.push(i);
        }
    }
    let mut next = 0;
    while next < queue.len() {
        let i = queue[next];
        next += 1;
        for (j, theirs) in ends.iter().enumerate() {
            if taken[j] {
                continue;
            }
            let meets = ends[i].iter().any(|a| {
                theirs.iter().any(|b| {
                    let turn = a.leaving.dot(b.leaving);
                    a.vertex == b.vertex
                        && a.at.distance(b.at) <= a.reach.max(b.reach)
                        && (-turn >= cos
                            || (rules.through_cusps && turn >= cos)
                            || (a.fades && b.fades))
                })
            });
            if meets {
                taken[j] = true;
                chain.push(all[j].clone());
                queue.push(j);
            }
        }
    }
    Ok(chain)
}

/// Marks the ends of `edge` where its crease fades: its two faces' outward
/// normals agree within the angle of cosine `cos` at the end, and disagree
/// somewhere inside the edge. An edge with other than two faces, or with
/// one face on both sides (a seam), has no crease to fade.
fn mark_fading(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    ends: &mut [End],
    cos: f64,
    tol: Tolerances,
) -> OgeomResult<()> {
    if ends.is_empty() {
        return Ok(());
    }
    let faces = ogeom_topo::ancestors_of(model, solid, edge, ShapeType::Face)?;
    let [f, g] = faces.as_slice() else {
        return Ok(());
    };
    if f.is_same(g) {
        return Ok(());
    }
    let agree = |p: Point| -> OgeomResult<Option<bool>> {
        let (Some(n), Some(m)) = (outward(model, f, p, tol)?, outward(model, g, p, tol)?) else {
            return Ok(None);
        };
        Ok(Some(n.dot(m) >= cos))
    };
    let Ok((curve, (s, t))) = crate::support::edge_curve(model, edge, tol) else {
        return Ok(());
    };
    let mut creased = false;
    for k in 1..8 {
        let p = curve.point_at(s + (t - s) * f64::from(k) / 8.0, tol)?;
        if agree(p)? == Some(false) {
            creased = true;
            break;
        }
    }
    if !creased {
        return Ok(());
    }
    for end in ends {
        end.fades = agree(end.at)? == Some(true);
    }
    Ok(())
}

/// `face`'s outward normal at the point of it nearest `p`.
fn outward(model: &Model, face: &Shape, p: Point, tol: Tolerances) -> OgeomResult<Option<Vector>> {
    let n = crate::marched::face_normal_near(model, face, p, tol)?;
    Ok(n.map(|n| {
        if face.orientation() == ogeom_topo::Orientation::Reversed {
            -n
        } else {
            n
        }
    }))
}

/// An end of an edge: the vertex there, its point, the unit direction the
/// edge leaves the vertex in, and how far apart two ends at one vertex may
/// stand (the vertex's tolerance), and whether the edge's crease fades
/// there (set only under [`ChainRules::through_fading`]).
struct End {
    vertex: u32,
    at: Point,
    leaving: Vector,
    reach: f64,
    fades: bool,
}

/// `edge`'s ends, placed in the world; none for an edge with no curve, a
/// degenerate one, or one that leaves a vertex with no direction.
fn edge_ends(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Vec<End>> {
    let Some((first, last)) = ogeom_algo::edge_vertices(model, edge)? else {
        return Ok(Vec::new());
    };
    let Ok((curve, (s, t))) = crate::support::edge_curve(model, edge, tol) else {
        return Ok(Vec::new());
    };
    let vertex = |v: &Shape| -> Option<(u32, Point, f64)> {
        match model.node(v)?.data() {
            NodeData::Vertex(data) => Some((
                v.node().index(),
                v.transform(model.datums()).ok()?.apply(data.point),
                data.tolerance.get().max(tol.confusion()),
            )),
            _ => None,
        }
    };
    let (Some(a), Some(b)) = (vertex(&first), vertex(&last)) else {
        return Ok(Vec::new());
    };
    let mut out = Vec::with_capacity(2);
    for (u, sign) in [(s, 1.0), (t, -1.0)] {
        let at = curve.point_at(u, tol)?;
        let d = curve.d1_at(u, tol)? * sign;
        let length = d.magnitude();
        if length <= tol.confusion() {
            continue;
        }
        // The vertex this end stands at: the nearer of the two.
        let (vertex, _, reach) = if at.distance(a.1) <= at.distance(b.1) {
            a
        } else {
            b
        };
        out.push(End {
            vertex,
            at,
            leaving: d / length,
            reach,
            fades: false,
        });
    }
    Ok(out)
}
