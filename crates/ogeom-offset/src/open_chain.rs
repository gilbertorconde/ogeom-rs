//! An open planar chain thickened into a band, and a rib grown from one.
//!
//! The band runs between two offsets of the chain, both signed to the
//! chain's left about the plane's normal: `(2, 0)` keeps the chain as the
//! right side and puts the left side 2 to its left, `(1, -1)` straddles it.
//! Each side is the chain's offset at its own distance, its corners closed
//! as a closed wire's offset closes them (a gap by the join, an overlap
//! trimmed to the crossing), and the two ends are closed by a cap.
//!
//! A rib is that band (or the chain thickened across its plane) grown from
//! the chain until the solid stops it on every side, and fused on. The
//! chain's ends are first carried along their tangents into the walls they
//! point at; the band then sweeps away from the chain, and the piece of the
//! sweep outside the solid that holds the chain is the web. A side of the
//! chain whose sweep reaches its far end unstopped meets no wall there.

use ogeom_algo::{Built, History, SolidBoundary, edge_vertices, make_prism};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Frame, Plane, Point, Point2, Transform, Vector, Vector2};
use ogeom_topo::{Filter, Model, Shape, ShapeType, explore};

use crate::wire2d::{
    Join, Piece, at_angle, build_ring, crossings, extended_to_meeting, lies_inside,
    nearest_crossing, read_pieces, same_support, trim_end, trim_start, within,
};

/// How an open chain's band is closed at the chain's two ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Cap {
    /// A straight line square to the chain through its end point.
    #[default]
    Flat,
    /// A half circle across the band's width, bulging past the end.
    Round,
    /// Both sides carried straight on past the end by half the band's
    /// width, then joined square across.
    Extended,
}

/// Where a piece of the band came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Origin {
    /// The offset of the chain's edge `i`.
    Side(usize),
    /// The join at the vertex after the chain's edge `i`.
    Corner(usize),
    /// A cap: at the chain's end when `true`, at its start otherwise.
    Cap(bool),
}

/// The chain read into the plane: the frame (right-handed about the
/// plane's normal), the edges in order, and their pieces.
struct Chain {
    frame: Frame,
    edges: Vec<Shape>,
    pieces: Vec<Piece>,
}

impl Chain {
    fn lift(&self, p: Point2) -> Point {
        self.frame.origin() + self.frame.x().vector() * p.x + self.frame.y().vector() * p.y
    }

    fn lift_vector(&self, v: Vector2) -> Vector {
        self.frame.x().vector() * v.x + self.frame.y().vector() * v.y
    }
}

/// Read an open chain lying in `plane`.
fn read_chain(model: &Model, wire: &Shape, plane: Plane, tol: Tolerances) -> OgeomResult<Chain> {
    if model.kind_of(wire)? != ShapeType::Wire {
        ogeom_bail!(Construction, "thickening an open chain starts from a wire");
    }
    if ogeom_algo::is_wire_closed(model, wire, tol)? {
        ogeom_bail!(
            Construction,
            "the chain is closed and has no ends to cap; offset it with offset_wire"
        );
    }
    let normal = plane.frame().z();
    let frame = Frame::new(plane.frame().origin(), normal, plane.frame().x(), tol)?;
    let edges = explore(model, wire, Filter::OfType(ShapeType::Edge))?;
    if edges.is_empty() {
        ogeom_bail!(Construction, "the chain has no edges");
    }
    for edge in &edges {
        let Some((curve, range)) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
            ogeom_bail!(Construction, "an edge of the chain has no curve");
        };
        for k in 0..=4 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 4.0;
            let p = curve.point_at(t, tol)?;
            if (p - frame.origin()).dot(normal.vector()).abs() > tol.confusion() * 10.0 {
                ogeom_bail!(
                    Construction,
                    "the chain leaves its plane: {p:?} stands off it"
                );
            }
        }
    }
    let flat = |p: Point| {
        let v = p - frame.origin();
        Point2::new(v.dot(frame.x().vector()), v.dot(frame.y().vector()))
    };
    let pieces = read_pieces(model, &edges, &flat, tol)?;
    for (i, pair) in pieces.windows(2).enumerate() {
        if pair[0].end_point().distance(pair[1].start_point()) > tol.confusion() * 10.0 {
            ogeom_bail!(
                Construction,
                "the chain's edges {i} and {} do not follow one another",
                i + 1
            );
        }
        let (a, b) = (pair[0].tangent(true), pair[1].tangent(false));
        if a.cross(b).abs() <= 1e-9 && a.dot(b) < 0.0 {
            ogeom_bail!(
                Construction,
                "the chain turns straight back on itself after edge {i}; \
                 its band has no side there"
            );
        }
    }
    Ok(Chain {
        frame,
        edges,
        pieces,
    })
}

/// The chain's left at a tangent: the tangent turned a quarter turn
/// counter-clockwise about the plane's normal.
fn left_of(tangent: Vector2) -> Vector2 {
    Vector2::new(-tangent.y, tangent.x)
}

/// One side of the band: the chain offset `o` to its left, corners closed.
fn side(
    pieces: &[Piece],
    o: f64,
    join: Join,
    tol: Tolerances,
) -> OgeomResult<Vec<(Piece, Origin)>> {
    let mut out: Vec<(Piece, Origin)> = Vec::with_capacity(pieces.len() * 2);
    let still = o.abs() <= tol.confusion();
    for (i, piece) in pieces.iter().enumerate() {
        let moved = if still {
            piece.clone()
        } else {
            match piece {
                Piece::Seg { from, to } => {
                    let shift = left_of(piece.tangent(false)) * o;
                    Piece::Seg {
                        from: *from + shift,
                        to: *to + shift,
                    }
                }
                Piece::Arc {
                    centre,
                    radius,
                    start,
                    end,
                } => {
                    // Whether the chain's left is radially out depends on
                    // the arc's sweep: out for clockwise, in for
                    // counter-clockwise.
                    let outward = if end > start { -1.0 } else { 1.0 };
                    let grown = radius + o * outward;
                    if grown <= tol.confusion() {
                        ogeom_bail!(
                            Construction,
                            "the band's side {o} from the chain reaches past the \
                             centre of arc {i} (radius {radius}); the thickness \
                             on that side must stay under the radius"
                        );
                    }
                    Piece::Arc {
                        centre: *centre,
                        radius: grown,
                        start: *start,
                        end: *end,
                    }
                }
            }
        };
        out.push((moved, Origin::Side(i)));
    }
    if still {
        return Ok(out);
    }
    for i in 0..pieces.len() - 1 {
        let Some(at_i) = out.iter().position(|(_, g)| *g == Origin::Side(i)) else {
            continue;
        };
        let at_j = at_i + 1;
        let e = out[at_i].0.end_point();
        let s = out[at_j].0.start_point();
        if e.distance(s) <= tol.confusion() * 10.0 {
            continue; // Tangent-continuous: nothing to do.
        }
        let turn = pieces[i].tangent(true).cross(pieces[i + 1].tangent(false));
        let corner = pieces[i].end_point();
        if turn * o < 0.0 {
            // A side on the outside of the turn: the pieces separate.
            match join {
                Join::Arc => {
                    let a0 = (e - corner).y.atan2((e - corner).x);
                    let mut a1 = (s - corner).y.atan2((s - corner).x);
                    let tau = core::f64::consts::TAU;
                    while a1 - a0 > core::f64::consts::PI {
                        a1 -= tau;
                    }
                    while a0 - a1 > core::f64::consts::PI {
                        a1 += tau;
                    }
                    out.insert(
                        at_j,
                        (
                            Piece::Arc {
                                centre: corner,
                                radius: o.abs(),
                                start: a0,
                                end: a1,
                            },
                            Origin::Corner(i),
                        ),
                    );
                }
                Join::Intersection => {
                    let (a, b) = extended_to_meeting(&out[at_i].0, &out[at_j].0, corner, tol)?;
                    out[at_i].0 = a;
                    out[at_j].0 = b;
                }
            }
        } else {
            // On the inside: the pieces overlap and are trimmed to their
            // crossing.
            let met = nearest_crossing(&out[at_i].0, &out[at_j].0, corner, tol)?;
            trim_end(&mut out[at_i].0, met, tol)?;
            trim_start(&mut out[at_j].0, met, tol)?;
        }
    }
    Ok(out)
}

fn reversed(piece: &Piece) -> Piece {
    match piece {
        Piece::Seg { from, to } => Piece::Seg {
            from: *to,
            to: *from,
        },
        Piece::Arc {
            centre,
            radius,
            start,
            end,
        } => Piece::Arc {
            centre: *centre,
            radius: *radius,
            start: *end,
            end: *start,
        },
    }
}

/// The cap from `a` to `b` across an end whose outward direction is `out`:
/// `a` on the band's left, `b` on its right as seen walking out.
fn cap_pieces(a: Point2, b: Point2, out: Vector2, cap: Cap, at_end: bool) -> Vec<(Piece, Origin)> {
    let origin = Origin::Cap(at_end);
    match cap {
        Cap::Flat => vec![(Piece::Seg { from: a, to: b }, origin)],
        Cap::Round => {
            let centre = a.midpoint(b);
            let start = (a - centre).y.atan2((a - centre).x);
            // Left to right through `out` is clockwise.
            vec![(
                Piece::Arc {
                    centre,
                    radius: a.distance(b) / 2.0,
                    start,
                    end: start - core::f64::consts::PI,
                },
                origin,
            )]
        }
        Cap::Extended => {
            let reach = out * (a.distance(b) / 2.0);
            vec![
                (
                    Piece::Seg {
                        from: a,
                        to: a + reach,
                    },
                    origin,
                ),
                (
                    Piece::Seg {
                        from: a + reach,
                        to: b + reach,
                    },
                    origin,
                ),
                (
                    Piece::Seg {
                        from: b + reach,
                        to: b,
                    },
                    origin,
                ),
            ]
        }
    }
}

/// The band's closed outline in the plane, walked clockwise: the left side
/// forward, the end cap, the right side back, the start cap.
fn band_ring(
    pieces: &[Piece],
    hi: f64,
    lo: f64,
    cap: Cap,
    join: Join,
    tol: Tolerances,
) -> OgeomResult<Vec<(Piece, Origin)>> {
    let left = side(pieces, hi, join, tol)?;
    let right = side(pieces, lo, join, tol)?;
    let (Some(first), Some(last)) = (pieces.first(), pieces.last()) else {
        ogeom_bail!(Construction, "the chain has no edges");
    };
    let (Some((l_first, _)), Some((l_last, _))) = (left.first(), left.last()) else {
        ogeom_bail!(Construction, "the band's side is empty");
    };
    let (Some((r_first, _)), Some((r_last, _))) = (right.first(), right.last()) else {
        ogeom_bail!(Construction, "the band's side is empty");
    };
    let end_cap = cap_pieces(
        l_last.end_point(),
        r_last.end_point(),
        last.tangent(true),
        cap,
        true,
    );
    let start_cap = cap_pieces(
        r_first.start_point(),
        l_first.start_point(),
        -first.tangent(false),
        cap,
        false,
    );
    let mut ring = left;
    ring.extend(end_cap);
    ring.extend(right.iter().rev().map(|(p, g)| (reversed(p), *g)));
    ring.extend(start_cap);

    // The band is refused where its outline meets itself: a thickness too
    // great for the chain's turns, or a chain that crosses itself.
    let n = ring.len();
    for i in 0..n {
        for j in i + 1..n {
            let (a, b) = (&ring[i].0, &ring[j].0);
            let crossed = if same_support(a, b, tol) {
                [b.start_point(), b.end_point()]
                    .iter()
                    .any(|p| lies_inside(a, *p, tol))
                    || [a.start_point(), a.end_point()]
                        .iter()
                        .any(|p| lies_inside(b, *p, tol))
            } else {
                crossings(a, b, tol)?
                    .into_iter()
                    .any(|p| within(a, p, tol) && within(b, p, tol))
            };
            if crossed {
                ogeom_bail!(
                    Construction,
                    "the band's outline crosses itself: the thickness is too \
                     great for the chain's turns, or the chain crosses itself"
                );
            }
        }
    }
    Ok(ring)
}

/// The closed outline of an open planar chain of lines and arcs thickened
/// from `left` to `right`.
///
/// Both are signed offsets to the chain's left about `plane`'s normal; 0
/// keeps the chain itself as one side, and their order does not matter.
/// The ends are closed by `cap` and each side's corners by `join`: a side
/// on the outside of a turn is joined across the gap, one on the inside is
/// trimmed to where its pieces cross. An arc's side is the concentric arc.
///
/// The history reports each chain edge generating its two sides, each
/// inner vertex generating its joins, each end vertex generating its cap,
/// and the chain generating the outline.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the wire is closed or leaves `plane`, an edge is neither straight nor
/// circular, the band has no width, a side reaches past an arc's centre, an
/// `Intersection` join finds no meeting of curved sides, a trim consumes a
/// side's piece, the chain turns straight back on itself, or the outline
/// crosses itself.
#[allow(
    clippy::too_many_arguments,
    reason = "the call the issue fixes: a chain, a plane, two offsets, two corner rules"
)]
pub fn thicken_open_wire(
    model: &mut Model,
    wire: &Shape,
    plane: Plane,
    left: f64,
    right: f64,
    cap: Cap,
    join: Join,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !left.is_finite() || !right.is_finite() {
        ogeom_bail!(Construction, "a band from {left} to {right} is not finite");
    }
    let (hi, lo) = if left >= right {
        (left, right)
    } else {
        (right, left)
    };
    if hi - lo <= tol.confusion() {
        ogeom_bail!(Construction, "a band from {left} to {right} has no width");
    }
    let chain = read_chain(model, wire, plane, tol)?;
    let ring = band_ring(&chain.pieces, hi, lo, cap, join, tol)?;
    let (outline, edges) = build_ring(model, &chain.frame, &ring, tol)?;

    let mut history = History::new();
    let count = chain.edges.len();
    for ((_, origin), edge) in ring.iter().zip(edges) {
        match origin {
            Origin::Side(i) => history.generate(&chain.edges[*i], edge),
            Origin::Corner(i) => {
                if let Some((_, vertex)) = edge_vertices(model, &chain.edges[*i])? {
                    history.generate(&vertex, edge);
                }
            }
            Origin::Cap(at_end) => {
                let (index, end) = if *at_end {
                    (count - 1, true)
                } else {
                    (0, false)
                };
                if let Some((first, last)) = edge_vertices(model, &chain.edges[index])? {
                    history.generate(if end { &last } else { &first }, edge);
                }
            }
        }
    }
    history.generate(wire, outline.clone());
    Ok(Built::new(outline, history))
}

/// Where a ray from `from` along `along` first meets the solid, within
/// `reach`, and how far it then runs inside: `(entry, inside)`. `None`
/// when the ray stays clear.
fn ray_entry(
    model: &Model,
    boundary: &SolidBoundary,
    from: Point,
    along: Vector,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(f64, f64)>> {
    let at = |s: f64| from + along * s;
    let first_meeting = |start: f64| -> OgeomResult<Option<f64>> {
        if boundary.clear_from(model, at(start), at(reach), &[], tol)? {
            return Ok(None);
        }
        let (mut clear, mut met) = (start, reach);
        for _ in 0..60 {
            let mid = f64::midpoint(clear, met);
            if boundary.clear_from(model, at(start), at(mid), &[], tol)? {
                clear = mid;
            } else {
                met = mid;
            }
            if met - clear <= tol.confusion() {
                break;
            }
        }
        Ok(Some(met))
    };
    if boundary.holds(model, from, tol)? != ogeom_algo::Containment::Out {
        return Ok(Some((0.0, 0.0)));
    }
    let Some(entry) = first_meeting(0.0)? else {
        return Ok(None);
    };
    // Past the entry, the run inside ends at the next meeting; a step
    // just inside starts the search clear of the face entered.
    let step = (tol.confusion() * 1e5)
        .min(reach * 1e-6)
        .max(tol.confusion() * 100.0);
    let inside_start = entry + step;
    if boundary.holds(model, at(inside_start), tol)? != ogeom_algo::Containment::In {
        return Ok(Some((entry, 0.0)));
    }
    let exit = first_meeting(inside_start)?.unwrap_or(inside_start);
    Ok(Some((entry, exit - entry)))
}

/// A rib: the chain thickened, grown until the solid stops it on every
/// side, and fused onto the solid.
///
/// `chain` is an open chain of lines and arcs in `plane`. With
/// `along_normal` false the web lies in the plane: it is thickened across
/// the plane, `thickness.0` along the plane's normal and `thickness.1`
/// against it, and grows within the plane square to the chain's chord
/// (from its first point to its last). With `along_normal` true the chain
/// is thickened within the plane, `thickness.0` to its left and
/// `thickness.1` to its right, with flat ends and sharp corners, and the
/// band grows along the plane's normal.
///
/// Each end of the chain is first carried straight on along its tangent
/// into the wall it points at; an end that points at nothing stays where
/// it is, and the web ends flat there, square to the growth, so it never
/// runs round the end of a wall shorter than the chain. Of the two ways the
/// web can grow, the ones the solid closes are kept: a way whose growth
/// runs out past the solid meets no wall. Where both close, the smaller
/// web is the rib.
///
/// The history reports the chain generating the result.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `solid` is not a solid, a thickness is negative or both are zero, the
/// chain is refused as [`thicken_open_wire`] refuses it, an in-plane chain
/// runs back along its chord or has its ends together, the chain meets no
/// wall (grown either way it runs clear of the solid), or the chain lies
/// wholly inside the solid; and whatever the booleans refuse.
pub fn make_rib(
    model: &mut Model,
    solid: &Shape,
    chain: &Shape,
    plane: Plane,
    thickness: (f64, f64),
    along_normal: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (t0, t1) = thickness;
    if !t0.is_finite() || !t1.is_finite() || t0 < 0.0 || t1 < 0.0 {
        ogeom_bail!(
            Construction,
            "a rib's thicknesses ({t0}, {t1}) must be finite and not negative"
        );
    }
    if t0 + t1 <= tol.confusion() {
        ogeom_bail!(Construction, "a rib of no thickness holds nothing");
    }
    if model.kind_of(solid)? != ShapeType::Solid {
        ogeom_bail!(Construction, "a rib grows onto a solid");
    }
    let read = read_chain(model, chain, plane, tol)?;
    let normal = read.frame.z().vector();
    let bounds = ogeom_algo::shape_bounds(model, solid, tol)?;
    let (Some(low), Some(high)) = (bounds.low(), bounds.high()) else {
        ogeom_bail!(Construction, "the solid has no extent for a rib to meet");
    };
    let centre = low.midpoint(high);
    let mut far = low.distance(high);
    for piece in &read.pieces {
        for p in [piece.start_point(), piece.mid_point(), piece.end_point()] {
            far = far.max(read.lift(p).distance(centre) * 2.0);
        }
    }
    // Past everything: twice the solid's diagonal and the chain's distance.
    let reach = far.mul_add(2.0, t0 + t1 + 1.0);

    // Each end carried on into the wall it points at, halfway through the
    // wall so the web's end is the wall's own face once the solid is cut
    // away.
    let boundary = SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)?;
    let (Some(first), Some(last)) = (read.pieces.first(), read.pieces.last()) else {
        ogeom_bail!(Construction, "the chain has no edges");
    };
    let mut pieces = read.pieces.clone();
    let start = first.start_point();
    let back = -first.tangent(false);
    if let Some((entry, inside)) = ray_entry(
        model,
        &boundary,
        read.lift(start),
        read.lift_vector(back),
        reach,
        tol,
    )? {
        let to = entry + inside / 2.0;
        if to > tol.confusion() {
            pieces.insert(
                0,
                Piece::Seg {
                    from: start + back * to,
                    to: start,
                },
            );
        }
    }
    let end = last.end_point();
    let ahead = last.tangent(true);
    if let Some((entry, inside)) = ray_entry(
        model,
        &boundary,
        read.lift(end),
        read.lift_vector(ahead),
        reach,
        tol,
    )? {
        let to = entry + inside / 2.0;
        if to > tol.confusion() {
            pieces.push(Piece::Seg {
                from: end,
                to: end + ahead * to,
            });
        }
    }

    // The ways the web may grow, each with the band it sweeps.
    let mut ways: Vec<(Vector, Shape)> = Vec::with_capacity(2);
    if along_normal {
        let ring = band_ring(&pieces, t0, -t1, Cap::Flat, Join::Intersection, tol)?;
        let (outline, _) = build_ring(model, &read.frame, &ring, tol)?;
        let face = planar_face(model, &read, &outline, tol)?;
        for sign in [1.0, -1.0] {
            let d = normal * sign;
            let tool = make_prism(model, &face, d * reach, tol)?.shape;
            ways.push((d, tool));
        }
    } else {
        let (Some(head), Some(tail)) = (pieces.first(), pieces.last()) else {
            ogeom_bail!(Construction, "the chain has no edges");
        };
        let chord = tail.end_point() - head.start_point();
        if chord.magnitude() <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "an in-plane rib grows square to its chain's chord, and this \
                 chain's ends are together"
            );
        }
        let u = chord / chord.magnitude();
        for (i, piece) in pieces.iter().enumerate() {
            if !runs_along(piece, u) {
                ogeom_bail!(
                    Construction,
                    "an in-plane rib grows square to its chain's chord, and \
                     piece {i} runs back along it; the swept web would fold"
                );
            }
        }
        for sign in [1.0, -1.0] {
            let d2 = left_of(u) * sign;
            let shift = d2 * reach;
            let mut ring: Vec<(Piece, ())> = pieces.iter().map(|p| (p.clone(), ())).collect();
            let (a, b) = (head.start_point(), tail.end_point());
            ring.push((
                Piece::Seg {
                    from: b,
                    to: b + shift,
                },
                (),
            ));
            for p in pieces.iter().rev() {
                let moved = match reversed(p) {
                    Piece::Seg { from, to } => Piece::Seg {
                        from: from + shift,
                        to: to + shift,
                    },
                    Piece::Arc {
                        centre,
                        radius,
                        start,
                        end,
                    } => Piece::Arc {
                        centre: centre + shift,
                        radius,
                        start,
                        end,
                    },
                };
                ring.push((moved, ()));
            }
            ring.push((
                Piece::Seg {
                    from: a + shift,
                    to: a,
                },
                (),
            ));
            let (outline, _) = build_ring(model, &read.frame, &ring, tol)?;
            let face = planar_face(model, &read, &outline, tol)?;
            let lowered =
                ogeom_algo::transformed(model, &face, Transform::translation(-normal * t1))?.shape;
            let tool = make_prism(model, &lowered, normal * (t0 + t1), tol)?.shape;
            ways.push((read.lift_vector(d2), tool));
        }
    }

    // Probes just off the chain toward each way's growth, at the middle of
    // the web's thickness: the web is the piece of the cut that holds them.
    let nudge = (tol.confusion() * 1e4).max(reach * 1e-7);
    let mut webs: Vec<(f64, Vec<Shape>)> = Vec::new();
    let mut inside_everywhere = true;
    for (d, tool) in &ways {
        let mut probes: Vec<Point> = Vec::new();
        for piece in &read.pieces {
            for k in 1..4 {
                let p = point_along(piece, f64::from(k) / 4.0);
                let across = if along_normal {
                    let left = left_of(tangent_along(piece, f64::from(k) / 4.0));
                    read.lift_vector(left * ((t0 - t1) / 2.0))
                } else {
                    normal * ((t0 - t1) / 2.0)
                };
                probes.push(read.lift(p) + across + *d * nudge);
            }
        }
        let cut = ogeom_bool::cut(model, tool, solid, tol)?.shape;
        let pieces_of_cut = explore(model, &cut, Filter::OfType(ShapeType::Solid))?;
        let mut kept: Vec<Shape> = Vec::new();
        let mut leaks = false;
        for piece in pieces_of_cut {
            let piece_boundary = SolidBoundary::of(model, &piece, tol.confusion() * 1e4, tol)?;
            let mut holds = false;
            for p in &probes {
                if piece_boundary.holds(model, *p, tol)? == ogeom_algo::Containment::In {
                    holds = true;
                    break;
                }
            }
            if !holds {
                continue;
            }
            // A piece reaching the sweep's far end was never stopped.
            let base = read.lift(read.pieces[0].start_point());
            for v in explore(model, &piece, Filter::OfType(ShapeType::Vertex))? {
                if let Some(data) = model.node(&v).and_then(|n| n.data().as_vertex()) {
                    let at = v.transform(model.datums())?.apply(data.point);
                    if (at - base).dot(*d) > reach / 2.0 {
                        leaks = true;
                    }
                }
            }
            kept.push(piece);
        }
        if !kept.is_empty() {
            inside_everywhere = false;
        }
        if leaks || kept.is_empty() {
            continue;
        }
        let mut volume = 0.0;
        for piece in &kept {
            volume += ogeom_algo::volume_properties(
                model,
                piece,
                ogeom_mesh::Deflection::for_tolerances(tol),
                tol,
            )?
            .mass;
        }
        webs.push((volume, kept));
    }
    if webs.is_empty() {
        if inside_everywhere {
            ogeom_bail!(
                Construction,
                "the rib's chain lies inside the solid; there is no web to add"
            );
        }
        ogeom_bail!(
            Construction,
            "the rib meets no wall: grown either way from the chain, the web \
             runs clear of the solid"
        );
    }
    webs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let (_, web) = webs.remove(0);
    let mut result = ogeom_bool::fuse(model, solid, &web[0], tol)?;
    for more in &web[1..] {
        let next = ogeom_bool::fuse(model, &result.shape, more, tol)?;
        result = next;
    }
    result.history.generate(chain, result.shape.clone());
    Ok(result)
}

/// The face a planar outline bounds, in the chain's plane.
fn planar_face(
    model: &mut Model,
    read: &Chain,
    outline: &Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let plane = Plane::new(read.frame);
    let faces =
        ogeom_algo::faces_from_planar_wires(model, plane, std::slice::from_ref(outline), tol)?;
    faces
        .faces
        .into_iter()
        .next()
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "the band's outline bounds no face"))
}

/// The point a fraction `f` of the way along a piece.
fn point_along(piece: &Piece, f: f64) -> Point2 {
    match piece {
        Piece::Seg { from, to } => *from + (*to - *from) * f,
        Piece::Arc {
            centre,
            radius,
            start,
            end,
        } => at_angle(*centre, *radius, start + (end - start) * f),
    }
}

/// The unit tangent a fraction `f` of the way along a piece.
fn tangent_along(piece: &Piece, f: f64) -> Vector2 {
    match piece {
        Piece::Seg { .. } => piece.tangent(false),
        Piece::Arc {
            centre,
            radius,
            start,
            end,
        } => {
            let a = start + (end - start) * f;
            let radial = (at_angle(*centre, *radius, a) - *centre) / *radius;
            if end > start {
                Vector2::new(-radial.y, radial.x)
            } else {
                Vector2::new(radial.y, -radial.x)
            }
        }
    }
}

/// Whether the piece never runs back along `u`.
fn runs_along(piece: &Piece, u: Vector2) -> bool {
    (0..=16).all(|k| tangent_along(piece, f64::from(k) / 16.0).dot(u) >= -1e-9)
}
