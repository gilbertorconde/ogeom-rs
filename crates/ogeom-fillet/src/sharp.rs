//! A setback corner with one edge left sharp: the bands of the vertex's
//! other edges stopped at their setbacks, and the opening between them
//! closed by a patch that runs into a point on the sharp edge.
//!
//! The opening crosses the sharp edge, and where it does, at `Q`, the
//! patch is to be tangent to both of that edge's faces: two tangent planes
//! at one point, which no patch with a tangent plane there has. The patch
//! at `Q` is three-sided, its fourth side shrunk to `Q`, and its normal
//! turns along that shrunk side from one face's to the other's, as a
//! cone's turns round its apex: tangent to both faces up to `Q` and,
//! locally, a cone at `Q` fanning between them.
//!
//! The sharp edge is read as a band of no width: its end on its section
//! plane is the point `Q`, and the trims on its two faces run from `Q`
//! into the bands' rails, each heading from `Q` for its band's rail half
//! way back to the vertex, so the two trims at `Q` part and the cone there
//! opens.
//!
//! The patch at `Q` runs out to the middle of each of those trims, where a
//! curve between them bounds it; it is fitted by
//! [`ogeom_algo::fit_network`], its halves of the trims held to them and
//! held across the faces, so it is tangent to both faces exactly up to the
//! held curves' own projection. The rest of the opening, six sides with
//! that curve, is filled by one patch G1 to every side, the patch at `Q`
//! among its supports.

use ogeom_algo::{Across, Built, CurveDef, History, NetworkOptions, PatchDef, PatchSide};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Continuity, Curve, Curve3d as _};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::fill_n::{FillBoundary, make_filling_n};
use crate::setback::{
    CornerEnd, Setbacks, bevel_end, corner_edges, filler_inside, point_of, round_end, trim_conic,
};

/// How many points along the curve between the two patches.
const ALONG: u32 = 32;
/// How many arcs, and points along each, seed the patch at `Q`.
const SEED: u32 = 40;
/// The lengths of each end tangent of an arc of the patch at `Q`, over the
/// arc's chord, tried in turn until the rest of the opening meets the
/// bands' target: fuller than a circle's (1.17 for a right angle), so the
/// arc leaves each trim nearly straight and the patch's normal turns
/// little where the rest meets it, and no fuller than keeps the patch's
/// normal at `Q` turning one way. On the block's corner at ten even and
/// uneven setbacks and radii, 1.6 met the target in every case, the worst
/// at 0.75 of it, its normal at `Q` turning one way; 4/3 missed the target
/// once and 2 folded the normal back twice, by up to 0.02 rad, and tangents
/// as long as the chord missed by up to 1.7 times.
const FULL: [f64; 3] = [1.6, 4.0 / 3.0, 2.0];
/// How far out from `Q` the patch's normal is read as its limit there, and
/// at how many stations along the shrunk side.
const AT_Q: f64 = 1e-7;
const TURNS: u32 = 64;

/// Round the edges of `vertex` other than `sharp` with bands at `radius`,
/// each stopped at its setback, and close the opening between them with a
/// patch that runs into the point `sharp`'s setback names on it.
///
/// The vertex is one of three planar faces and three straight edges, a
/// convex corner one ball touches every face of (a box's corner). The two
/// bands are cut square to their edges at their setbacks; the sharp edge
/// is cut at its own setback, at `Q`. On the face between the two bands
/// the trim is the conic tangent to both rails, as
/// [`round_vertex_setback`](crate::round_vertex_setback) trims; on each face
/// of the sharp edge the trim is the conic from `Q` tangent to the band's
/// rail, heading from `Q` for that rail half way back to the vertex.
///
/// The opening is closed by two patches. At `Q`, a three-sided one bounded
/// by the halves of the two trims from `Q` and a curve between their
/// middles, its fourth side shrunk to `Q`: fitted by
/// [`ogeom_algo::fit_network`] to arcs swept from `Q` that leave each trim
/// square across it in its face, its halves of the trims held to them and
/// held across the two faces, so it is tangent to each face exactly up to
/// the held curves' own projection, its normal turning at `Q` from one
/// face's to the other's. Over the rest, one patch from
/// [`make_filling_n`](crate::make_filling_n) meeting the two bands' end
/// arcs, the trims and the patch at `Q` at G1, to the bands' fit target.
/// Where it misses, or the patch at `Q` folds its normal back at `Q`, the
/// arcs are drawn fuller or flatter and the two patches fitted again.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// `radius` or a setback is not a positive length, a setback is shorter
/// than the ball's own along a rounded edge or reaches past an edge's far
/// end, the sharp edge's setback lies inside the bands' own corner, the
/// continuity is not G1, `sharp` does not leave the vertex, or the vertex
/// is not a convex planar corner of three edges one ball touches every
/// face of; [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where,
/// however its arcs are drawn, the rest misses the bands' fit target or
/// the patch at `Q` folds its normal back.
pub fn round_vertex_setback_sharp(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    sharp: &Shape,
    radius: f64,
    setbacks: &Setbacks,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut own: Vec<Shape> = vec![vertex.clone(), sharp.clone()];
    own.extend(setbacks.per_edge.iter().map(|(e, _)| e.clone()));
    ogeom_algo::on_own_nodes(model, solid, &own, |model, solid, parts| {
        let named: Vec<(Shape, f64)> = parts[2..]
            .iter()
            .zip(&setbacks.per_edge)
            .map(|(e, (_, s))| (e.clone(), *s))
            .collect();
        let own = Setbacks {
            per_edge: named,
            default: setbacks.default,
            continuity: setbacks.continuity,
        };
        sharp_own(model, solid, &parts[0], &parts[1], radius, &own, tol)
    })
}

/// A curve of the opening, as the patches take it: its curve over a range,
/// run the loop's way.
struct Side {
    curve: Curve,
    range: (f64, f64),
}

impl Side {
    fn at(&self, f: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (lo, hi) = self.range;
        self.curve.point_at((hi - lo).mul_add(f, lo), tol)
    }

    /// The unit tangent a share `f` along, the loop's way.
    fn tangent(&self, f: f64, tol: Tolerances) -> OgeomResult<Vector> {
        let (lo, hi) = self.range;
        let d = self.curve.d1_at((hi - lo).mul_add(f, lo), tol)? * (hi - lo).signum();
        d.normalized(tol)
    }

    fn def(&self) -> CurveDef {
        CurveDef::Curve {
            curve: self.curve.clone(),
            range: self.range,
        }
    }
}

/// [`round_vertex_setback_sharp`] on a shape every node below which is its
/// own.
#[allow(clippy::too_many_lines, reason = "one corner, assembled end to end")]
fn sharp_own(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    sharp: &Shape,
    radius: f64,
    setbacks: &Setbacks,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a blend radius must be a positive distance");
    }
    if setbacks.continuity != Continuity::G1 {
        ogeom_bail!(
            Construction,
            "a setback corner with an edge left sharp meets its bands and faces \
             at G1, not {:?}",
            setbacks.continuity
        );
    }
    let frame = crate::corner::corner_frame(model, solid, vertex, radius, tol)?;
    let n = frame.d.len();
    let (edges, reach) = corner_edges(model, solid, vertex, &frame, tol)?;
    let Some(ks) = edges.iter().position(|e| e.is_same(sharp)) else {
        ogeom_bail!(
            Construction,
            "the edge left sharp does not leave the vertex"
        );
    };
    if n != 3 {
        ogeom_bail!(
            Construction,
            "a setback corner with an edge left sharp closes a vertex of three \
             edges; one of {n} is still owed; see docs/CAPABILITIES.md, \
             fillet.setback-corners"
        );
    }
    let setback: Vec<f64> = edges
        .iter()
        .map(|e| {
            setbacks
                .per_edge
                .iter()
                .find(|(named, _)| named.is_same(e))
                .map_or(setbacks.default, |(_, s)| *s)
        })
        .collect();
    for k in 0..n {
        let s = setback[k];
        if !s.is_finite() || s <= tol.confusion() {
            ogeom_bail!(Construction, "a setback of {s} is not a positive distance");
        }
        let own = (frame.far - frame.corner).dot(frame.d[k].vector());
        if k == ks && s <= own + tol.confusion() {
            ogeom_bail!(
                Construction,
                "the edge left sharp is cut at {s}, inside the bands' own \
                 corner, which reaches {own} along it; cut it farther from \
                 the vertex"
            );
        }
        if k != ks && s < own - tol.confusion() {
            ogeom_bail!(
                Construction,
                "a setback of {s} stops the band inside the ball's own corner, \
                 which reaches {own} along its edge"
            );
        }
        if s >= reach[k] - tol.confusion() {
            ogeom_bail!(
                Construction,
                "a setback of {s} reaches past its edge, {} long",
                reach[k]
            );
        }
    }
    let rounded: Vec<Shape> = (0..n)
        .filter(|&k| k != ks)
        .map(|k| edges[k].clone())
        .collect();
    let banded = crate::fillet::fillet_edges(model, solid, &rounded, radius, tol)?;
    let kinds = vec![CornerEnd::Round; n];
    close_sharp_corner(model, banded, vertex, &frame, &setback, ks, &kinds, tol)
}

/// The bands of `banded` meeting at `vertex`, whose corner reads as
/// `frame`, each stopped at its `setback`, edge `ks` left sharp and cut at
/// its own, and the opening closed by the patch at `Q` and one patch over
/// the rest, as [`round_vertex_setback_sharp`] closes it. Each other edge
/// ends as `kinds` says: a round's end arc, or a bevel's straight end, which
/// the rest meets in position only, as it meets the trims from that end:
/// the bevel meets those faces in a crease. The patch at `Q` is held
/// across both faces of the sharp edge whatever the bands are, so its
/// normal turns there from one face's to the other's.
#[allow(clippy::too_many_lines, reason = "one corner, assembled end to end")]
#[allow(clippy::too_many_arguments, reason = "one corner and its reading")]
pub(crate) fn close_sharp_corner(
    model: &mut Model,
    banded: Built,
    vertex: &Shape,
    frame: &crate::corner::CornerFrame,
    setback: &[f64],
    ks: usize,
    kinds: &[CornerEnd],
    tol: Tolerances,
) -> OgeomResult<Built> {
    let n = frame.d.len();
    // The cell the patches replace: in from every plane through the vertex
    // and short of each edge's section, the sharp edge's included.
    let mut walls: Vec<(Vector, Point)> = frame.m.iter().map(|mk| (*mk, frame.corner)).collect();
    let sections: Vec<(Vector, Point)> = (0..n)
        .map(|k| {
            let dk = frame.d[k].vector();
            (-dk, frame.corner + dk * setback[k])
        })
        .collect();
    walls.extend(sections.iter().copied());
    let cell = crate::corner::convex_block(model, &walls, tol)?;
    let cut = ogeom_bool::cut(model, &banded.shape, &cell, tol)?;
    let history = banded.history.then(&cut.history);
    let near = tol.confusion() * 100.0;
    let on_plane = |p: Point, (normal, at): (Vector, Point)| normal.dot(p - at).abs() <= near;
    let mut kept: Vec<Shape> = Vec::new();
    for face in explore_unique(model, &cut.shape, ShapeType::Face)? {
        let data = model.node(&face).and_then(|n| n.data().as_face());
        let planar = data
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Plane(_)));
        let on_section = planar && {
            let (origin, outward) = ogeom_algo::face_normal(model, &face, tol)?;
            sections.iter().any(|section| {
                outward.cross(section.0).magnitude() <= tol.angular() * 10.0
                    && on_plane(origin, *section)
            })
        };
        if !on_section {
            kept.push(face);
        }
    }

    // `Q`, where the cut leaves the sharp edge, and each band's end arc.
    let q_point = sections[ks].1;
    let Some(vq) = explore_unique(model, &cut.shape, ShapeType::Vertex)?
        .into_iter()
        .find(|v| point_of(model, v).is_ok_and(|p| p.distance(q_point) <= near))
    else {
        ogeom_bail!(
            Construction,
            "the cut does not leave the sharp edge at its setback"
        );
    };
    let mut arcs: Vec<Option<(Shape, Shape)>> = Vec::with_capacity(n);
    for (k, section) in sections.iter().enumerate() {
        arcs.push(if k == ks {
            None
        } else if kinds[k] == CornerEnd::Bevel {
            Some(bevel_end(model, &kept, frame, *section, &walls, tol)?)
        } else {
            Some(round_end(model, &kept, k, *section, &walls, tol)?)
        });
    }
    let host = |t: usize| -> (Vector, Point) { (frame.m[frame.ring[t]], frame.corner) };
    let end_on = |model: &Model, k: usize, plane: (Vector, Point)| -> OgeomResult<(Shape, Point)> {
        let Some((arc, _)) = &arcs[k] else {
            return Ok((vq.clone(), q_point));
        };
        let Some((a, b)) = ogeom_algo::edge_vertices(model, arc)? else {
            ogeom_bail!(Construction, "a band's end arc has no vertices");
        };
        for v in [a, b] {
            let p = point_of(model, &v)?;
            if on_plane(p, plane) {
                return Ok((v, p));
            }
        }
        ogeom_bail!(Construction, "a band's end arc does not reach its face");
    };

    // The trim on the face after each edge, from its end there to the next
    // edge's: tangent to a band's rail at its end, and from `Q` heading for
    // the band's rail half way back to the vertex. A trim from `Q` is split
    // at its middle, where the patch at `Q` ends.
    let mut trims: Vec<Trim> = Vec::with_capacity(n);
    let mut fillers: Vec<Shape> = Vec::with_capacity(n);
    for k in 0..n {
        let next = (k + 1) % n;
        let plane = host(next);
        let (v0, p0) = end_on(model, k, plane)?;
        let (v1, p1) = end_on(model, next, plane)?;
        let Some(q) = crate::corner::planes_meet(&plane, &sections[k], &sections[next]) else {
            ogeom_bail!(
                Construction,
                "two section planes do not meet on the face between them"
            );
        };
        let into = |from: (usize, Point), to: (usize, Point)| -> OgeomResult<Vector> {
            if from.0 == ks {
                let aim = to.1 - frame.d[to.0].vector() * (setback[to.0] / 2.0);
                Ok((aim - from.1).normalized(tol)?)
            } else {
                Ok(-frame.d[from.0].vector())
            }
        };
        let (conic, apex) = trim_conic(
            (p0, into((k, p0), (next, p1))?),
            (p1, into((next, p1), (k, p0))?),
            tol,
        )?;
        let pieces: &[(f64, f64)] = if k == ks || next == ks {
            &[(0.0, 0.5), (0.5, 1.0)]
        } else {
            &[(0.0, 1.0)]
        };
        let mut ends = vec![v0.clone()];
        if pieces.len() == 2 {
            ends.push(ogeom_algo::make_vertex(model, conic.point_at(0.5, tol)?).shape);
        }
        ends.push(v1.clone());
        let mut edges = Vec::with_capacity(pieces.len());
        for (piece, range) in pieces.iter().enumerate() {
            edges.push(
                ogeom_algo::make_edge_between(
                    model,
                    conic.clone(),
                    *range,
                    &ends[piece],
                    &ends[piece + 1],
                    tol,
                )?
                .shape,
            );
        }
        let filler = filler_inside(
            model,
            edges.clone(),
            apex,
            (&v0, p0),
            (&v1, p1),
            q,
            plane,
            tol,
        )?;
        trims.push(Trim {
            curve: conic,
            edges,
            normal: plane.0,
            filler: filler.clone(),
        });
        fillers.push(filler);
    }
    let (r1, r2) = ((ks + 1) % n, (ks + 2) % n);

    // The patch at `Q`: the first half of the trim from `Q`, a curve from
    // its middle to the middle of the trim back to `Q`, and that trim's
    // second half, its fourth side shrunk to `Q`.
    let out_half = Side {
        curve: trims[ks].curve.clone(),
        range: (0.0, 0.5),
    };
    let back_half = Side {
        curve: trims[r2].curve.clone(),
        range: (0.5, 1.0),
    };
    let across = |normal: Vector, along: Vector, towards: Vector| -> OgeomResult<Vector> {
        let d = normal.cross(along).normalized(tol)?;
        Ok(if d.dot(towards) >= 0.0 { d } else { -d })
    };
    let (Some(arc1), Some(arc2)) = (&arcs[r1], &arcs[r2]) else {
        ogeom_bail!(Construction, "a rounded edge has no end arc");
    };
    let fit = crate::marched::band_fit_target(tol);
    // The patch at `Q` seeded on arcs of fullness `full`, and the rest of
    // the opening filled against it.
    let attempt = |model: &mut Model, full: f64| -> OgeomResult<(Shape, Shape)> {
        // Its arc a share `d` of the way out from `Q` along both halves:
        // the cubic from one to the other leaving each square across it in
        // its face, so the patch seeded on them is tangent to both faces.
        // The last is the curve between the patch and the rest.
        let arc = |d: f64, w: f64| -> OgeomResult<Point> {
            let (from, to) = (out_half.at(d, tol)?, back_half.at(1.0 - d, tol)?);
            let leave = across(trims[ks].normal, out_half.tangent(d, tol)?, to - from)?;
            let arrive = -across(
                trims[r2].normal,
                back_half.tangent(1.0 - d, tol)?,
                from - to,
            )?;
            let reach = from.distance(to) * full;
            let (w2, w3) = (w * w, w * w * w);
            let (h00, h10) = (2.0 * w3 - 3.0 * w2 + 1.0, w3 - 2.0 * w2 + w);
            let (h01, h11) = (-2.0 * w3 + 3.0 * w2, w3 - w2);
            Ok(Point::ORIGIN
                + (from - Point::ORIGIN) * h00
                + leave * (reach * h10)
                + (to - Point::ORIGIN) * h01
                + arrive * (reach * h11))
        };
        let between: Vec<Point> = (0..=ALONG)
            .map(|i| arc(1.0, f64::from(i) / f64::from(ALONG)))
            .collect::<OgeomResult<_>>()?;
        let mut cloud: Vec<Point> = Vec::new();
        for i in 1..=SEED {
            let d = f64::from(i) / f64::from(SEED);
            for j in 0..=SEED {
                cloud.push(arc(d, f64::from(j) / f64::from(SEED))?);
            }
        }
        let index = ogeom_mesh::PointIndex::build(&cloud)?;
        let curves = vec![out_half.def(), CurveDef::Polyline(between), back_half.def()];
        let side = |curve: usize| {
            Some(PatchSide {
                curve,
                reversed: false,
            })
        };
        let patches = vec![PatchDef {
            sides: [side(0), side(1), side(2), None],
        }];
        // The halves are held to the trims and held across the faces, which
        // makes the patch tangent to both up to `Q`; the points only shape
        // it inside, so how near it comes to them is not asked.
        let options = NetworkOptions {
            tolerance: fit,
            smoothing: fit * fit,
            creases: Vec::new(),
            spans: 2,
            most: 16,
            held: vec![0, 2],
            across: vec![
                (0, Across::Plane(trims[ks].normal)),
                (2, Across::Plane(trims[r2].normal)),
            ],
        };
        let network = ogeom_algo::fit_network(model, &patches, &curves, &index, &options, tol)?;
        // Its normal at `Q`, the limit along the shrunk side from one face's
        // to the other's, turns one way: the patch is a cone there that
        // does not fold back on itself.
        let surface = &network.surfaces[0];
        let limit = |w: f64| -> OgeomResult<Vector> {
            let (du, dw) = ogeom_geom::Surface::d1_at(surface, AT_Q, w, tol)?;
            du.cross(dw).normalized(tol)
        };
        let first = limit(0.0)?;
        let mut turned = 0.0_f64;
        for i in 1..=TURNS {
            let now = limit(f64::from(i) / f64::from(TURNS))?
                .dot(first)
                .clamp(-1.0, 1.0)
                .acos();
            if now < turned - fit {
                ogeom_bail!(
                    NotDone,
                    "the patch at the sharp edge turns its normal back by {} there",
                    turned - now
                );
            }
            turned = turned.max(now);
        }
        let at_q = network.faces[0].clone();
        // The rest of the opening: the trim from `Q`'s second half, band
        // `r1`'s arc, the trim between the bands, band `r2`'s arc, the trim
        // back to `Q`'s first half, and the curve along the patch at `Q`,
        // filled by one patch G1 to every side.
        // G1 to a round's end and to a trim between rounds; a bevel's end,
        // and a trim from it, in position only.
        let side = |edge: &Shape, support: &Shape, tangent: bool| FillBoundary {
            edge: edge.clone(),
            support: Some(support.clone()),
            continuity: if tangent {
                Continuity::G1
            } else {
                Continuity::C0
            },
        };
        let round = |k: usize| kinds[k] == CornerEnd::Round;
        let boundary = vec![
            side(&trims[ks].edges[1], &trims[ks].filler, round(r1)),
            side(&arc1.0, &arc1.1, round(r1)),
            side(
                &trims[r1].edges[0],
                &trims[r1].filler,
                round(r1) && round(r2),
            ),
            side(&arc2.0, &arc2.1, round(r2)),
            side(&trims[r2].edges[0], &trims[r2].filler, round(r2)),
            side(&network.edges[1], &at_q, true),
        ];
        let filled = make_filling_n(model, &boundary, &[], fit, tol)?;
        Ok((at_q, filled.built.shape))
    };
    let mut made = None;
    for (round, full) in FULL.iter().enumerate() {
        match attempt(model, *full) {
            Ok(both) => {
                made = Some(both);
                break;
            }
            Err(ogeom_core::OgeomError::NotDone(_)) if round + 1 < FULL.len() => {}
            Err(e) => return Err(e),
        }
    }
    let Some((at_q, patch)) = made else {
        ogeom_bail!(NotDone, "no fullness was tried");
    };
    let mut faces = kept;
    faces.extend(fillers);
    faces.push(at_q.clone());
    faces.push(patch.clone());
    let sewn = ogeom_algo::sew_within(model, &faces, fit, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(
            Construction,
            "the setback corner's faces did not close: {} shells",
            sewn.shells.len()
        );
    }
    let rounded = model.add_solid(std::slice::from_ref(&sewn.shells[0]))?;
    // The two patches as the solid holds them: sewing copies a face it
    // shares with another shape, the patch at `Q` its network's shell.
    let surface_of = |model: &Model, face: &Shape| {
        model
            .node(face)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface)
    };
    let mut history: History = history;
    for made in [&at_q, &patch] {
        let wanted = surface_of(model, made);
        let Some(face) = explore_unique(model, &rounded, ShapeType::Face)?
            .into_iter()
            .find(|f| wanted.is_some() && surface_of(model, f) == wanted)
        else {
            ogeom_bail!(
                Construction,
                "a patch of the corner is not in the sewn solid"
            );
        };
        history.generate(vertex, face);
    }
    Ok(Built::new(rounded, history))
}

/// A face's trim: its conic, run from the end on the edge before it to the
/// end on the edge after, the edges it is built as, the face's normal and
/// its filler inside the cell.
struct Trim {
    curve: Curve,
    edges: Vec<Shape>,
    normal: Vector,
    filler: Shape,
}
