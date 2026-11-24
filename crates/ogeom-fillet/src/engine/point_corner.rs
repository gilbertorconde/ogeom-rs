//! The two-patch G1 opening at the endpoint of an edge left sharp.

use ogeom_algo::{
    FillBoundary, FillRegion as Region, make_filling_n_region, make_filling_n_within,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_fit::{Across, CurveDef, NetworkOptions, PatchDef, PatchSide};
use ogeom_geom::{Continuity, Curve, Curve3d as _};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape};

/// How many points along the curve between the two patches.
const ALONG: u32 = 32;
/// How many arcs, and points along each, seed the patch at `Q`.
const SEED: u32 = 40;
/// The lengths of each end tangent of an arc of the patch at `Q`, over the
/// arc's chord, tried in order until the rest meets the bands' fit target.
/// Fuller than a circular arc, each leaves its trim nearly straight; the
/// resulting patch must still turn its normal monotonically at `Q`.
const FULL: [f64; 3] = [1.6, 4.0 / 3.0, 2.0];
/// How far out from `Q` the patch's normal is read as its limit there, and
/// at how many stations along the shrunk side.
const AT_Q: f64 = 1e-7;
const TURNS: u32 = 64;

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

/// Close an opening that runs into `Q`, a point on a sharp edge: `out` the
/// trim from `Q` to the first band's end, `middle` the trim from that band
/// to the second, `back` the trim from the second band's end to `Q`, each
/// built by [`trim_edges`], those from `Q` split at their middles; `arcs`
/// the two bands' end arcs and their faces. The patch at `Q`, three-sided,
/// its fourth side shrunk to `Q`, tangent to both faces of the sharp edge
/// exactly up to the held trims' projection, and the patch over the rest,
/// meeting the second half of `out`, the first arc, `middle`, the second
/// arc and the first half of `back` at G1 where `rest` says and in
/// position elsewhere, and the patch at `Q` at G1. Neither patch may stand
/// outside what `held` holds them in.
pub(crate) fn close_through_q(
    model: &mut Model,
    [out, middle, back]: [&Trim; 3],
    [arc1, arc2]: [&(Shape, Shape); 2],
    rest: [bool; 5],
    held: &Held<'_>,
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape)> {
    let Held {
        within,
        seam,
        region,
    } = *held;
    // The patch at `Q`: the first half of the trim from `Q`, a curve from
    // its middle to the middle of the trim back to `Q`, and that trim's
    // second half, its fourth side shrunk to `Q`.
    let out_half = Side {
        curve: out.curve.clone(),
        range: (0.0, 0.5),
    };
    let back_half = Side {
        curve: back.curve.clone(),
        range: (0.5, 1.0),
    };
    let across = |normal: Vector, along: Vector, towards: Vector| -> OgeomResult<Vector> {
        let d = normal.cross(along).normalized(tol)?;
        Ok(if d.dot(towards) >= 0.0 { d } else { -d })
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
            let reach0 = from.distance(to) * full;
            let leave = across(out.normal, out_half.tangent(d, tol)?, to - from)?;
            let leave = turned_within(leave, out.normal, from, reach0, seam, tol)?;
            let arrive = -across(back.normal, back_half.tangent(1.0 - d, tol)?, from - to)?;
            let arrive = -turned_within(-arrive, back.normal, to, reach0, seam, tol)?;
            // Each arc lies within the hull of its Bezier points, the ends
            // and a third of the way along each tangent: held inside the
            // planes `seam` asks (the convex piece about the sharp edge,
            // which the solid holds), the arc stays inside them.
            let mut reach = reach0;
            for (normal, on) in seam {
                let lift = normal.dot(leave);
                if lift > 0.0 {
                    reach = reach.min(3.0 * (-normal.dot(from - *on)).max(0.0) / lift);
                }
                let lift = -normal.dot(arrive);
                if lift > 0.0 {
                    reach = reach.min(3.0 * (-normal.dot(to - *on)).max(0.0) / lift);
                }
            }
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
                (0, Across::Plane(out.normal)),
                (2, Across::Plane(back.normal)),
            ],
        };
        let network = ogeom_fit::fit_network(model, &patches, &curves, &index, &options, tol)?;
        // Its normal at `Q`, the limit along the shrunk side from one face's
        // to the other's, turns one way: the patch is a cone there that
        // does not fold back on itself.
        let surface = &network.surfaces[0];
        let limit = |w: f64| -> OgeomResult<Vector> {
            let (du, dw) = ogeom_geom::Surface::d1_at(surface, AT_Q, w, tol)?;
            // `dw` shrinks with the distance from `Q`; scaled back, the
            // cross product is the limit normal's own size.
            du.cross(dw * (1.0 / AT_Q)).normalized(tol)
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
        // Neither patch may stand past a plane of the vertex: the corner is
        // convex, so what stands past one lies outside the solid.
        let past = past_planes(model, &at_q, within, tol)?;
        if past > fit {
            ogeom_bail!(
                NotDone,
                "the patch at the sharp edge stands {past} past a face of the vertex"
            );
        }
        // The rest of the opening: the trim from `Q`'s second half, band
        // `r1`'s arc, the trim between the bands, band `r2`'s arc, the trim
        // back to `Q`'s first half, and the curve along the patch at `Q`,
        // filled by one patch G1 to every side.
        let side = |edge: &Shape, support: &Shape, tangent: bool| FillBoundary {
            edge: edge.clone(),
            support: Some(support.clone()),
            continuity: if tangent {
                Continuity::G1
            } else {
                Continuity::C0
            },
        };
        let boundary = vec![
            side(&out.edges[1], &out.filler, rest[0]),
            side(&arc1.0, &arc1.1, rest[1]),
            side(&middle.edges[0], &middle.filler, rest[2]),
            side(&arc2.0, &arc2.1, rest[3]),
            side(&back.edges[0], &back.filler, rest[4]),
            side(&network.edges[1], &at_q, true),
        ];
        // Held inside the region the caller names where the vertex is not
        // convex, and inside the vertex's planes where it is.
        let filled = match region {
            Some(region) => make_filling_n_region(model, &boundary, &[], fit, region, tol)?,
            None => make_filling_n_within(model, &boundary, &[], fit, within, tol)?,
        };
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
    let Some(both) = made else {
        ogeom_bail!(NotDone, "no fullness was tried");
    };
    Ok(both)
}

/// A trim's edges on `conic` from `v0` to `v1`: one, or two split at its
/// middle where the patch at a sharp edge's point ends.
pub(crate) fn trim_edges(
    model: &mut Model,
    conic: &Curve,
    (v0, v1): (&Shape, &Shape),
    split: bool,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let pieces: &[(f64, f64)] = if split {
        &[(0.0, 0.5), (0.5, 1.0)]
    } else {
        &[(0.0, 1.0)]
    };
    let mut ends = vec![v0.clone()];
    if split {
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
    Ok(edges)
}

/// What [`close_through_q`] holds its patches in.
#[derive(Clone, Copy)]
pub(crate) struct Held<'a> {
    /// Half-spaces, each an outward unit normal and a point on its plane,
    /// neither patch may stand past: a convex vertex's planes, past which
    /// lies nothing of the solid; none where the vertex is not convex.
    pub within: &'a [(Vector, Point)],
    /// Half-spaces the arcs of the patch at `Q` are kept inside: the convex
    /// piece about the sharp edge, where the vertex is not convex.
    pub seam: &'a [(Vector, Point)],
    /// The region the patch over the rest is held inside, where the vertex
    /// is not convex; the half-spaces `within` hold it otherwise.
    pub region: Option<Region<'a>>,
}

/// A face trim, oriented around the opening, and its supporting face.
pub(crate) struct Trim {
    pub curve: Curve,
    pub edges: Vec<Shape>,
    pub normal: Vector,
    pub filler: Shape,
}

/// How far `face` stands past the planes `within` (each an outward unit
/// normal and a point), read at the vertices of its mesh drawn at a quarter
/// of the bands' fit target.
fn past_planes(
    model: &Model,
    face: &Shape,
    within: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<f64> {
    // No plane to stand past (a vertex that is not convex): nothing to mesh.
    if within.is_empty() {
        return Ok(f64::NEG_INFINITY);
    }
    let deflection = ogeom_mesh::Deflection {
        chord: crate::marched::band_fit_target(tol) * 0.25,
        ..ogeom_mesh::Deflection::default()
    };
    let mesh = ogeom_mesh::triangulate_face(model, face, deflection, tol)?;
    Ok(mesh
        .positions
        .iter()
        .flat_map(|p| within.iter().map(move |(n, at)| n.dot(*p - *at)))
        .fold(f64::NEG_INFINITY, f64::max))
}

/// `dir`, a unit direction in the plane with unit normal `face` leaving
/// `from`, turned within that plane as little as keeps the point a third of
/// `reach` along it inside every plane of `seam` (an outward unit normal
/// and a point): its Bezier point stays in the convex piece they bound.
/// Where a quarter turn does not reach, `dir` is returned and the reach is
/// held instead.
fn turned_within(
    dir: Vector,
    face: Vector,
    from: Point,
    reach: f64,
    seam: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<Vector> {
    let mut d = dir;
    if reach <= tol.confusion() {
        return Ok(d);
    }
    for (normal, on) in seam {
        let room = 3.0 * (-normal.dot(from - *on)).max(0.0) / reach;
        let a = d.dot(*normal);
        if a <= room {
            continue;
        }
        let side = face.cross(d).normalized(tol)?;
        let side = if side.dot(*normal) < 0.0 { side } else { -side };
        let b = side.dot(*normal);
        let r = a.hypot(b);
        let turn = b.atan2(a) + (room / r).min(1.0).acos();
        if !(0.0..=core::f64::consts::FRAC_PI_2).contains(&turn) {
            return Ok(dir);
        }
        d = d * turn.cos() + side * turn.sin();
    }
    Ok(d)
}
