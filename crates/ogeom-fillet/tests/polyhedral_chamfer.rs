//! Chamfers of solids bounded by planes, built by setbacks in one pass:
//! every face shrunk to its setback corners, every bevelled edge a quad,
//! every vertex of three bevels a facet. Checked against closed forms that
//! do not share the construction, and against every order of the request.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, ChamferCorner};
use ogeom_math::Point;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A prism: the polygon `profile` in (y, z), wound counter-clockwise seen
/// from +x, swept from x = 0 to x = `length`.
fn prism(m: &mut Model, profile: &[(f64, f64)], length: f64) -> Shape {
    let n = profile.len();
    let mut points: Vec<Point> = profile
        .iter()
        .map(|&(y, z)| Point::new(0.0, y, z))
        .collect();
    points.extend(profile.iter().map(|&(y, z)| Point::new(length, y, z)));
    let mut rings: Vec<Vec<usize>> = vec![(0..n).rev().collect(), (n..2 * n).collect()];
    for i in 0..n {
        let j = (i + 1) % n;
        rings.push(vec![i, j, n + j, n + i]);
    }
    ogeom_algo::make_faceted_solid(m, &points, &rings, &[], T)
        .unwrap()
        .0
        .shape
}

/// A part of the corpus: its name, its profile and its length.
type Part = (&'static str, Vec<(f64, f64)>, f64);

/// The corpus: each part's name, its profile (rectilinear) and its length.
fn corpus() -> Vec<Part> {
    vec![
        (
            "box",
            vec![(0.0, 0.0), (40.0, 0.0), (40.0, 20.0), (0.0, 20.0)],
            30.0,
        ),
        (
            "L",
            vec![
                (0.0, 0.0),
                (40.0, 0.0),
                (40.0, 8.0),
                (8.0, 8.0),
                (8.0, 30.0),
                (0.0, 30.0),
            ],
            50.0,
        ),
        (
            "T",
            vec![
                (12.0, 0.0),
                (18.0, 0.0),
                (18.0, 24.0),
                (30.0, 24.0),
                (30.0, 30.0),
                (0.0, 30.0),
                (0.0, 24.0),
                (12.0, 24.0),
            ],
            40.0,
        ),
        (
            "U",
            vec![
                (0.0, 0.0),
                (40.0, 0.0),
                (40.0, 25.0),
                (34.0, 25.0),
                (34.0, 6.0),
                (6.0, 6.0),
                (6.0, 25.0),
                (0.0, 25.0),
            ],
            45.0,
        ),
        (
            "stepped block",
            vec![
                (0.0, 0.0),
                (50.0, 0.0),
                (50.0, 30.0),
                (35.0, 30.0),
                (35.0, 20.0),
                (20.0, 20.0),
                (20.0, 10.0),
                (0.0, 10.0),
            ],
            35.0,
        ),
        (
            "stepped bracket",
            vec![
                (0.0, 0.0),
                (60.0, 0.0),
                (60.0, 30.0),
                (45.0, 30.0),
                (45.0, 16.0),
                (30.0, 16.0),
                (30.0, 6.0),
                (0.0, 6.0),
            ],
            80.0,
        ),
    ]
}

/// The polygon's area, perimeter, and its convex and concave corners, for
/// a counter-clockwise ring.
fn profile_measures(profile: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    let n = profile.len();
    let (mut area, mut perimeter, mut convex, mut concave) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..n {
        let (a, b, c) = (profile[(i + n - 1) % n], profile[i], profile[(i + 1) % n]);
        area += b.0 * c.1 - c.0 * b.1;
        perimeter += (c.0 - b.0).hypot(c.1 - b.1);
        let turn = (b.0 - a.0) * (c.1 - b.1) - (b.1 - a.1) * (c.0 - b.0);
        if turn > 0.0 {
            convex += 1.0;
        } else {
            concave += 1.0;
        }
    }
    (area / 2.0, perimeter, convex, concave)
}

/// The volume of a rectilinear prism with every edge chamfered `d` and a
/// facet at every vertex, by sections across x. Farther than `d` from both
/// ends each section is the profile with every corner bevelled `d`: less
/// `d²/2` at a convex corner, more at a concave one. Within `d` of an end,
/// at `x` from it, the end's bevels have set the section back `d - x` all
/// round (area `A - P s + (convex - concave) s²` for a rectilinear profile
/// set back `s`), and the side edges' bevels with the corner facets bevel
/// its corners `x`.
fn chamfered_prism_volume(profile: &[(f64, f64)], length: f64, d: f64) -> f64 {
    let (a, p, convex, concave) = profile_measures(profile);
    let (cx, cc) = (convex, concave);
    let middle = a - cx * d * d / 2.0 + cc * d * d / 2.0;
    // The integral over 0..d of A - P (d - x) + (cx - cc)(d - x)² - (cx - cc) x²/2.
    let end = a * d - p * d * d / 2.0 + (cx - cc) * d * d * d / 3.0 - (cx - cc) * d * d * d / 6.0;
    middle * (length - 2.0 * d) + 2.0 * end
}

/// Whether two faces of `s` meet along an edge folded flat: one plane,
/// facing opposite ways, the one lying over the other. A volume, and a
/// check of each face alone, would not see it.
fn folded(m: &Model, s: &Shape) -> bool {
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let faces = ogeom_topo::ancestors_of(m, s, &e, ShapeType::Face).unwrap();
        if let [a, b] = faces.as_slice() {
            let (_, na) = ogeom_algo::face_normal(m, a, T).unwrap();
            let (_, nb) = ogeom_algo::face_normal(m, b, T).unwrap();
            if na.dot(nb) < -1.0 + 1e-9 {
                return true;
            }
        }
    }
    false
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// Every vertex's point, sorted, as bits; with the volume's bits.
fn key(m: &Model, s: &Shape) -> (u64, Vec<[u64; 3]>) {
    let mut out: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
        })
        .collect();
    out.sort_unstable();
    (volume(m, s).to_bits(), out)
}

fn all_edges(m: &Model, s: &Shape, d: f64) -> Vec<(Shape, Chamfer)> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .map(|e| (e, Chamfer::Symmetric(d)))
        .collect()
}

/// A permutation of `0..n` from a seed, without a random number crate.
fn shuffled(n: usize, seed: u64) -> Vec<usize> {
    let mut order: Vec<usize> = (0..n).collect();
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    for i in (1..n).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let j = usize::try_from(state % (i as u64 + 1)).unwrap();
        order.swap(i, j);
    }
    order
}

/// Every edge of every part of the corpus, at three distances: valid, a
/// face for every face, edge and vertex of the part, the volume of the
/// sections' closed form, and the same solid to the bit under five orders
/// of the request.
#[test]
fn every_edge_of_a_prismatic_part_chamfers_exactly_in_any_order() {
    for (name, profile, length) in corpus() {
        for d in [0.5, 1.0, 2.0] {
            let want = chamfered_prism_volume(&profile, length, d);
            let mut first: Option<(u64, Vec<[u64; 3]>)> = None;
            for seed in 0..5_u64 {
                let mut m = Model::new();
                let s = prism(&mut m, &profile, length);
                let (f, e, v) = (
                    explore_unique(&m, &s, ShapeType::Face).unwrap().len(),
                    explore_unique(&m, &s, ShapeType::Edge).unwrap().len(),
                    explore_unique(&m, &s, ShapeType::Vertex).unwrap().len(),
                );
                let specs = all_edges(&m, &s, d);
                let order = shuffled(specs.len(), seed);
                let specs: Vec<(Shape, Chamfer)> =
                    order.iter().map(|i| specs[*i].clone()).collect();
                let built = ogeom_fillet::chamfer_edges_cornered(
                    &mut m,
                    &s,
                    &specs,
                    ChamferCorner::Facet,
                    T,
                )
                .unwrap_or_else(|err| panic!("{name}, d {d}: {err}"));
                let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
                assert!(found.is_valid(), "{name}, d {d}: {:?}", found.problems);
                assert!(!folded(&m, &built.shape), "folded");
                let faces = explore_unique(&m, &built.shape, ShapeType::Face)
                    .unwrap()
                    .len();
                assert_eq!(faces, f + e + v, "{name}, d {d}: faces");
                let got = volume(&m, &built.shape);
                assert!(
                    (got - want).abs() <= 1e-9 * want,
                    "{name}, d {d}: volume {got} against {want}"
                );
                let k = key(&m, &built.shape);
                match &first {
                    None => first = Some(k),
                    Some(held) => assert!(*held == k, "{name}, d {d}: order {seed} differs"),
                }
            }
        }
    }
}

/// The stepped bracket, every edge at 2: every edge of the part has its
/// bevel, the side profiles' edges among them, and every vertex a
/// triangular facet, the mixed vertices at the step's inside corners too.
/// The concave bevel at the step ends on the side faces' facets, within
/// the part.
#[test]
fn the_stepped_bracket_bevels_every_edge_and_closes_every_vertex_with_a_facet() {
    let (_, profile, length) = corpus().pop().unwrap();
    let mut m = Model::new();
    let s = prism(&mut m, &profile, length);
    let edges = explore_unique(&m, &s, ShapeType::Edge).unwrap();
    let vertices = explore_unique(&m, &s, ShapeType::Vertex).unwrap();
    let specs = all_edges(&m, &s, 2.0);
    let built =
        ogeom_fillet::chamfer_edges_cornered(&mut m, &s, &specs, ChamferCorner::Facet, T).unwrap();
    for e in &edges {
        assert_eq!(
            built.history.generated(e).len(),
            1,
            "an edge without its bevel"
        );
    }
    for v in &vertices {
        let facet = built.history.generated(v);
        assert_eq!(facet.len(), 1, "a vertex without its facet");
        let corners = explore_unique(&m, &facet[0], ShapeType::Vertex).unwrap();
        assert_eq!(corners.len(), 3, "a facet that is not a triangle");
    }
    // Nothing stands outside the part's box: no bevel hangs past a side.
    for v in explore_unique(&m, &built.shape, ShapeType::Vertex).unwrap() {
        let p = m.node(&v).unwrap().data().as_vertex().unwrap().point;
        assert!(
            (-1e-9..=length + 1e-9).contains(&p.x)
                && (-1e-9..=60.0 + 1e-9).contains(&p.y)
                && (-1e-9..=30.0 + 1e-9).contains(&p.z),
            "a corner outside the part: {p:?}"
        );
    }
}

/// The concave edge at the step's inside corner alone: its bevel adds its
/// triangle's section all along and ends flush with both side faces.
#[test]
fn a_concave_bevel_alone_ends_flush_with_the_side_faces() {
    let (_, profile, length) = corpus().pop().unwrap();
    let mut m = Model::new();
    let s = prism(&mut m, &profile, length);
    let before = volume(&m, &s);
    let inside = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
            use ogeom_geom::Curve3d as _;
            let (p, q) = (c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap());
            (p.y - 45.0).abs() < 1e-9
                && (q.y - 45.0).abs() < 1e-9
                && (p.z - 16.0).abs() < 1e-9
                && (q.z - 16.0).abs() < 1e-9
        })
        .unwrap();
    let d = 2.0;
    let built = ogeom_fillet::chamfer_edges_cornered(
        &mut m,
        &s,
        &[(inside, Chamfer::Symmetric(d))],
        ChamferCorner::Facet,
        T,
    )
    .unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert!(!folded(&m, &built.shape), "folded");
    let got = volume(&m, &built.shape);
    let want = before + d * d / 2.0 * length;
    assert!((got - want).abs() <= 1e-9 * want, "{got} against {want}");
    for v in explore_unique(&m, &built.shape, ShapeType::Vertex).unwrap() {
        let p = m.node(&v).unwrap().data().as_vertex().unwrap().point;
        assert!((-1e-9..=length + 1e-9).contains(&p.x), "past a side: {p:?}");
    }
}

/// The side edges (along x) and the end outlines (in the end planes) of a
/// prism, by where their two vertices stand.
fn sides_and_ends(m: &Model, s: &Shape, length: f64) -> (Vec<Shape>, Vec<Shape>, Vec<Shape>) {
    let (mut sides, mut near, mut far) = (Vec::new(), Vec::new(), Vec::new());
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let (a, b) = ogeom_algo::edge_vertices(m, &e).unwrap().unwrap();
        let p = m.node(&a).unwrap().data().as_vertex().unwrap().point;
        let q = m.node(&b).unwrap().data().as_vertex().unwrap().point;
        if (p.x - q.x).abs() > 1e-9 {
            sides.push(e);
        } else if p.x.abs() < 1e-9 {
            near.push(e);
        } else {
            assert!((p.x - length).abs() < 1e-9);
            far.push(e);
        }
    }
    (sides, near, far)
}

/// Partial requests on every part of the corpus, against closed forms that
/// do not share the construction:
/// - one end's outline: within `d` of that end each section is the profile
///   set back `d - x` all round, mitred at every corner, reflex ones too
///   (`A - P s + (convex - concave) s²`), so the volume is
///   `A L - P d²/2 + (convex - concave) d³/3`;
/// - both ends' outlines: twice that taken off;
/// - the side edges only: every section the profile with its corners
///   bevelled, `L (A - convex d²/2 + concave d²/2)`.
///
/// Each in two orders of the request, the same solid to the bit.
#[test]
fn partial_requests_chamfer_to_the_sections_closed_forms() {
    let d = 1.5;
    for (name, profile, length) in corpus() {
        let (a, p, convex, concave) = profile_measures(&profile);
        let k = convex - concave;
        let one_end = a * length - p * d * d / 2.0 + k * d * d * d / 3.0;
        let cases = [
            ("one end", one_end),
            ("both ends", one_end - (a * length - one_end)),
            (
                "the sides",
                length * (a - convex * d * d / 2.0 + concave * d * d / 2.0),
            ),
        ];
        for (case, want) in cases {
            let mut first: Option<(u64, Vec<[u64; 3]>)> = None;
            for flip in [false, true] {
                let mut m = Model::new();
                let s = prism(&mut m, &profile, length);
                let (sides, near, far) = sides_and_ends(&m, &s, length);
                let mut edges = match case {
                    "one end" => near,
                    "both ends" => [near, far].concat(),
                    _ => sides,
                };
                if flip {
                    edges.reverse();
                }
                let specs: Vec<(Shape, Chamfer)> = edges
                    .into_iter()
                    .map(|e| (e, Chamfer::Symmetric(d)))
                    .collect();
                let built = ogeom_fillet::chamfer_edges_cornered(
                    &mut m,
                    &s,
                    &specs,
                    ChamferCorner::Facet,
                    T,
                )
                .unwrap_or_else(|err| panic!("{name}, {case}: {err}"));
                let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
                assert!(found.is_valid(), "{name}, {case}: {:?}", found.problems);
                assert!(!folded(&m, &built.shape), "folded");
                let got = volume(&m, &built.shape);
                assert!(
                    (got - want).abs() <= 1e-9 * want,
                    "{name}, {case}: volume {got} against {want}"
                );
                let k = key(&m, &built.shape);
                match &first {
                    None => first = Some(k),
                    Some(held) => assert!(*held == k, "{name}, {case}: the orders differ"),
                }
            }
        }
    }
}

/// On a box, whose faces have no reflex corner, the setback chamfer and
/// the wedge chamfer (built by booleans, mitred) agree wherever no vertex
/// holds three bevels, for many subsets of its edges.
#[test]
fn on_a_box_the_setback_chamfer_agrees_with_the_wedge_chamfer_below_three_bevels() {
    let (_, profile, length) = corpus().remove(0);
    for seed in 0..12_u64 {
        let mut m = Model::new();
        let s = prism(&mut m, &profile, length);
        let edges = explore_unique(&m, &s, ShapeType::Edge).unwrap();
        let mut held: ogeom_core::FastMap<u32, usize> = ogeom_core::FastMap::default();
        let mut taken = Vec::new();
        for i in shuffled(edges.len(), seed + 11) {
            let (a, b) = ogeom_algo::edge_vertices(&m, &edges[i]).unwrap().unwrap();
            let (ka, kb) = (a.node().index(), b.node().index());
            if held.get(&ka).copied().unwrap_or(0) < 2 && held.get(&kb).copied().unwrap_or(0) < 2 {
                *held.entry(ka).or_default() += 1;
                *held.entry(kb).or_default() += 1;
                taken.push((edges[i].clone(), Chamfer::Symmetric(1.5)));
            }
        }
        let setback =
            ogeom_fillet::chamfer_edges_cornered(&mut m, &s, &taken, ChamferCorner::Facet, T)
                .unwrap_or_else(|err| panic!("seed {seed}: {err}"));
        let found = ogeom_algo::check(&m, &setback.shape, T).unwrap();
        assert!(found.is_valid(), "seed {seed}: {:?}", found.problems);
        assert!(!folded(&m, &setback.shape), "folded");
        let wedged =
            ogeom_fillet::chamfer_edges_cornered(&mut m, &s, &taken, ChamferCorner::Point, T)
                .unwrap_or_else(|err| panic!("seed {seed}, wedges: {err}"));
        let (a, b) = (volume(&m, &setback.shape), volume(&m, &wedged.shape));
        assert!(
            (a - b).abs() <= 1e-6 * b,
            "seed {seed}: setbacks {a}, wedges {b}"
        );
    }
}

/// A distance past what the part holds (the bracket's plate is 6 thick, so
/// its two bevels meet at 3) is refused by name, not built: it lies
/// outside the engine's setback domain, and the facet path refuses it.
#[test]
fn a_chamfer_past_what_the_part_holds_is_refused_by_name() {
    let (_, profile, length) = corpus().pop().unwrap();
    let mut m = Model::new();
    let s = prism(&mut m, &profile, length);
    let specs = all_edges(&m, &s, 4.9);
    let refused = ogeom_fillet::chamfer_edges_cornered(&mut m, &s, &specs, ChamferCorner::Facet, T);
    let Err(why) = refused else {
        panic!("built");
    };
    assert!(
        matches!(why, ogeom_core::OgeomError::Construction(_)),
        "{why}"
    );
}

/// The L's end face, its two outline edges that meet at the face's reflex
/// corner (y = 8, z = 8), chamfered at `d`: the bevels mitre across the
/// corner, the end face's outline set back round it.
///
/// At `x` from that end each section loses the two strips the edges set
/// back by `s = d - x`: along the edge at z = 8 (32 long) the strip from
/// y = 8 - s, along the edge at y = 8 (22 long) the strip from z = 8 - s,
/// both running into the corner square `s x s` they share, so their union
/// is `54 s + s²` and the volume lost `27 d² + d³/3`. A wedge stopped at
/// each edge's end loses only `27 d²` and leaves the corner's `d³/3`
/// standing as a notch: a wrong solid, which the wedge chamfer built.
#[test]
fn bevels_meeting_at_a_face_s_reflex_corner_mitre_without_a_notch() {
    let (_, profile, length) = corpus().remove(1);
    let d = 1.5;
    let mut m = Model::new();
    let s = prism(&mut m, &profile, length);
    let before = volume(&m, &s);
    let (_, near, _) = sides_and_ends(&m, &s, length);
    let at_reflex: Vec<(Shape, Chamfer)> = near
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            [a, b].iter().any(|v| {
                let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
                (p.y - 8.0).abs() < 1e-9 && (p.z - 8.0).abs() < 1e-9
            })
        })
        .map(|e| (e, Chamfer::Symmetric(d)))
        .collect();
    assert_eq!(at_reflex.len(), 2, "the reflex corner's two outline edges");
    let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &at_reflex, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert!(!folded(&m, &built.shape), "folded");
    let got = volume(&m, &built.shape);
    let want = before - 27.0 * d * d - d * d * d / 3.0;
    assert!(
        (got - want).abs() <= 1e-9 * want,
        "volume {got} against {want} (a notch would leave {})",
        want + d * d * d / 3.0
    );
}

/// A box at `(x, y, z)` of `size`.
fn block(m: &mut Model, (x, y, z): (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
    let at = ogeom_math::Frame::new(
        Point::new(x, y, z),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    ogeom_algo::make_box(m, at, size, T).unwrap().shape
}

/// The stepped bracket as an application models it: a plate, a step and a
/// block fused, which leaves the end faces and the back face split along
/// seams where the boxes met, and the back's two upright edges each split
/// in three.
fn fused_bracket(m: &mut Model) -> Shape {
    let plate = block(m, (0.0, 0.0, 0.0), (80.0, 60.0, 6.0));
    let step = block(m, (0.0, 30.0, 6.0), (80.0, 30.0, 10.0));
    let top = block(m, (0.0, 45.0, 16.0), (80.0, 15.0, 14.0));
    let s = ogeom_bool::fuse(m, &plate, &step, T).unwrap().shape;
    ogeom_bool::fuse(m, &s, &top, T).unwrap().shape
}

/// The fused stepped bracket with every edge chamfered, the seams the
/// fuse left among the request, is the bracket built as one prism chamfered
/// alike: the same corners to the bit, under three orders of the request,
/// and the closed form's volume. A seam names nothing: it leaves no bevel
/// and no facet, and its faces come out one face, as the part's form has
/// them. Every edge of the form is bevelled, the split upright edges as
/// one bevel each, and every corner of the form closed by a triangle.
#[test]
fn the_fused_stepped_bracket_chamfers_as_the_prism_with_its_seams_naming_nothing() {
    let (_, profile, length) = corpus().pop().unwrap();
    for d in [1.0, 2.0] {
        let mut m = Model::new();
        let p = prism(&mut m, &profile, length);
        let prism_edges = explore_unique(&m, &p, ShapeType::Edge).unwrap();
        let from_prism = ogeom_fillet::chamfer_edges(&mut m, &p, &prism_edges, d, T).unwrap();
        let (want_faces, want) = (
            explore_unique(&m, &from_prism.shape, ShapeType::Face)
                .unwrap()
                .len(),
            key(&m, &from_prism.shape).1,
        );
        let closed = chamfered_prism_volume(&profile, length, d);
        for seed in 0..3_u64 {
            let mut m = Model::new();
            let s = fused_bracket(&mut m);
            let edges = explore_unique(&m, &s, ShapeType::Edge).unwrap();
            let vertices = explore_unique(&m, &s, ShapeType::Vertex).unwrap();
            let faces = explore_unique(&m, &s, ShapeType::Face).unwrap();
            assert_eq!(
                (faces.len(), edges.len(), vertices.len()),
                (16, 34, 20),
                "the fuse left no seams to test"
            );
            let specs = all_edges(&m, &s, d);
            let specs: Vec<(Shape, Chamfer)> = shuffled(specs.len(), seed)
                .iter()
                .map(|i| specs[*i].clone())
                .collect();
            let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T)
                .unwrap_or_else(|err| panic!("d {d}, order {seed}: {err}"));
            let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(found.is_valid(), "d {d}: {:?}", found.problems);
            assert!(!folded(&m, &built.shape), "folded");
            let got = volume(&m, &built.shape);
            assert!(
                (got - closed).abs() <= 1e-9 * closed,
                "d {d}: volume {got} against {closed}"
            );
            assert_eq!(
                explore_unique(&m, &built.shape, ShapeType::Face)
                    .unwrap()
                    .len(),
                want_faces,
                "d {d}: faces"
            );
            assert!(
                key(&m, &built.shape).1 == want,
                "d {d}, order {seed}: corners differ from the prism's"
            );
            // A seam is an edge whose two faces lie on one plane facing one
            // way: read here from the faces' normals, not the kernel's test.
            let seam = |e: &Shape| {
                let around = ogeom_topo::ancestors_of(&m, &s, e, ShapeType::Face).unwrap();
                let (_, a) = ogeom_algo::face_normal(&m, &around[0], T).unwrap();
                let (_, b) = ogeom_algo::face_normal(&m, &around[1], T).unwrap();
                a.dot(b) > 1.0 - 1e-12
            };
            let (seams, form): (Vec<&Shape>, Vec<&Shape>) = edges.iter().partition(|e| seam(e));
            assert_eq!(seams.len(), 6, "seams");
            let mut bevels: Vec<Shape> = Vec::new();
            for e in &form {
                let made = built.history.generated(e);
                assert_eq!(made.len(), 1, "an edge of the form without its bevel");
                if !bevels.iter().any(|b| b.is_same(&made[0])) {
                    bevels.push(made[0].clone());
                }
            }
            assert_eq!(bevels.len(), 24, "one bevel per edge of the form");
            for e in &seams {
                assert!(built.history.generated(e).is_empty(), "a seam bevelled");
                assert!(built.history.is_deleted(e), "a seam left standing");
            }
            let mut triangles = 0;
            for v in &vertices {
                match built.history.generated(v) {
                    [] => assert!(built.history.is_deleted(v), "a seam's vertex left"),
                    [facet] => {
                        let corners = explore_unique(&m, facet, ShapeType::Vertex).unwrap();
                        assert_eq!(corners.len(), 3, "a facet that is not a triangle");
                        triangles += 1;
                    }
                    more => panic!("a vertex with {} facets", more.len()),
                }
            }
            assert_eq!(triangles, 16, "a facet at every corner of the form");
            let mut kept: Vec<Shape> = Vec::new();
            for f in &faces {
                let into = built.history.modified(f);
                assert_eq!(into.len(), 1, "a face not carried to one face");
                if !kept.iter().any(|k| k.is_same(&into[0])) {
                    kept.push(into[0].clone());
                }
            }
            assert_eq!(kept.len(), 10, "the split faces come out one face each");
        }
    }
}

/// A request naming nothing but seams leaves the solid as it is: alone,
/// in a set, and in each corner's spelling.
#[test]
fn a_chamfer_naming_only_seams_leaves_the_solid_as_it_is() {
    let mut m = Model::new();
    let s = fused_bracket(&mut m);
    let seams: Vec<Shape> = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let around = ogeom_topo::ancestors_of(&m, &s, e, ShapeType::Face).unwrap();
            let (_, a) = ogeom_algo::face_normal(&m, &around[0], T).unwrap();
            let (_, b) = ogeom_algo::face_normal(&m, &around[1], T).unwrap();
            a.dot(b) > 1.0 - 1e-12
        })
        .collect();
    assert_eq!(seams.len(), 6);
    // The solid as it stood: its corners and its volume to the bit, its
    // faces as many (a solid shared with the fuse's inputs comes back a
    // copy).
    let as_given = (
        key(&m, &s),
        explore_unique(&m, &s, ShapeType::Face).unwrap().len(),
    );
    let same = |m: &Model, built: &Shape| {
        (
            key(m, built),
            explore_unique(m, built, ShapeType::Face).unwrap().len(),
        ) == as_given
    };
    let alone = ogeom_fillet::chamfer_edge(&mut m, &s, &seams[0], 1.0, T).unwrap();
    assert!(same(&m, &alone.shape), "a seam alone changed the solid");
    let specs: Vec<(Shape, Chamfer)> = seams
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(1.0)))
        .collect();
    for corner in [ChamferCorner::Facet, ChamferCorner::Point] {
        let built = ogeom_fillet::chamfer_edges_cornered(&mut m, &s, &specs, corner, T).unwrap();
        assert!(
            same(&m, &built.shape),
            "{corner:?}: seams changed the solid"
        );
    }
}

/// Two plates fused where they touch along one line, the upper one's
/// underside on the lower one's top plane: across that line the plane
/// holds material below on one side and above on the other. The faces are
/// coplanar but the edge is real, so a chamfer naming it is judged as one
/// (refused: four faces meet there), never passed over as a seam.
#[test]
fn coplanar_faces_with_the_material_on_opposite_sides_meet_at_a_real_edge() {
    let mut m = Model::new();
    let lower = block(&mut m, (0.0, 0.0, 0.0), (80.0, 30.0, 6.0));
    let upper = block(&mut m, (0.0, 30.0, 6.0), (80.0, 30.0, 6.0));
    let s = ogeom_bool::fuse(&mut m, &lower, &upper, T).unwrap().shape;
    let touching = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            [a, b].iter().all(|v| {
                let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
                (p.y - 30.0).abs() < 1e-9 && (p.z - 6.0).abs() < 1e-9
            })
        })
        .unwrap();
    let around = ogeom_topo::ancestors_of(&m, &s, &touching, ShapeType::Face).unwrap();
    let flat: Vec<&Shape> = around
        .iter()
        .filter(|f| {
            let (_, n) = ogeom_algo::face_normal(&m, f, T).unwrap();
            n.z.abs() > 1.0 - 1e-12
        })
        .collect();
    let [top, underside] = flat.as_slice() else {
        panic!("{} faces on the plane at the line", flat.len());
    };
    assert!(!ogeom_heal::same_domain(&m, top, underside, T).unwrap());
    let err = ogeom_fillet::chamfer_edge(&mut m, &s, &touching, 1.0, T).unwrap_err();
    assert!(err.to_string().contains("exactly two faces"), "{err}");
}

/// The volume of the convex polytope `{x : n·x <= c}` over `planes`, each
/// an outward normal and its offset, with the count of its faces: every
/// point where three planes meet inside all the others is a corner; each
/// face's corners, in turn about its centroid, give its area; the volume
/// is the sum over faces of a third of area times offset. Shares nothing
/// with the engine.
fn convex_polytope(planes: &[(ogeom_math::Vector, f64)]) -> (f64, usize) {
    use ogeom_math::Vector;
    let mut corners: Vec<Vector> = Vec::new();
    for i in 0..planes.len() {
        for j in i + 1..planes.len() {
            for k in j + 1..planes.len() {
                let (a, b, c) = (planes[i], planes[j], planes[k]);
                let det = a.0.dot(b.0.cross(c.0));
                if det.abs() < 1e-12 {
                    continue;
                }
                let x = (b.0.cross(c.0) * a.1 + c.0.cross(a.0) * b.1 + a.0.cross(b.0) * c.1) / det;
                if planes.iter().all(|(n, o)| n.dot(x) <= o + 1e-9)
                    && !corners.iter().any(|y| (*y - x).magnitude() < 1e-9)
                {
                    corners.push(x);
                }
            }
        }
    }
    let (mut volume, mut faces) = (0.0, 0);
    for (n, o) in planes {
        let on: Vec<Vector> = corners
            .iter()
            .copied()
            .filter(|x| (n.dot(*x) - o).abs() < 1e-9)
            .collect();
        if on.len() < 3 {
            continue;
        }
        let (sum, count) = on
            .iter()
            .fold((Vector::ZERO, 0.0), |(s, c), x| (s + *x, c + 1.0));
        let centre = sum / count;
        let u = (on[0] - centre) / (on[0] - centre).magnitude();
        let w = n.cross(u);
        let mut turn: Vec<(f64, Vector)> = on
            .iter()
            .map(|x| ((*x - centre).dot(w).atan2((*x - centre).dot(u)), *x))
            .collect();
        turn.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut area = 0.0;
        for k in 0..turn.len() {
            let (p, q) = (turn[k].1 - centre, turn[(k + 1) % turn.len()].1 - centre);
            area += p.cross(q).dot(*n) / 2.0;
        }
        if area > 1e-9 {
            faces += 1;
            volume += area * o / 3.0;
        }
    }
    (volume, faces)
}

/// A plate 20 by 20 by 2, its front edges top and bottom chamfered 1.9
/// and every edge of its back, its sides' top and bottom and its back
/// uprights 0.5; the front uprights stay sharp. The front face's setback
/// lines cross: it goes, and its two bevels meet along a ridge at half
/// height, the side faces notched up one bevel and down the other. Each
/// back corner of three bevels is a facet, as a convex plate cut by every
/// one of its planes is: the result is the polytope of the plate's, the
/// bevels' and the facets' planes, its volume and its face count, under
/// three orders of the request. The front face and its sharp uprights are
/// gone from the history. Before, the band sent the whole request down the
/// wedge path, its back corners mitred.
#[test]
fn a_thin_wall_chamfered_on_both_edges_meets_in_a_ridge() {
    use ogeom_math::Vector;
    let (front, rest) = (1.9, 0.5);
    let mut planes: Vec<(Vector, f64)> = vec![
        (Vector::new(-1.0, 0.0, 0.0), 0.0),
        (Vector::new(1.0, 0.0, 0.0), 20.0),
        (Vector::new(0.0, -1.0, 0.0), 0.0),
        (Vector::new(0.0, 1.0, 0.0), 20.0),
        (Vector::new(0.0, 0.0, -1.0), 0.0),
        (Vector::new(0.0, 0.0, 1.0), 2.0),
    ];
    // A bevel of `d` on both faces along the edge where the outward
    // normals `a` and `b` meet, `at` a point of the edge: through the line
    // `d` in along the face of normal `a`, away from the other face.
    let bevel = |a: Vector, b: Vector, at: Vector, d: f64| {
        let n = (a + b) / (a + b).magnitude();
        (n, n.dot(at - b * d))
    };
    let (x, y, z) = (Vector::X, Vector::Y, Vector::Z);
    planes.push(bevel(-y, -z, Vector::ZERO, front));
    planes.push(bevel(-y, z, 2.0 * z, front));
    for (a, b, at) in [
        (y, -z, 20.0 * y),
        (y, z, 20.0 * y + 2.0 * z),
        (-x, -z, Vector::ZERO),
        (-x, z, 2.0 * z),
        (x, -z, 20.0 * x),
        (x, z, 20.0 * x + 2.0 * z),
        (-x, y, 20.0 * y),
        (x, y, 20.0 * x + 20.0 * y),
    ] {
        planes.push(bevel(a, b, at, rest));
    }
    // Each back corner's facet, through its three setback corners: on each
    // face, `rest` in from both bevelled edges.
    for (corner, sx, sz) in [
        (20.0 * y, -1.0, -1.0),
        (20.0 * y + 2.0 * z, -1.0, 1.0),
        (20.0 * x + 20.0 * y, 1.0, -1.0),
        (20.0 * x + 20.0 * y + 2.0 * z, 1.0, 1.0),
    ] {
        let (ix, iy, iz) = (-sx * rest * x, -rest * y, -sz * rest * z);
        let (p, q, r) = (corner + iy + iz, corner + ix + iz, corner + ix + iy);
        let mut n = (q - p).cross(r - p);
        n = n / n.magnitude();
        if n.dot(Vector::new(sx, 1.0, sz)) < 0.0 {
            n = -n;
        }
        planes.push((n, n.dot(p)));
    }
    let (want, want_faces) = convex_polytope(&planes);
    let mut first: Option<Vec<[u64; 3]>> = None;
    for seed in 0..3_u64 {
        let mut m = Model::new();
        let s = block(&mut m, (0.0, 0.0, 0.0), (20.0, 20.0, 2.0));
        let front_face = explore_unique(&m, &s, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let (_, n) = ogeom_algo::face_normal(&m, f, T).unwrap();
                n.y < -0.5
            })
            .unwrap();
        let mut sharp: Vec<Shape> = Vec::new();
        let specs: Vec<(Shape, Chamfer)> = explore_unique(&m, &s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter_map(|e| {
                let (a, b) = ogeom_algo::edge_vertices(&m, &e).unwrap().unwrap();
                let point = |v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;
                let (p, q) = (point(&a), point(&b));
                let on_front = p.y.max(q.y) < 1e-9;
                let upright = (p.z - q.z).abs() > 1.0;
                match (on_front, upright) {
                    (true, true) => {
                        sharp.push(e);
                        None
                    }
                    (true, false) => Some((e, Chamfer::Symmetric(front))),
                    _ => Some((e, Chamfer::Symmetric(rest))),
                }
            })
            .collect();
        assert_eq!(specs.len(), 10);
        let specs: Vec<(Shape, Chamfer)> = shuffled(specs.len(), seed)
            .iter()
            .map(|i| specs[*i].clone())
            .collect();
        let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        assert!(!folded(&m, &built.shape), "folded");
        let got = volume(&m, &built.shape);
        assert!(
            (got - want).abs() <= 1e-9 * want,
            "volume {got} against {want}"
        );
        let faces = explore_unique(&m, &built.shape, ShapeType::Face)
            .unwrap()
            .len();
        assert_eq!(faces, want_faces, "faces");
        assert!(built.history.is_deleted(&front_face), "the front face left");
        for e in &sharp {
            assert!(built.history.is_deleted(e), "a front upright left");
        }
        let k = key(&m, &built.shape).1;
        match &first {
            None => first = Some(k),
            Some(held) => assert!(*held == k, "order {seed} differs"),
        }
    }
}

/// A plate 60 by 40 by 6 with a 20 by 10 hole through it, cut by a
/// boolean: its top and bottom faces each hold the hole as an inner ring.
fn holed_plate(m: &mut Model) -> Shape {
    let plate = block(m, (0.0, 0.0, 0.0), (60.0, 40.0, 6.0));
    let hole = block(m, (20.0, 15.0, -1.0), (20.0, 10.0, 8.0));
    ogeom_bool::cut(m, &plate, &hole, T).unwrap().shape
}

/// [`chamfered_prism_volume`] for a profile with holes: the outer ring
/// counter-clockwise and each hole clockwise, so their areas and corners
/// sum as one region's (a hole's corners are the region's concave ones).
fn chamfered_holed_prism_volume(loops: &[Vec<(f64, f64)>], length: f64, d: f64) -> f64 {
    let (mut a, mut p, mut cx, mut cc) = (0.0, 0.0, 0.0, 0.0);
    for ring in loops {
        let (ra, rp, rx, rc) = profile_measures(ring);
        (a, p, cx, cc) = (a + ra, p + rp, cx + rx, cc + rc);
    }
    let middle = a - cx * d * d / 2.0 + cc * d * d / 2.0;
    let end = a * d - p * d * d / 2.0 + (cx - cc) * d * d * d / 3.0 - (cx - cc) * d * d * d / 6.0;
    middle * (length - 2.0 * d) + 2.0 * end
}

/// The holed plate with every edge chamfered, the hole's too: valid,
/// nothing folded, the volume of the holed profile's sections in closed
/// form, a face for every face, edge and vertex of the plate, and the same
/// corners to the bit under three orders. The top and bottom faces come
/// out holed, the hole's rim set back outward.
#[test]
fn a_plate_with_a_hole_chamfers_every_edge_exactly_in_any_order() {
    // The profile across z: the plate's outline and the hole's, the hole
    // walked the other way.
    let loops = vec![
        vec![(0.0, 0.0), (60.0, 0.0), (60.0, 40.0), (0.0, 40.0)],
        vec![(20.0, 15.0), (20.0, 25.0), (40.0, 25.0), (40.0, 15.0)],
    ];
    for d in [0.5, 1.0, 2.0] {
        let want = chamfered_holed_prism_volume(&loops, 6.0, d);
        let mut first: Option<Vec<[u64; 3]>> = None;
        for seed in 0..3_u64 {
            let mut m = Model::new();
            let s = holed_plate(&mut m);
            let (f, e, v) = (
                explore_unique(&m, &s, ShapeType::Face).unwrap().len(),
                explore_unique(&m, &s, ShapeType::Edge).unwrap().len(),
                explore_unique(&m, &s, ShapeType::Vertex).unwrap().len(),
            );
            assert_eq!(
                (f, e, v),
                (10, 24, 16),
                "the cut left no holed faces to test"
            );
            let specs = all_edges(&m, &s, d);
            let specs: Vec<(Shape, Chamfer)> = shuffled(specs.len(), seed)
                .iter()
                .map(|i| specs[*i].clone())
                .collect();
            let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T)
                .unwrap_or_else(|err| panic!("d {d}: {err}"));
            let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(found.is_valid(), "d {d}: {:?}", found.problems);
            assert!(!folded(&m, &built.shape), "folded");
            let got = volume(&m, &built.shape);
            assert!(
                (got - want).abs() <= 1e-9 * want,
                "d {d}: volume {got} against {want}"
            );
            let faces = explore_unique(&m, &built.shape, ShapeType::Face).unwrap();
            assert_eq!(faces.len(), f + e + v, "d {d}: faces");
            let holed = faces
                .iter()
                .filter(|face| m.children_of(face).unwrap().len() == 2)
                .count();
            assert_eq!(holed, 2, "d {d}: the top and bottom keep their hole");
            let k = key(&m, &built.shape).1;
            match &first {
                None => first = Some(k),
                Some(held) => assert!(*held == k, "d {d}: order {seed} differs"),
            }
        }
    }
}

/// The hole's top rim alone, bevelled `d`: the four bevels mitre at the
/// top face's reflex corners round the hole, so at depth `z` the hole is
/// widened by `d - z` all round, and the plate loses
/// `d²(a + b) + 4d³/3` for a hole `a` by `b`.
#[test]
fn a_hole_s_rim_alone_bevels_to_its_closed_form() {
    let d = 1.5;
    let mut m = Model::new();
    let s = holed_plate(&mut m);
    let rim: Vec<(Shape, Chamfer)> = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            [a, b].iter().all(|v| {
                let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
                (p.z - 6.0).abs() < 1e-9
                    && (20.0..=40.0).contains(&p.x)
                    && (15.0..=25.0).contains(&p.y)
            })
        })
        .map(|e| (e, Chamfer::Symmetric(d)))
        .collect();
    assert_eq!(rim.len(), 4);
    let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &rim, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert!(!folded(&m, &built.shape), "folded");
    let lost = 60.0 * 40.0 * 6.0 - 20.0 * 10.0 * 6.0 - volume(&m, &built.shape);
    let want = d * d * (20.0 + 10.0) + 4.0 * d * d * d / 3.0;
    assert!(
        (lost - want).abs() <= 1e-9 * want,
        "lost {lost} against {want}"
    );
    for (e, _) in &rim {
        assert_eq!(
            built.history.generated(e).len(),
            1,
            "a rim edge without its bevel"
        );
    }
}

/// The plane of a bevel along a convex edge through `at` with direction
/// `u`, between faces of outward normals `na` and `nb`, set back `da` on
/// the first face and `db` on the second, measured in each face square to
/// the edge: an outward normal and its offset. Independent of the engine:
/// each contact line is the edge moved into its face, away from the other.
fn bevel_plane(
    at: ogeom_math::Vector,
    u: ogeom_math::Vector,
    (na, da): (ogeom_math::Vector, f64),
    (nb, db): (ogeom_math::Vector, f64),
) -> (ogeom_math::Vector, f64) {
    let inward = |n: ogeom_math::Vector, other: ogeom_math::Vector| {
        let m = n.cross(u);
        let m = m / m.magnitude();
        if m.dot(other) > 0.0 { -m } else { m }
    };
    let (ca, cb) = (at + inward(na, nb) * da, at + inward(nb, na) * db);
    let mut n = u.cross(cb - ca);
    n = n / n.magnitude();
    if n.dot(na + nb) < 0.0 {
        n = -n;
    }
    (n, n.dot(ca))
}

/// A box's six planes, outward, for a box from the origin of `size`.
fn box_planes((x, y, z): (f64, f64, f64)) -> Vec<(ogeom_math::Vector, f64)> {
    use ogeom_math::Vector;
    vec![
        (-Vector::X, 0.0),
        (Vector::X, x),
        (-Vector::Y, 0.0),
        (Vector::Y, y),
        (-Vector::Z, 0.0),
        (Vector::Z, z),
    ]
}

/// Two bottom edges of a box meeting at a corner beside its sharp upright,
/// bevelled unequally: the bevels mitre along the line their planes share,
/// the larger running on over the corner to the side face, with no facet
/// between them, as `chamfer_edges_with` states. The result is the box cut
/// by both bevels' planes: its volume and its face count, for symmetric
/// and per-face spellings, either distance the larger.
#[test]
fn two_unequal_bevels_at_a_convex_corner_mitre_without_a_cap() {
    use ogeom_math::Vector;
    let size = (20.0, 20.0, 2.0);
    for (front, side) in [(1.9, 0.5), (0.5, 1.9), (1.2, 0.3), (1.0, 1.0)] {
        let mut planes = box_planes(size);
        planes.push(bevel_plane(
            Vector::ZERO,
            Vector::X,
            (-Vector::Y, front),
            (-Vector::Z, front),
        ));
        planes.push(bevel_plane(
            Vector::ZERO,
            Vector::Y,
            (-Vector::X, side),
            (-Vector::Z, side),
        ));
        let (want, want_faces) = convex_polytope(&planes);
        let mut m = Model::new();
        let s = block(&mut m, (0.0, 0.0, 0.0), size);
        let point = |m: &Model, v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;
        let specs: Vec<(Shape, Chamfer)> = explore_unique(&m, &s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter_map(|e| {
                let (a, b) = ogeom_algo::edge_vertices(&m, &e).unwrap().unwrap();
                let (p, q) = (point(&m, &a), point(&m, &b));
                let bottom = p.z.max(q.z) < 1e-9;
                if bottom && p.y.max(q.y) < 1e-9 {
                    Some((e, Chamfer::Symmetric(front)))
                } else if bottom && p.x.max(q.x) < 1e-9 {
                    Some((e, Chamfer::Symmetric(side)))
                } else {
                    None
                }
            })
            .collect();
        assert_eq!(specs.len(), 2);
        let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{front}, {side}: {:?}", found.problems);
        let got = volume(&m, &built.shape);
        assert!(
            (got - want).abs() <= 1e-9 * want,
            "{front}, {side}: volume {got} against {want}"
        );
        let faces = explore_unique(&m, &built.shape, ShapeType::Face)
            .unwrap()
            .len();
        assert_eq!(faces, want_faces, "{front}, {side}: faces");
    }
}

/// A plate 20 by 20 by 2 with every edge bevelled past half its thickness:
/// each side face's setback lines cross, so all four go, their top and
/// bottom bevels meeting in ridges that end on the upright bevels, and
/// each corner of three bevels keeps the part of its facet inside them.
/// The plate is convex and every bevel convex, so the result is the plate
/// cut by the planes of its faces, its bevels and its corner facets (each
/// facet through the setback corners `d` in from both edges on each face),
/// its volume and face count, the same corners to the bit in any order.
#[test]
fn a_plate_bevelled_all_round_past_half_its_thickness_keeps_the_polytope() {
    use ogeom_math::Vector;
    let size = (20.0, 20.0, 2.0);
    for d in [1.2, 1.9] {
        let mut planes = box_planes(size);
        let (x, y, z) = (Vector::X, Vector::Y, Vector::Z);
        let (sx, sy, sz) = (size.0 * x, size.1 * y, size.2 * z);
        // The twelve edges: a point, a direction and the two faces.
        let edges = [
            (Vector::ZERO, x, -y, -z),
            (sy, x, y, -z),
            (sz, x, -y, z),
            (sy + sz, x, y, z),
            (Vector::ZERO, y, -x, -z),
            (sx, y, x, -z),
            (sz, y, -x, z),
            (sx + sz, y, x, z),
            (Vector::ZERO, z, -x, -y),
            (sx, z, x, -y),
            (sy, z, -x, y),
            (sx + sy, z, x, y),
        ];
        for (at, u, na, nb) in edges {
            planes.push(bevel_plane(at, u, (na, d), (nb, d)));
        }
        for (cx, cy, cz) in [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (1.0, 1.0, 0.0),
        ]
        .into_iter()
        .flat_map(|(a, b, _)| [(a, b, 0.0), (a, b, 1.0)])
        {
            let corner = Vector::new(cx * size.0, cy * size.1, cz * size.2);
            let inward = Vector::new(1.0 - 2.0 * cx, 1.0 - 2.0 * cy, 1.0 - 2.0 * cz) * d;
            let (ix, iy, iz) = (inward.x * x, inward.y * y, inward.z * z);
            let (p, q, r) = (corner + iy + iz, corner + ix + iz, corner + ix + iy);
            let mut n = (q - p).cross(r - p);
            n = n / n.magnitude();
            if n.dot(-inward) < 0.0 {
                n = -n;
            }
            planes.push((n, n.dot(p)));
        }
        let (want, want_faces) = convex_polytope(&planes);
        let mut first: Option<Vec<[u64; 3]>> = None;
        for seed in 0..3_u64 {
            let mut m = Model::new();
            let s = block(&mut m, (0.0, 0.0, 0.0), size);
            let specs = all_edges(&m, &s, d);
            let specs: Vec<(Shape, Chamfer)> = shuffled(specs.len(), seed)
                .iter()
                .map(|i| specs[*i].clone())
                .collect();
            let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T)
                .unwrap_or_else(|err| panic!("d {d}: {err}"));
            let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(found.is_valid(), "d {d}: {:?}", found.problems);
            assert!(!folded(&m, &built.shape), "folded");
            let got = volume(&m, &built.shape);
            assert!(
                (got - want).abs() <= 1e-9 * want,
                "d {d}: volume {got} against {want}"
            );
            let faces = explore_unique(&m, &built.shape, ShapeType::Face)
                .unwrap()
                .len();
            assert_eq!(faces, want_faces, "d {d}: faces");
            let k = key(&m, &built.shape).1;
            match &first {
                None => first = Some(k),
                Some(held) => assert!(*held == k, "d {d}: order {seed} differs"),
            }
        }
    }
}

/// A plate tapered towards its top, its front a trapezoid between two
/// slanted edges. Both slanted front edges bevelled 4: their setback lines
/// cross on the front face below its top, so the front keeps the triangle
/// below the crossing and the bevels meet in a ridge above it up to the
/// top. The tapered plate is convex, so the result is it cut by the two
/// bevels' planes: its volume and face count.
#[test]
fn bevels_crossing_on_a_face_at_an_angle_leave_it_the_part_below() {
    use ogeom_math::Vector;
    let mut points: Vec<Point> = [
        (0.0, 0.0, 0.0),
        (20.0, 0.0, 0.0),
        (20.0, 10.0, 0.0),
        (0.0, 10.0, 0.0),
        (5.0, 0.0, 4.0),
        (15.0, 0.0, 4.0),
        (15.0, 10.0, 4.0),
        (5.0, 10.0, 4.0),
    ]
    .iter()
    .map(|&(x, y, z)| Point::new(x, y, z))
    .collect();
    let rings: Vec<Vec<usize>> = vec![
        vec![3, 2, 1, 0],
        vec![4, 5, 6, 7],
        vec![0, 1, 5, 4],
        vec![2, 3, 7, 6],
        vec![3, 0, 4, 7],
        vec![1, 2, 6, 5],
    ];
    let plane_of = |ring: &[usize], points: &[Point]| {
        let (a, b, c) = (points[ring[0]], points[ring[1]], points[ring[2]]);
        let n = (b - a).cross(c - b);
        let n = n / n.magnitude();
        (n, n.dot(a.to_vector()))
    };
    let mut planes: Vec<(Vector, f64)> = rings.iter().map(|r| plane_of(r, &points)).collect();
    let d = 4.0;
    let (front, left, right) = (planes[2].0, planes[4].0, planes[5].0);
    let along = |a: usize, b: usize, points: &[Point]| {
        let u = points[b] - points[a];
        u / u.magnitude()
    };
    planes.push(bevel_plane(
        points[0].to_vector(),
        along(0, 4, &points),
        (front, d),
        (left, d),
    ));
    planes.push(bevel_plane(
        points[1].to_vector(),
        along(1, 5, &points),
        (front, d),
        (right, d),
    ));
    let (want, want_faces) = convex_polytope(&planes);
    let mut m = Model::new();
    let s = ogeom_algo::make_faceted_solid(&mut m, &points, &rings, &[], T)
        .unwrap()
        .0
        .shape;
    points.clear();
    let specs: Vec<(Shape, Chamfer)> = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            let pa = m.node(&a).unwrap().data().as_vertex().unwrap().point;
            let pb = m.node(&b).unwrap().data().as_vertex().unwrap().point;
            pa.y.max(pb.y) < 1e-9 && (pa.z - pb.z).abs() > 1.0
        })
        .map(|e| (e, Chamfer::Symmetric(d)))
        .collect();
    assert_eq!(specs.len(), 2);
    let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert!(!folded(&m, &built.shape), "folded");
    let got = volume(&m, &built.shape);
    assert!(
        (got - want).abs() <= 1e-9 * want,
        "volume {got} against {want}"
    );
    let faces = explore_unique(&m, &built.shape, ShapeType::Face)
        .unwrap()
        .len();
    assert_eq!(faces, want_faces, "faces");
}

/// A slab: the polygon `outline` in (x, y), wound counter-clockwise seen
/// from +z, from z = 0 to z = `height`.
fn slab(m: &mut Model, outline: &[(f64, f64)], height: f64) -> Shape {
    let n = outline.len();
    let mut points: Vec<Point> = outline
        .iter()
        .map(|&(x, y)| Point::new(x, y, 0.0))
        .collect();
    points.extend(outline.iter().map(|&(x, y)| Point::new(x, y, height)));
    let mut rings: Vec<Vec<usize>> = vec![(0..n).rev().collect(), (n..2 * n).collect()];
    for i in 0..n {
        let j = (i + 1) % n;
        rings.push(vec![i, j, n + j, n + i]);
    }
    ogeom_algo::make_faceted_solid(m, &points, &rings, &[], T)
        .unwrap()
        .0
        .shape
}

/// A web's section once both its top edges are bevelled past half its
/// width: `(u, z)` with `u` across the web from its outer face, the web
/// `w` wide and `t` thick, each bevel `d` on both its faces, meeting in a
/// ridge at mid-web. Its area and its first moment about `u = 0`.
fn ridged_web(w: f64, t: f64, d: f64) -> (f64, f64) {
    assert!(
        2.0 * d > w && d < t,
        "the bevels meet above the web's sides"
    );
    let ring = [
        (0.0, 0.0),
        (w, 0.0),
        (w, t - d),
        (w / 2.0, t - d + w / 2.0),
        (0.0, t - d),
    ];
    let (mut area, mut moment) = (0.0, 0.0);
    for k in 0..ring.len() {
        let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
        let cross = a.0 * b.1 - b.0 * a.1;
        area += cross / 2.0;
        moment += (a.0 + b.0) * cross / 6.0;
    }
    (area, moment)
}

/// The top edges of `s` at height `top`, each bevelled `d`.
fn top_edges(m: &Model, s: &Shape, top: f64, d: f64) -> Vec<(Shape, Chamfer)> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
            [a, b].iter().all(|v| {
                let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
                (p.z - top).abs() < 1e-9
            })
        })
        .map(|e| (e, Chamfer::Symmetric(d)))
        .collect()
}

/// A frame 24 by 16 by 4, its 18 by 10 hole leaving webs 3 wide all
/// round, its top edges outside and in the hole bevelled 2: past half the
/// web, so the bevels consume the holed top face and meet in a ridge. At
/// each corner the outer corner and the hole's stand on one diagonal, so
/// the frame is four webs mitred there, each the ridged section swept
/// along its length: a web `L` long outside is `L - 2u` long at `u` in,
/// and holds `A L - 2 ∫u dA`.
#[test]
fn a_frame_bevelled_past_half_its_webs_meets_in_a_mitred_ridge() {
    let (w, t, d) = (3.0, 4.0, 2.0);
    let (area, moment) = ridged_web(w, t, d);
    let want = 2.0 * (area * 24.0 - 2.0 * moment) + 2.0 * (area * 16.0 - 2.0 * moment);
    let mut first: Option<Vec<[u64; 3]>> = None;
    for seed in 0..3_u64 {
        let mut m = Model::new();
        let outer = block(&mut m, (0.0, 0.0, 0.0), (24.0, 16.0, t));
        let hole = block(&mut m, (w, w, -1.0), (18.0, 10.0, t + 2.0));
        let s = ogeom_bool::cut(&mut m, &outer, &hole, T).unwrap().shape;
        let top = explore_unique(&m, &s, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| ogeom_algo::face_normal(&m, f, T).unwrap().1.z > 0.5)
            .unwrap();
        assert_eq!(
            m.children_of(&top).unwrap().len(),
            2,
            "the top holds the hole"
        );
        let specs = top_edges(&m, &s, t, d);
        assert_eq!(specs.len(), 8);
        let specs: Vec<(Shape, Chamfer)> = shuffled(specs.len(), seed)
            .iter()
            .map(|i| specs[*i].clone())
            .collect();
        let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        assert!(!folded(&m, &built.shape), "folded");
        let got = volume(&m, &built.shape);
        assert!(
            (got - want).abs() <= 1e-9 * want,
            "volume {got} against {want}"
        );
        // The bottom, the eight walls and the eight bevels.
        let faces = explore_unique(&m, &built.shape, ShapeType::Face)
            .unwrap()
            .len();
        assert_eq!(faces, 17, "faces");
        assert!(built.history.is_deleted(&top), "the top face left");
        for (e, _) in &specs {
            assert_eq!(
                built.history.generated(e).len(),
                1,
                "an edge without its bevel"
            );
        }
        let k = key(&m, &built.shape).1;
        match &first {
            None => first = Some(k),
            Some(held) => assert!(*held == k, "order {seed} differs"),
        }
    }
}

/// An L-shaped slab 4 thick, legs 20 long and 3 wide, its long top edges
/// bevelled 2 on both sides of each leg, its two end edges left sharp: the
/// bevels consume the top face, reflex corner and all, and meet in a
/// ridge that turns on the diagonal through the outline's corner and its
/// reflex corner. Each leg is the ridged section swept from that diagonal
/// to its square end: `20 - u` long at `u` in from its outer side.
#[test]
fn an_l_slab_bevelled_past_half_its_legs_turns_its_ridge_at_the_reflex_corner() {
    let (w, t, d) = (3.0, 4.0, 2.0);
    let (area, moment) = ridged_web(w, t, d);
    let want = 2.0 * (area * 20.0 - moment);
    let outline = [
        (0.0, 0.0),
        (20.0, 0.0),
        (20.0, w),
        (w, w),
        (w, 20.0),
        (0.0, 20.0),
    ];
    let mut first: Option<Vec<[u64; 3]>> = None;
    for seed in 0..3_u64 {
        let mut m = Model::new();
        let s = slab(&mut m, &outline, t);
        let top = explore_unique(&m, &s, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| ogeom_algo::face_normal(&m, f, T).unwrap().1.z > 0.5)
            .unwrap();
        let specs: Vec<(Shape, Chamfer)> = top_edges(&m, &s, t, d)
            .into_iter()
            .filter(|(e, _)| {
                let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                let point = |v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;
                let (p, q) = (point(&a), point(&b));
                // Not an end: the ends lie at x = 20 or y = 20.
                !((p.x - 20.0).abs() < 1e-9 && (q.x - 20.0).abs() < 1e-9
                    || (p.y - 20.0).abs() < 1e-9 && (q.y - 20.0).abs() < 1e-9)
            })
            .collect();
        assert_eq!(specs.len(), 4);
        let specs: Vec<(Shape, Chamfer)> = shuffled(specs.len(), seed)
            .iter()
            .map(|i| specs[*i].clone())
            .collect();
        let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        assert!(!folded(&m, &built.shape), "folded");
        let got = volume(&m, &built.shape);
        assert!(
            (got - want).abs() <= 1e-9 * want,
            "volume {got} against {want}"
        );
        // The bottom, the six sides and the four bevels.
        let faces = explore_unique(&m, &built.shape, ShapeType::Face)
            .unwrap()
            .len();
        assert_eq!(faces, 11, "faces");
        assert!(built.history.is_deleted(&top), "the top face left");
        let k = key(&m, &built.shape).1;
        match &first {
            None => first = Some(k),
            Some(held) => assert!(*held == k, "order {seed} differs"),
        }
    }
}

/// Where three planes meet.
fn meet(
    (a, da): (ogeom_math::Vector, f64),
    (b, db): (ogeom_math::Vector, f64),
    (c, dc): (ogeom_math::Vector, f64),
) -> ogeom_math::Vector {
    (b.cross(c) * da + c.cross(a) * db + a.cross(b) * dc) / a.dot(b.cross(c))
}

/// The pyramid over the quadrilateral `base` (counter-clockwise seen from
/// above, on z = 0) up to `apex`, its four slant edges bevelled `d`, where
/// the apex's four corners (each where a face's two bevels' contacts
/// cross) lie off one plane, so no facet closes the apex. It closes on the
/// bilinear patch through them, its sides the bevels' ends: valid,
/// nothing folded, a face for each face, bevel and the patch, one of them
/// a B-spline, the apex generating the patch. The bilinear patch through a
/// tetrahedron's four corners halves it, so the solid holds the corner
/// closed by the two triangles that keep it convex (the polytope of the
/// faces', the bevels' and those triangles' planes) less half the
/// tetrahedron of the corners. The solid's centroid, returned.
fn chamfered_apex(
    base: [ogeom_math::Vector; 4],
    apex: ogeom_math::Vector,
    d: f64,
) -> ogeom_math::Point {
    use ogeom_math::Vector;
    let h = apex.z;
    let centre = (base.iter().fold(apex, |s, p| s + *p)) / 5.0;
    // Side `k` through base corners `k` and `k + 1` and the apex.
    let side = |k: usize| {
        let (p, q) = (base[k], base[(k + 1) % 4]);
        let mut n = (q - p).cross(apex - p);
        n = n / n.magnitude();
        if n.dot(p - centre) < 0.0 {
            n = -n;
        }
        (n, n.dot(p))
    };
    // Slant edge `k`, from the apex to base corner `k`, between sides
    // `k - 1` and `k`.
    let bevel = |k: usize| {
        let u = base[k] - apex;
        bevel_plane(
            apex,
            u / u.magnitude(),
            (side((k + 3) % 4).0, d),
            (side(k).0, d),
        )
    };
    let mut planes: Vec<(Vector, f64)> = vec![(Vector::new(0.0, 0.0, -1.0), 0.0)];
    planes.extend((0..4).map(side));
    planes.extend((0..4).map(bevel));
    let corners: Vec<Vector> = (0..4)
        .map(|k| meet(side(k), bevel(k), bevel((k + 1) % 4)))
        .collect();
    let low = corners.iter().map(|c| c.z).fold(f64::INFINITY, f64::min);
    let high = corners
        .iter()
        .map(|c| c.z)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        high - low > 1e-3,
        "the corners lie off one plane: {corners:?}"
    );
    // Closed by the triangles on either diagonal, each facing away from
    // the solid: the convex closure holds more.
    let outward = |p: Vector, q: Vector, r: Vector| {
        let mut n = (q - p).cross(r - p);
        n = n / n.magnitude();
        if n.dot(p - centre) < 0.0 {
            n = -n;
        }
        (n, n.dot(p))
    };
    let closed = |i: usize| {
        let c = |k: usize| corners[(i + k) % 4];
        let mut with = planes.clone();
        with.push(outward(c(0), c(1), c(2)));
        with.push(outward(c(0), c(2), c(3)));
        convex_polytope(&with).0
    };
    let convex = closed(0).max(closed(1));
    let tetra = (corners[1] - corners[0])
        .dot((corners[2] - corners[0]).cross(corners[3] - corners[0]))
        .abs()
        / 6.0;
    let want = convex - tetra / 2.0;

    let mut m = Model::new();
    let points: Vec<Point> = base
        .iter()
        .chain([&apex])
        .map(|v| Point::new(v.x, v.y, v.z))
        .collect();
    let rings = vec![
        vec![0, 3, 2, 1],
        vec![0, 1, 4],
        vec![1, 2, 4],
        vec![2, 3, 4],
        vec![3, 0, 4],
    ];
    let s = ogeom_algo::make_polyhedron(&mut m, &points, &rings, T)
        .unwrap()
        .shape;
    let top = explore_unique(&m, &s, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| (m.node(v).unwrap().data().as_vertex().unwrap().point.z - h).abs() < 1e-9)
        .unwrap();
    let specs: Vec<(Shape, Chamfer)> = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (p, q) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            p.is_same(&top) || q.is_same(&top)
        })
        .map(|e| (e, Chamfer::Symmetric(d)))
        .collect();
    assert_eq!(specs.len(), 4);
    let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert!(!folded(&m, &built.shape), "folded");
    let faces = explore_unique(&m, &built.shape, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 10, "faces");
    let patches = faces
        .iter()
        .filter(|f| {
            let data = m.node(f).unwrap().data().as_face().unwrap();
            matches!(
                m.geometry().surface(data.surface),
                Some(ogeom_geom::SurfaceGeometry::BSpline(_))
            )
        })
        .count();
    assert_eq!(patches, 1, "one patch");
    assert_eq!(built.history.generated(&top).len(), 1, "the apex's patch");
    let fine = ogeom_mesh::Deflection {
        chord: 1e-5,
        ..ogeom_mesh::Deflection::default()
    };
    let props = ogeom_algo::volume_properties(&m, &built.shape, fine, T).unwrap();
    assert!(
        (props.mass - want).abs() <= 1e-7 * want,
        "volume {} against {want}",
        props.mass
    );
    props.centre
}

/// The pyramid over the rectangle 8 by 6, height 4, its slant edges
/// bevelled 1/2 ([`chamfered_apex`]): the long faces' corners and the
/// short faces' stand at two heights. It keeps the pyramid's two mirrors:
/// its centroid stands on the axis.
#[test]
fn a_rectangular_pyramid_s_apex_chamfered_closes_on_one_patch() {
    use ogeom_math::Vector;
    let (a, b) = (4.0, 3.0);
    let base = [
        Vector::new(-a, -b, 0.0),
        Vector::new(a, -b, 0.0),
        Vector::new(a, b, 0.0),
        Vector::new(-a, b, 0.0),
    ];
    let centre = chamfered_apex(base, Vector::new(0.0, 0.0, 4.0), 0.5);
    assert!(centre.x.abs() < 1e-6 && centre.y.abs() < 1e-6, "{centre:?}");
}

/// The pyramid over the square 20 by 20 with its apex moved off the middle
/// to (5, 8, 12), its slant edges bevelled 2 ([`chamfered_apex`]): no two
/// of the apex's corners mirror each other.
#[test]
fn a_leaning_pyramid_s_apex_chamfered_closes_on_one_patch() {
    use ogeom_math::Vector;
    let base = [
        Vector::new(0.0, 0.0, 0.0),
        Vector::new(20.0, 0.0, 0.0),
        Vector::new(20.0, 20.0, 0.0),
        Vector::new(0.0, 20.0, 0.0),
    ];
    chamfered_apex(base, Vector::new(5.0, 8.0, 12.0), 2.0);
}

/// The volume of the pyramid over the polygon `base` (counter-clockwise
/// seen from above, on z = 0) up to `apex`, its slant edges bevelled `d`
/// and the apex closed by the corners' Wachspress blend over the regular
/// polygon of as many corners: by the divergence theorem about the solid's
/// middle, each planar face's polygon from the planes meeting at its
/// corners, the patch's part summed over `grid` by `grid` cells of each of
/// its sectors.
fn wachspress_apex_volume(
    base: &[ogeom_math::Vector],
    apex: ogeom_math::Vector,
    d: f64,
    grid: usize,
) -> f64 {
    use ogeom_math::Vector;
    let n = base.len();
    #[expect(clippy::cast_precision_loss, reason = "a handful of corners")]
    let centre = base.iter().fold(apex, |s, p| s + *p) / (n + 1) as f64;
    let side = |k: usize| {
        let (p, q) = (base[k], base[(k + 1) % n]);
        let mut m = (q - p).cross(apex - p);
        m = m / m.magnitude();
        if m.dot(p - centre) < 0.0 {
            m = -m;
        }
        (m, m.dot(p))
    };
    let bevel = |k: usize| {
        let u = base[k] - apex;
        bevel_plane(
            apex,
            u / u.magnitude(),
            (side((k + n - 1) % n).0, d),
            (side(k).0, d),
        )
    };
    let floor = (Vector::new(0.0, 0.0, -1.0), 0.0);
    let corner = |k: usize| meet(side(k), bevel(k), bevel((k + 1) % n));
    let low = |k: usize| meet(floor, side(k), bevel(k));
    let high = |k: usize| meet(floor, side(k), bevel((k + 1) % n));
    // A polygon's share, turned to face away from the middle.
    let share = |ring: Vec<Vector>| {
        let mut newell = Vector::ZERO;
        let mut mid = Vector::ZERO;
        for k in 0..ring.len() {
            newell += ring[k].cross(ring[(k + 1) % ring.len()]);
            mid += ring[k];
        }
        #[expect(clippy::cast_precision_loss, reason = "a few corners")]
        let mid = mid / ring.len() as f64;
        let sign = if newell.dot(mid - centre) < 0.0 {
            -1.0
        } else {
            1.0
        };
        let mut sum = 0.0;
        for k in 1..ring.len() - 1 {
            sum += (ring[0] - centre).dot((ring[k] - centre).cross(ring[k + 1] - centre));
        }
        sign * sum / 6.0
    };
    let mut volume = share((0..n).flat_map(|k| [low(k), high(k)]).collect());
    for k in 0..n {
        volume += share(vec![low(k), high(k), corner(k)]);
        volume += share(vec![
            high((k + n - 1) % n),
            low(k),
            corner(k),
            corner((k + n - 1) % n),
        ]);
    }
    // The patch: the corners blended by Wachspress's coordinates.
    #[expect(clippy::cast_precision_loss, reason = "a handful of corners")]
    let polygon: Vec<(f64, f64)> = (0..n)
        .map(|k| {
            let t = core::f64::consts::TAU * k as f64 / n as f64;
            (t.cos(), t.sin())
        })
        .collect();
    let corners: Vec<Vector> = (0..n).map(corner).collect();
    let patch = |x: (f64, f64)| -> Vector {
        let area = |j: usize| {
            let (a, b) = (polygon[j], polygon[(j + 1) % n]);
            ((a.0 - x.0) * (b.1 - x.1) - (a.1 - x.1) * (b.0 - x.0)).max(0.0)
        };
        let w: Vec<f64> = (0..n)
            .map(|i| {
                (0..n)
                    .filter(|j| *j != i && *j != (i + n - 1) % n)
                    .map(area)
                    .product()
            })
            .collect();
        let total: f64 = w.iter().sum();
        corners
            .iter()
            .zip(&w)
            .fold(Vector::ZERO, |s, (c, wi)| s + *c * (wi / total))
    };
    let mut cells = 0.0;
    #[expect(clippy::cast_precision_loss, reason = "a fine grid")]
    let g = grid as f64;
    for k in 0..n {
        let (a, b) = (polygon[k], polygon[(k + 1) % n]);
        // The sector from the polygon's middle to side `k`, in barycentric
        // cells: `x = (i a + j b) / grid` over the triangle.
        let at = |i: usize, j: usize| {
            #[expect(clippy::cast_precision_loss, reason = "a fine grid")]
            let (s, t) = (i as f64 / g, j as f64 / g);
            patch((s * a.0 + t * b.0, s * a.1 + t * b.1))
        };
        for i in 0..grid {
            for j in 0..grid - i {
                let tri = |p: Vector, q: Vector, r: Vector| {
                    (p - centre).dot((q - centre).cross(r - centre)) / 6.0
                };
                cells += tri(at(i, j), at(i + 1, j), at(i, j + 1));
                if i + j + 1 < grid {
                    cells += tri(at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
                }
            }
        }
    }
    // The patch faces away from the base, toward the apex.
    let up = (patch((0.0, 0.0)) - centre).dot(apex - centre);
    volume
        + if up * cells < 0.0 {
            -cells.abs()
        } else {
            cells.abs()
        }
}

/// [`wachspress_apex_volume`] on a square's corners is the bilinear
/// patch's, so it agrees with [`chamfered_apex`]'s exact oracle on the
/// rectangular pyramid; the oracle checked, the pyramid over an uneven
/// pentagon with its apex off the middle, its five slant edges bevelled
/// 1/2, closes on the five corners' patch, a face on each side's sector:
/// valid, nothing folded, a face for each face and bevel and five for the
/// patch, all five B-splines generated by the apex, the volume the
/// oracle's.
#[test]
fn a_pentagonal_pyramid_s_apex_chamfered_closes_on_one_patch_of_five_sectors() {
    use ogeom_math::Vector;
    // The oracle against the exact one, on four corners.
    let rectangle = [
        Vector::new(-4.0, -3.0, 0.0),
        Vector::new(4.0, -3.0, 0.0),
        Vector::new(4.0, 3.0, 0.0),
        Vector::new(-4.0, 3.0, 0.0),
    ];
    let exact = {
        let mut m = Model::new();
        let _ = chamfered_apex(rectangle, Vector::new(0.0, 0.0, 4.0), 0.5);
        let points: Vec<Point> = rectangle
            .iter()
            .chain([&Vector::new(0.0, 0.0, 4.0)])
            .map(|v| Point::new(v.x, v.y, v.z))
            .collect();
        let rings = vec![
            vec![0, 3, 2, 1],
            vec![0, 1, 4],
            vec![1, 2, 4],
            vec![2, 3, 4],
            vec![3, 0, 4],
        ];
        let s = ogeom_algo::make_polyhedron(&mut m, &points, &rings, T)
            .unwrap()
            .shape;
        let specs: Vec<(Shape, Chamfer)> = explore_unique(&m, &s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let (p, q) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                let z = |v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point.z;
                z(&p) > 1.0 || z(&q) > 1.0
            })
            .map(|e| (e, Chamfer::Symmetric(0.5)))
            .collect();
        let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
        volume(&m, &built.shape)
    };
    let oracle = wachspress_apex_volume(&rectangle, Vector::new(0.0, 0.0, 4.0), 0.5, 400);
    assert!(
        (oracle - exact).abs() <= 1e-7 * exact,
        "oracle {oracle} against {exact}"
    );

    let base = [
        Vector::new(-3.0, -2.0, 0.0),
        Vector::new(3.0, -2.0, 0.0),
        Vector::new(4.0, 1.5, 0.0),
        Vector::new(0.0, 4.0, 0.0),
        Vector::new(-4.0, 1.5, 0.0),
    ];
    let apex = Vector::new(0.5, 0.3, 5.0);
    let d = 0.5;
    let want = wachspress_apex_volume(&base, apex, d, 400);
    let mut m = Model::new();
    let points: Vec<Point> = base
        .iter()
        .chain([&apex])
        .map(|v| Point::new(v.x, v.y, v.z))
        .collect();
    let rings = vec![
        vec![0, 4, 3, 2, 1],
        vec![0, 1, 5],
        vec![1, 2, 5],
        vec![2, 3, 5],
        vec![3, 4, 5],
        vec![4, 0, 5],
    ];
    let s = ogeom_algo::make_polyhedron(&mut m, &points, &rings, T)
        .unwrap()
        .shape;
    let top = explore_unique(&m, &s, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| (m.node(v).unwrap().data().as_vertex().unwrap().point.z - apex.z).abs() < 1e-9)
        .unwrap();
    let specs: Vec<(Shape, Chamfer)> = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (p, q) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            p.is_same(&top) || q.is_same(&top)
        })
        .map(|e| (e, Chamfer::Symmetric(d)))
        .collect();
    assert_eq!(specs.len(), 5);
    let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert!(!folded(&m, &built.shape), "folded");
    let faces = explore_unique(&m, &built.shape, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 16, "faces");
    let splines: Vec<&Shape> = faces
        .iter()
        .filter(|f| {
            let data = m.node(f).unwrap().data().as_face().unwrap();
            matches!(
                m.geometry().surface(data.surface),
                Some(ogeom_geom::SurfaceGeometry::BSpline(_))
            )
        })
        .collect();
    assert_eq!(splines.len(), 5, "the patch's sectors");
    assert_eq!(built.history.generated(&top).len(), 5, "the apex's patch");
    let fine = ogeom_mesh::Deflection {
        chord: 1e-5,
        ..ogeom_mesh::Deflection::default()
    };
    let got = ogeom_algo::volume_properties(&m, &built.shape, fine, T)
        .unwrap()
        .mass;
    assert!(
        (got - want).abs() <= 1e-7 * want,
        "volume {got} against {want}"
    );
}

/// The prism of [`prism`]'s profile built as a sweep: the profile's face
/// in x = 0 swept by `make_prism`, unbaked, so its sides are the profile's
/// lines swept along x (flat faces on extrusion surfaces) and its far cap
/// is the near one placed.
fn swept_prism(m: &mut Model, profile: &[(f64, f64)], length: f64) -> Shape {
    let corners: Vec<Shape> = profile
        .iter()
        .map(|&(y, z)| ogeom_algo::make_vertex(m, Point::new(0.0, y, z)).shape)
        .collect();
    let n = corners.len();
    let edges: Vec<Shape> = (0..n)
        .map(|k| {
            let (p, q) = (profile[k], profile[(k + 1) % n]);
            let (a, b) = (Point::new(0.0, p.0, p.1), Point::new(0.0, q.0, q.1));
            let line = ogeom_geom::LineCurve::segment(a, b, T).unwrap();
            let range = ogeom_geom::Curve3d::domain(&line);
            ogeom_algo::make_edge_between(
                m,
                line.into(),
                range,
                &corners[k],
                &corners[(k + 1) % n],
                T,
            )
            .unwrap()
            .shape
        })
        .collect();
    let wire = ogeom_algo::make_wire(m, &edges, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        ogeom_math::Plane::through(Point::ORIGIN, ogeom_math::Direction::X),
        (-1000.0, 1000.0),
        (-1000.0, 1000.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom_algo::make_prism(m, &face, ogeom_math::Vector::new(length, 0.0, 0.0), T)
        .unwrap()
        .shape
}

/// The corpus's parts built as unbaked sweeps, every edge chamfered:
/// valid, a face for every face, edge and vertex, the sections' closed
/// form, every edge of the sweep as named (the far cap's through its
/// placement) generating its bevel.
#[test]
fn a_swept_prism_with_its_far_cap_placed_chamfers_every_edge_exactly() {
    for (name, profile, length) in corpus() {
        let d = 1.0;
        let want = chamfered_prism_volume(&profile, length, d);
        let mut m = Model::new();
        let s = swept_prism(&mut m, &profile, length);
        let (f, e, v) = (
            explore_unique(&m, &s, ShapeType::Face).unwrap().len(),
            explore_unique(&m, &s, ShapeType::Edge).unwrap().len(),
            explore_unique(&m, &s, ShapeType::Vertex).unwrap().len(),
        );
        let specs = all_edges(&m, &s, d);
        let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T)
            .unwrap_or_else(|err| panic!("{name}: {err}"));
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{name}: {:?}", found.problems);
        assert!(!folded(&m, &built.shape), "{name}: folded");
        let faces = explore_unique(&m, &built.shape, ShapeType::Face)
            .unwrap()
            .len();
        assert_eq!(faces, f + e + v, "{name}: faces");
        let got = volume(&m, &built.shape);
        assert!(
            (got - want).abs() <= 1e-9 * want,
            "{name}: volume {got} against {want}"
        );
        for (edge, _) in &specs {
            assert_eq!(
                built.history.generated(edge).len(),
                1,
                "{name}: an edge without its bevel"
            );
        }
    }
}

/// A box's every edge chamfered with the mitre corner: at each corner the
/// three bevels run on to the point their planes share, so the box keeps,
/// past the facet chamfer's closed form, each corner's tetrahedron between
/// that point and the facet's three corners, `d³ / 12`; no corner has a
/// facet, a face for every face and edge of the box. The same solid to the
/// bit under three orders.
#[test]
fn a_box_mitred_on_every_edge_keeps_each_corner_s_tetrahedron() {
    let profile = vec![(0.0, 0.0), (40.0, 0.0), (40.0, 20.0), (0.0, 20.0)];
    for d in [0.5, 2.0] {
        let want = chamfered_prism_volume(&profile, 30.0, d) + 8.0 * d * d * d / 12.0;
        let mut first: Option<(u64, Vec<[u64; 3]>)> = None;
        for seed in 0..3_u64 {
            let mut m = Model::new();
            let s = prism(&mut m, &profile, 30.0);
            let specs = all_edges(&m, &s, d);
            let specs: Vec<(Shape, Chamfer)> = shuffled(specs.len(), seed)
                .iter()
                .map(|i| specs[*i].clone())
                .collect();
            let built =
                ogeom_fillet::chamfer_edges_cornered(&mut m, &s, &specs, ChamferCorner::Point, T)
                    .unwrap();
            let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(found.is_valid(), "d {d}: {:?}", found.problems);
            assert!(!folded(&m, &built.shape), "d {d}: folded");
            let faces = explore_unique(&m, &built.shape, ShapeType::Face)
                .unwrap()
                .len();
            assert_eq!(faces, 6 + 12, "d {d}: faces");
            let got = volume(&m, &built.shape);
            assert!(
                (got - want).abs() <= 1e-9 * want,
                "d {d}: volume {got} against {want}"
            );
            let k = key(&m, &built.shape);
            match &first {
                None => first = Some(k),
                Some(held) => assert!(*held == k, "d {d}: order {seed} differs"),
            }
        }
    }
}

/// A compound of two solids apart, a box and an L prism, every edge
/// chamfered: each solid chamfered as it would be alone, the result a
/// compound of the two, valid, a face for every face, edge and vertex of
/// each, the two closed forms summed, every edge generating its bevel.
#[test]
fn a_compound_of_two_solids_chamfers_each_as_alone() {
    let parts = corpus();
    let (_, box_profile, box_length) = parts.iter().find(|p| p.0 == "box").unwrap().clone();
    let (_, l_profile, l_length) = parts.iter().find(|p| p.0 == "L").unwrap().clone();
    let d = 1.0;
    let want = chamfered_prism_volume(&box_profile, box_length, d)
        + chamfered_prism_volume(&l_profile, l_length, d);
    let mut m = Model::new();
    let a = prism(&mut m, &box_profile, box_length);
    // The L moved clear of the box, 100 up.
    let l: Vec<(f64, f64)> = l_profile.iter().map(|&(y, z)| (y, z + 100.0)).collect();
    let b = prism(&mut m, &l, l_length);
    let counts = |m: &Model, s: &Shape| {
        explore_unique(m, s, ShapeType::Face).unwrap().len()
            + explore_unique(m, s, ShapeType::Edge).unwrap().len()
            + explore_unique(m, s, ShapeType::Vertex).unwrap().len()
    };
    let expected_faces = counts(&m, &a) + counts(&m, &b);
    let compound = ogeom_algo::make_compound(&mut m, &[a, b]).unwrap().shape;
    let specs = all_edges(&m, &compound, d);
    let built = ogeom_fillet::chamfer_edges_with(&mut m, &compound, &specs, T).unwrap();
    assert_eq!(
        m.kind_of(&built.shape).unwrap(),
        ShapeType::Compound,
        "a compound back"
    );
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert!(!folded(&m, &built.shape), "folded");
    let faces = explore_unique(&m, &built.shape, ShapeType::Face)
        .unwrap()
        .len();
    assert_eq!(faces, expected_faces, "faces");
    let got = volume(&m, &built.shape);
    assert!(
        (got - want).abs() <= 1e-9 * want,
        "volume {got} against {want}"
    );
    for (edge, _) in &specs {
        assert_eq!(
            built.history.generated(edge).len(),
            1,
            "an edge without its bevel"
        );
    }
}

/// A box 30 by 40 by 20 holding a void 10 by 6 by 8 well inside it, cut by
/// the boolean (a solid of two shells), every edge of its outer shell
/// chamfered: the outer shell as the box alone, the void kept as it was.
/// Valid, two shells, a face for every face, edge and vertex of the box
/// and the void's six, the box's closed form less the void, every edge
/// generating its bevel.
#[test]
fn a_box_with_a_void_chamfers_its_outside_and_keeps_the_void() {
    let profile = vec![(0.0, 0.0), (40.0, 0.0), (40.0, 20.0), (0.0, 20.0)];
    let d = 1.0;
    let want = chamfered_prism_volume(&profile, 30.0, d) - 10.0 * 6.0 * 8.0;
    let mut m = Model::new();
    let outside = prism(&mut m, &profile, 30.0);
    let edges = explore_unique(&m, &outside, ShapeType::Edge).unwrap();
    let (f, e, v) = (
        explore_unique(&m, &outside, ShapeType::Face).unwrap().len(),
        edges.len(),
        explore_unique(&m, &outside, ShapeType::Vertex)
            .unwrap()
            .len(),
    );
    let void = block(&mut m, (10.0, 15.0, 7.0), (10.0, 6.0, 8.0));
    let s = ogeom_bool::cut(&mut m, &outside, &void, T).unwrap().shape;
    assert_eq!(m.children_of(&s).unwrap().len(), 2, "a solid of two shells");
    let outer_edges: Vec<Shape> = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|edge| {
            let b = ogeom_algo::shape_bounds(&m, edge, T).unwrap();
            let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
            // On the box's outside: some coordinate at its bounds.
            [lo.x, hi.x]
                .iter()
                .any(|x| x.abs() < 1e-6 || (x - 30.0).abs() < 1e-6)
                || [lo.y, hi.y]
                    .iter()
                    .any(|y| y.abs() < 1e-6 || (y - 40.0).abs() < 1e-6)
                || [lo.z, hi.z]
                    .iter()
                    .any(|z| z.abs() < 1e-6 || (z - 20.0).abs() < 1e-6)
        })
        .collect();
    assert_eq!(outer_edges.len(), e, "the box's own edges");
    let specs: Vec<(Shape, Chamfer)> = outer_edges
        .into_iter()
        .map(|edge| (edge, Chamfer::Symmetric(d)))
        .collect();
    let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert_eq!(
        m.children_of(&built.shape).unwrap().len(),
        2,
        "the void kept"
    );
    let faces = explore_unique(&m, &built.shape, ShapeType::Face)
        .unwrap()
        .len();
    assert_eq!(faces, f + e + v + 6, "faces");
    let got = volume(&m, &built.shape);
    assert!(
        (got - want).abs() <= 1e-9 * want,
        "volume {got} against {want}"
    );
    for (edge, _) in &specs {
        assert_eq!(
            built.history.generated(edge).len(),
            1,
            "an edge without its bevel"
        );
    }
}
