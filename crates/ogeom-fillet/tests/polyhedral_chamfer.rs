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
