//! Soundness probes round the corner fills (not committed).
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Blend, Chamfer, ChamferCorner, Setbacks};
use ogeom_geom::{Continuity, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn vol(m: &Model, s: &Shape) -> f64 {
    let d = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    ogeom_algo::volume_properties(m, s, d, T).unwrap().mass
}

fn point_of(m: &Model, v: &Shape) -> Point {
    let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
    v.transform(m.datums()).unwrap().apply(p)
}

/// A prism on a parallelogram of sides 20 (along x) and 20 (at `angle`
/// degrees from x), raised along (`tilt`, 0, 10). The corner probed is the
/// top of the origin's edge; returns the solid, that vertex, and its three
/// edges: along x, along the slanted side, and the raised edge.
fn skewed(m: &mut Model, angle: f64, tilt: f64) -> (Shape, Shape, [Shape; 3]) {
    let a = angle.to_radians();
    let side = Vector::new(a.cos(), a.sin(), 0.0) * 20.0;
    let pts = [
        Point::ORIGIN,
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 0.0) + side,
        Point::ORIGIN + side,
    ];
    let wire = ogeom_algo::make_polygon(m, &pts, true, T).unwrap().shape;
    let plane: SurfaceGeometry = ogeom_geom::PlaneSurface::over(
        ogeom_math::Plane::through(Point::ORIGIN, Direction::Z),
        (-1e3, 1e3),
        (-1e3, 1e3),
    )
    .unwrap()
    .into();
    let face = ogeom_algo::make_face(m, plane, &[wire], T).unwrap().shape;
    let up = Vector::new(tilt, 0.0, 10.0);
    let prism = ogeom_algo::make_prism(m, &face, up, T).unwrap().shape;
    let solid = explore_unique(m, &prism, ShapeType::Solid)
        .unwrap()
        .into_iter()
        .next()
        .unwrap_or(prism);
    let top = Point::ORIGIN + up;
    let vertex = explore_unique(m, &solid, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| point_of(m, v).distance(top) < 1e-9)
        .unwrap();
    let dirs = [Vector::new(1.0, 0.0, 0.0), side / 20.0, -up / up.magnitude()];
    let mut edges: Vec<Shape> = Vec::new();
    for d in dirs {
        let e = explore_unique(m, &solid, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let Some((p, q)) = ogeom_algo::edge_vertices(m, e).unwrap() else {
                    return false;
                };
                let (p, q) = (point_of(m, &p), point_of(m, &q));
                let other = if p.distance(top) < 1e-9 {
                    q
                } else if q.distance(top) < 1e-9 {
                    p
                } else {
                    return false;
                };
                let u = other - top;
                u.cross(d).magnitude() < 1e-9 * u.magnitude() && u.dot(d) > 0.0
            })
            .unwrap();
        edges.push(e);
    }
    (solid, vertex, [edges[0].clone(), edges[1].clone(), edges[2].clone()])
}

/// Whether case `k` of a probe is the one asked (`PROBE_CASE`), or all.
fn wanted(k: &mut usize) -> bool {
    let here = *k;
    *k += 1;
    std::env::var("PROBE_CASE").map_or(true, |v| v.parse::<usize>().unwrap() == here)
}

/// Findings, printed and gathered.
#[derive(Default)]
struct Log(Vec<String>);
impl Log {
    fn bad(&mut self, s: String) {
        eprintln!("FINDING {s}");
        self.0.push(s);
    }
}

/// The generic checks on a result: valid; where the original is a skewed
/// prism (`prism`), no point of the result's mesh outside it; where a
/// reference is given, the two differ only within `reach` of `at`; and no
/// fitted face's normal folds across a triangle of its mesh.
fn audit(
    log: &mut Log,
    name: &str,
    m: &mut Model,
    prism: Option<(f64, f64)>,
    result: &Shape,
    reference: Option<&Shape>,
    at: Point,
    reach: f64,
) -> f64 {
    let d = ogeom_algo::check(m, result, T).unwrap();
    if !d.is_valid() {
        log.bad(format!("{name}: invalid {:?}", d.problems.iter().take(3).collect::<Vec<_>>()));
    }
    let v = vol(m, result);
    let deflection = ogeom_mesh::Deflection { chord: 1e-3, ..ogeom_mesh::Deflection::default() };
    for f in explore_unique(m, result, ShapeType::Face).unwrap() {
        let Ok(tri) = ogeom_mesh::triangulate_face(m, &f, deflection, T) else {
            eprintln!("NOTE {name}: a face did not mesh");
            continue;
        };
        let data = m.node(&f).unwrap().data().as_face().unwrap();
        let surface = m.geometry().surface(data.surface).unwrap().clone();
        let kind = format!("{:?}", core::mem::discriminant(&surface));
        if let Some((angle, tilt)) = prism {
            let (worst, p) = tri
                .positions
                .iter()
                .map(|p| (outside_skewed(angle, tilt, *p), *p))
                .fold((f64::NEG_INFINITY, Point::ORIGIN), |a, b| if b.0 > a.0 { b } else { a });
            if worst > 2e-4 {
                log.bad(format!("{name}: a face ({kind}) stands {worst} outside the original at {p:?}"));
            }
        }
        if let SurfaceGeometry::BSpline(s) = &surface {
            let mut worst = 1.0_f64;
            for t in &tri.triangles {
                let ns: Vec<Option<Vector>> = t
                    .iter()
                    .map(|&i| {
                        let (u, w) = tri.parameters[i as usize];
                        let (du, dv) = s.d1_at(u, w, T).unwrap();
                        let c = du.cross(dv);
                        (c.magnitude() > 1e-9).then(|| c / c.magnitude())
                    })
                    .collect();
                for (a, b) in [(0, 1), (1, 2), (0, 2)] {
                    if let (Some(x), Some(y)) = (ns[a], ns[b]) {
                        worst = worst.min(x.dot(y));
                    }
                }
            }
            if worst < 0.0 {
                log.bad(format!("{name}: a fitted face's normal folds across a triangle (cos {worst})"));
            }
        }
    }
    if let Some(r) = reference {
        for (label, a, b) in [("result-ref", result, r), ("ref-result", r, result)] {
            match ogeom_bool::cut(m, a, b, T) {
                Ok(diff) => {
                    for vx in explore_unique(m, &diff.shape, ShapeType::Vertex).unwrap() {
                        let p = point_of(m, &vx);
                        if p.distance(at) > reach {
                            let dv = vol(m, &diff.shape);
                            if dv > 1e-6 {
                                log.bad(format!(
                                    "{name}: {label} differs at {p:?}, {} from the corner (reach {reach}), diff volume {dv}",
                                    p.distance(at)
                                ));
                            }
                            break;
                        }
                    }
                }
                Err(e) => eprintln!("NOTE {name}: {label} boolean refused: {e}"),
            }
        }
    }
    v
}

fn corners(m: &mut Model) -> Vec<(String, (f64, f64), Shape, Shape, [Shape; 3])> {
    let mut out = Vec::new();
    for (angle, tilt) in [(90.0, 0.0), (60.0, 0.0), (120.0, 0.0), (90.0, 3.0), (75.0, -2.5)] {
        let (s, v, e) = skewed(m, angle, tilt);
        out.push((format!("corner{angle}/{tilt}"), (angle, tilt), s, v, e));
    }
    out
}

#[test]
#[ignore = "probe"]
fn probe_sharp() {
    let mut log = Log::default();
    let cases: &[(usize, f64, [f64; 3])] = &[
        (2, 2.0, [5.0, 5.0, 5.0]),
        (2, 1.0, [3.0, 6.0, 4.0]),
        (2, 3.0, [7.0, 7.0, 9.0]),
        (0, 1.5, [4.0, 5.0, 3.5]),
        (1, 0.5, [2.0, 1.5, 3.0]),
        (2, 2.0, [4.2, 4.2, 4.5]),
        (0, 2.0, [8.0, 6.0, 5.0]),
    ];
    let mut m = Model::new();
    let mut k = 0;
    for (cname, prism, solid, vertex, edges) in corners(&mut m) {
        for &(sharp, radius, backs) in cases {
            if !wanted(&mut k) {
                continue;
            }
            let name = format!("sharp {cname} s{sharp} r{radius} {backs:?}");
            let setbacks = Setbacks {
                per_edge: (0..3).map(|k| (edges[k].clone(), backs[k])).collect(),
                default: 5.0,
                continuity: Continuity::G1,
            };
            let at = point_of(&m, &vertex);
            let r = ogeom_fillet::round_vertex_setback_sharp(
                &mut m, &solid, &vertex, &edges[sharp], radius, &setbacks, T,
            );
            let r = match r {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("REFUSED {name}: {e}");
                    continue;
                }
            };
            let rounded: Vec<(Shape, f64)> =
                (0..3).filter(|&k| k != sharp).map(|k| (edges[k].clone(), radius)).collect();
            let reference = ogeom_fillet::fillet_edges_with(&mut m, &solid, &rounded, T).ok();
            let reach = backs.iter().copied().fold(0.0, f64::max) + 3.0 * radius;
            let v = audit(&mut log, &name, &mut m, Some(prism), &r.shape, reference.as_ref().map(|b| &b.shape), at, reach);
            eprintln!("OK {name}: volume {v}");
        }
    }
    eprintln!("FINDINGS {}", log.0.len());
}

fn blends(e: &[Shape; 3], kinds: [usize; 3], sizes: [f64; 3], order: [usize; 3]) -> Vec<(Shape, Blend)> {
    order
        .iter()
        .filter(|&&i| kinds[i] != 2)
        .map(|&i| {
            let b = if kinds[i] == 0 {
                Blend::Chamfer(Chamfer::Symmetric(sizes[i]))
            } else {
                Blend::Round(sizes[i])
            };
            (e[i].clone(), b)
        })
        .collect()
}

const ORDERS: [[usize; 3]; 6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];

#[test]
#[ignore = "probe"]
fn probe_blend() {
    let mut log = Log::default();
    let mut k = 0;
    let kinds_list = [[0, 1, 1], [1, 0, 0], [0, 0, 1], [1, 1, 0], [0, 1, 2], [1, 0, 2], [2, 1, 0], [0, 2, 1]];
    let sizes_list = [[2.0, 2.0, 2.0], [1.0, 3.0, 2.0], [3.0, 1.0, 1.5], [0.5, 0.5, 2.5]];
    let mut m = Model::new();
    for (cname, prism, solid, vertex, edges) in corners(&mut m) {
        let at = point_of(&m, &vertex);
        for kinds in kinds_list {
            for sizes in sizes_list {
                if !wanted(&mut k) {
                    continue;
                }
                let name = format!("blend {cname} {kinds:?} {sizes:?}");
                let mut first: Option<f64> = None;
                for order in ORDERS {
                    let req = blends(&edges, kinds, sizes, order);
                    let r = match ogeom_fillet::blend_edges(&mut m, &solid, &req, ChamferCorner::Point, T) {
                        Ok(r) => r,
                        Err(e) => {
                            eprintln!("REFUSED {name} {order:?}: {e}");
                            break;
                        }
                    };
                    let v = vol(&m, &r.shape);
                    match first {
                        None => {
                            first = Some(v);
                            let reach = sizes.iter().copied().fold(0.0, f64::max) * 6.0;
                            audit(&mut log, &name, &mut m, Some(prism), &r.shape, None, at, reach);
                            eprintln!("OK {name}: volume {v}");
                        }
                        Some(w) if (w - v).abs() > 1e-9 * w => {
                            log.bad(format!("{name}: order {order:?} volume {v} against {w}"));
                        }
                        _ => {}
                    }
                }
            }
        }
    }
    eprintln!("FINDINGS {}", log.0.len());
}

#[test]
#[ignore = "probe"]
fn probe_radii() {
    let mut log = Log::default();
    let mut k = 0;
    let radii_list = [[1.0, 3.0, 2.0], [0.5, 3.0, 2.0], [1.0, 1.0, 3.0], [3.0, 1.0, 0.5], [1.0, 3.5, 2.0], [2.0, 2.5, 3.0]];
    let mut m = Model::new();
    for (cname, prism, solid, vertex, edges) in corners(&mut m) {
        let at = point_of(&m, &vertex);
        for radii in radii_list {
            if !wanted(&mut k) {
                continue;
            }
            let name = format!("radii {cname} {radii:?}");
            let mut first: Option<f64> = None;
            for order in ORDERS {
                let req: Vec<(Shape, f64)> = order.iter().map(|&k| (edges[k].clone(), radii[k])).collect();
                let r = match ogeom_fillet::fillet_edges_with(&mut m, &solid, &req, T) {
                    Ok(r) => r,
                    Err(e) => {
                        eprintln!("REFUSED {name} {order:?}: {e}");
                        break;
                    }
                };
                let v = vol(&m, &r.shape);
                match first {
                    None => {
                        first = Some(v);
                        audit(&mut log, &name, &mut m, Some(prism), &r.shape, None, at, 0.0);
                        let lo = radii.iter().copied().fold(f64::INFINITY, f64::min);
                        let hi = radii.iter().copied().fold(0.0, f64::max);
                        let all = |m: &mut Model, r: f64| {
                            let e: Vec<Shape> = edges.to_vec();
                            ogeom_fillet::fillet_edges(m, &solid, &e, r, T).map(|b| vol(m, &b.shape))
                        };
                        if let (Ok(small), Ok(big)) = (all(&mut m, lo), all(&mut m, hi)) {
                            if !(big < v && v < small) {
                                log.bad(format!("{name}: volume {v} not between {big} and {small}"));
                            }
                        }
                        for f in explore_unique(&m, &r.shape, ShapeType::Face).unwrap() {
                            let data = m.node(&f).unwrap().data().as_face().unwrap();
                            if !matches!(m.geometry().surface(data.surface), Some(SurfaceGeometry::BSpline(_))) {
                                continue;
                            }
                            if ogeom_algo::shape_bounds(&m, &f, T).unwrap().distance_to(at) > 8.0 {
                                continue;
                            }
                            let worst = ogeom_fillet::analyse_blend(&m, &r.shape, &f, 17, T)
                                .unwrap()
                                .iter()
                                .map(|c| c.tangency_error)
                                .fold(0.0, f64::max);
                            if worst > 2e-4 {
                                log.bad(format!("{name}: a patch meets a neighbour {worst} rad off"));
                            }
                        }
                        eprintln!("OK {name}: volume {v}");
                    }
                    Some(w) if (w - v).abs() > 1e-9 * w => {
                        log.bad(format!("{name}: order {order:?} volume {v} against {w}"));
                    }
                    _ => {}
                }
            }
        }
    }
    eprintln!("FINDINGS {}", log.0.len());
}

fn bracket(m: &mut Model, base: f64, up: f64) -> Shape {
    let bx = |m: &mut Model, s: (f64, f64, f64)| {
        ogeom_algo::make_box(m, Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap(), s, T)
            .unwrap()
            .shape
    };
    let a = bx(m, (20.0, 20.0, base));
    let b = bx(m, (up, 20.0, 20.0));
    ogeom_bool::fuse(m, &a, &b, T).unwrap().shape
}

fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
    let d = |e: &Shape| {
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
        (0..=64)
            .map(|i| c.point_at(a + (b - a) * f64::from(i) / 64.0, T).unwrap().distance(q))
            .fold(f64::INFINITY, f64::min)
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

#[test]
#[ignore = "probe"]
fn probe_mixed() {
    let mut log = Log::default();
    let mut k = 0;
    for (base, up) in [(5.0, 5.0), (3.0, 8.0), (8.0, 3.0)] {
        for radius in [0.5, 1.0, 2.0, 2.9] {
            if !wanted(&mut k) {
                continue;
            }
            let name = format!("mixed base{base} up{up} r{radius}");
            let mut first: Option<f64> = None;
            for order in ORDERS {
                let mut m = Model::new();
                let l = bracket(&mut m, base, up);
                let edges = [
                    edge_near(&m, &l, Point::new(12.5, 0.0, base)),
                    edge_near(&m, &l, Point::new(up, 10.0, base)),
                    edge_near(&m, &l, Point::new(up, 0.0, 12.5)),
                ];
                let asked: Vec<Shape> = order.iter().map(|&k| edges[k].clone()).collect();
                let r = match ogeom_fillet::fillet_edges(&mut m, &l, &asked, radius, T) {
                    Ok(r) => r,
                    Err(e) => {
                        eprintln!("REFUSED {name} {order:?}: {e}");
                        break;
                    }
                };
                let v = vol(&m, &r.shape);
                match first {
                    None => {
                        first = Some(v);
                        let at = Point::new(up, 0.0, base);
                        // Stepwise bands, one after another, as the reference.
                        let mut stepwise = l.clone();
                        let mut ok = true;
                        for k in order {
                            let e = [
                                edge_near(&m, &stepwise, Point::new(12.5, 0.0, base)),
                                edge_near(&m, &stepwise, Point::new(up, 10.0, base)),
                                edge_near(&m, &stepwise, Point::new(up, 0.0, 12.5)),
                            ][k]
                                .clone();
                            match ogeom_fillet::fillet_edge(&mut m, &stepwise, &e, radius, T) {
                                Ok(b) => stepwise = b.shape,
                                Err(_) => {
                                    ok = false;
                                    break;
                                }
                            }
                        }
                        let reference = ok.then_some(stepwise);
                        // Concave rounds add material: no inside check.
                        let d = ogeom_algo::check(&m, &r.shape, T).unwrap();
                        if !d.is_valid() {
                            log.bad(format!("{name}: invalid {:?}", d.problems.iter().take(3).collect::<Vec<_>>()));
                        }
                        if let Some(s) = &reference {
                            let reach = radius * 4.0 + 2.0;
                            for (label, a, b) in [("result-ref", &r.shape, s), ("ref-result", s, &r.shape)] {
                                match ogeom_bool::cut(&mut m, a, b, T) {
                                    Ok(diff) => {
                                        for vx in explore_unique(&m, &diff.shape, ShapeType::Vertex).unwrap() {
                                            let p = point_of(&m, &vx);
                                            if p.distance(at) > reach {
                                                log.bad(format!("{name}: {label} differs at {p:?} (reach {reach}), vol {}", vol(&m, &diff.shape)));
                                                break;
                                            }
                                        }
                                    }
                                    Err(e) => eprintln!("NOTE {name}: {label} refused {e}"),
                                }
                            }
                        }
                        audit(&mut log, &name, &mut m, None, &r.shape.clone(), None, at, 0.0);
                        eprintln!("OK {name}: volume {v}");
                    }
                    Some(w) if (w - v).abs() > 1e-9 * w => {
                        log.bad(format!("{name}: order {order:?} volume {v} against {w}"));
                    }
                    _ => {}
                }
            }
        }
    }
    eprintln!("FINDINGS {}", log.0.len());
}

/// How far `p` stands outside the skewed prism (negative inside).
fn outside_skewed(angle: f64, tilt: f64, p: Point) -> f64 {
    let a = angle.to_radians();
    let s = Vector::new(a.cos(), a.sin(), 0.0);
    let up = Vector::new(tilt, 0.0, 10.0);
    let inner = Point::ORIGIN + Vector::new(20.0, 0.0, 0.0) * 0.5 + s * 10.0 + up * 0.5;
    let mut planes: Vec<(Point, Vector)> = vec![
        (Point::ORIGIN, Vector::new(0.0, 0.0, -1.0)),
        (Point::ORIGIN + up, Vector::new(0.0, 0.0, 1.0)),
    ];
    let side_normal = |d: Vector| {
        let n = d.cross(up);
        n / n.magnitude()
    };
    for (on, d) in [
        (Point::ORIGIN, Vector::new(1.0, 0.0, 0.0)),
        (Point::ORIGIN + s * 20.0, Vector::new(1.0, 0.0, 0.0)),
        (Point::ORIGIN, s),
        (Point::new(20.0, 0.0, 0.0), s),
    ] {
        let mut n = side_normal(d);
        if (inner - on).dot(n) > 0.0 {
            n = -n;
        }
        planes.push((on, n));
    }
    planes.iter().map(|(o, n)| (p - *o).dot(*n)).fold(f64::NEG_INFINITY, f64::max)
}

#[test]
#[ignore = "probe"]
fn probe_verify_sharp() {
    let cases: &[(f64, f64, usize, f64, [f64; 3])] = &[
        (60.0, 0.0, 2, 2.0, [5.0, 5.0, 5.0]),
        (120.0, 0.0, 2, 3.0, [7.0, 7.0, 9.0]),
        (60.0, 0.0, 2, 3.0, [7.0, 7.0, 9.0]),
        (90.0, 0.0, 2, 2.0, [5.0, 5.0, 5.0]),
    ];
    for &(angle, tilt, sharp, radius, backs) in cases {
        let mut m = Model::new();
        let (solid, vertex, edges) = skewed(&mut m, angle, tilt);
        let setbacks = Setbacks {
            per_edge: (0..3).map(|k| (edges[k].clone(), backs[k])).collect(),
            default: 5.0,
            continuity: Continuity::G1,
        };
        let r = ogeom_fillet::round_vertex_setback_sharp(&mut m, &solid, &vertex, &edges[sharp], radius, &setbacks, T).unwrap();
        let d = ogeom_mesh::Deflection { chord: 1e-3, ..ogeom_mesh::Deflection::default() };
        for f in explore_unique(&m, &r.shape, ShapeType::Face).unwrap() {
            let data = m.node(&f).unwrap().data().as_face().unwrap();
            let kind = format!("{:?}", core::mem::discriminant(m.geometry().surface(data.surface).unwrap()));
            let tri = ogeom_mesh::triangulate_face(&m, &f, d, T).unwrap();
            let (worst, at) = tri.positions.iter().map(|p| (outside_skewed(angle, tilt, *p), *p)).fold((f64::NEG_INFINITY, Point::ORIGIN), |a, b| if b.0 > a.0 { b } else { a });
            if worst > 1e-6 {
                eprintln!("VERIFY {angle}/{tilt} r{radius} {backs:?}: face {kind} stands {worst} outside at {at:?}");
            }
            if let Some(SurfaceGeometry::BSpline(s)) = m.geometry().surface(data.surface) {
                // Normals along the mesh: any triangle whose surface normal at
                // its vertices points against the next.
                let mut worst_cos = 1.0_f64;
                let mut where_ = (0.0, 0.0);
                for t in &tri.triangles {
                    let ns: Vec<Option<Vector>> = t.iter().map(|&i| {
                        let (u, v) = tri.parameters[i as usize];
                        let (du, dv) = s.d1_at(u, v, T).unwrap();
                        let c = du.cross(dv);
                        (c.magnitude() > 1e-9).then(|| c / c.magnitude())
                    }).collect();
                    for (a, b) in [(0, 1), (1, 2), (0, 2)] {
                        if let (Some(x), Some(y)) = (ns[a], ns[b]) {
                            if x.dot(y) < worst_cos {
                                worst_cos = x.dot(y);
                                where_ = tri.parameters[t[a] as usize];
                            }
                        }
                    }
                }
                eprintln!("VERIFY {angle}/{tilt} r{radius} {backs:?}: bspline face domain {:?} worst adjacent normal cos {worst_cos} at {where_:?}", s.domain());
            }
        }
    }
}

#[test]
#[ignore = "probe"]
fn probe_setback() {
    let mut log = Log::default();
    let mut k = 0;
    let cases: &[(f64, [f64; 3])] = &[
        (2.0, [5.0, 5.0, 5.0]),
        (2.0, [2.5, 6.0, 4.0]),
        (1.0, [3.0, 3.0, 3.0]),
        (3.0, [8.0, 6.0, 7.0]),
        (0.5, [1.0, 4.0, 2.0]),
    ];
    let mut m = Model::new();
    for (cname, prism, solid, vertex, edges) in corners(&mut m) {
        for &(radius, backs) in cases {
            if !wanted(&mut k) {
                continue;
            }
            let name = format!("setback {cname} r{radius} {backs:?}");
            let setbacks = Setbacks {
                per_edge: (0..3).map(|k| (edges[k].clone(), backs[k])).collect(),
                default: 5.0,
                continuity: Continuity::G1,
            };
            let at = point_of(&m, &vertex);
            match ogeom_fillet::round_vertex_setback(&mut m, &solid, &vertex, radius, &setbacks, T) {
                Ok(r) => {
                    let reference = ogeom_fillet::fillet_edges(&mut m, &solid, &edges, radius, T).ok();
                    let reach = backs.iter().copied().fold(0.0, f64::max) + 3.0 * radius;
                    let v = audit(&mut log, &name, &mut m, Some(prism), &r.shape, reference.as_ref().map(|b| &b.shape), at, reach);
                    eprintln!("OK {name}: volume {v}");
                }
                Err(e) => eprintln!("REFUSED {name}: {e}"),
            }
        }
    }
    eprintln!("FINDINGS {}", log.0.len());
}

#[test]
#[ignore = "probe"]
fn probe_bulge_cases() {
    // (angle, tilt, sharp or 3 for none, radius, setbacks)
    let cases: &[(f64, f64, usize, f64, [f64; 3])] = &[
        (90.0, 0.0, 2, 2.0, [5.0, 5.0, 5.0]),
        (60.0, 0.0, 2, 2.0, [5.0, 5.0, 5.0]),
        (90.0, 0.0, 0, 1.5, [4.0, 5.0, 3.5]),
        (90.0, 3.0, 3, 0.5, [1.0, 4.0, 2.0]),
        (120.0, 0.0, 3, 3.0, [8.0, 6.0, 7.0]),
        (90.0, 0.0, 3, 2.0, [5.0, 5.0, 5.0]),
    ];
    let which: Option<usize> = std::env::var("PROBE_CASE").ok().map(|v| v.parse().unwrap());
    for (i, &(angle, tilt, sharp, radius, backs)) in cases.iter().enumerate() {
        if which.is_some_and(|w| w != i) {
            continue;
        }
        let mut m = Model::new();
        let (solid, vertex, edges) = skewed(&mut m, angle, tilt);
        let setbacks = Setbacks {
            per_edge: (0..3).map(|k| (edges[k].clone(), backs[k])).collect(),
            default: 5.0,
            continuity: Continuity::G1,
        };
        let r = if sharp == 3 {
            ogeom_fillet::round_vertex_setback(&mut m, &solid, &vertex, radius, &setbacks, T)
        } else {
            ogeom_fillet::round_vertex_setback_sharp(&mut m, &solid, &vertex, &edges[sharp], radius, &setbacks, T)
        };
        eprintln!("CASE {i} {angle}/{tilt} s{sharp} r{radius} {backs:?}: {:?}", r.err().map(|e| e.to_string()));
    }
}
