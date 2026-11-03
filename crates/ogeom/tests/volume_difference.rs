//! The volume between a scan and its solid: a boss, a pocket and an offset
//! face against their closed forms, holes in the scan, noise under the
//! grid's resolution, and the error as the cell shrinks.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::time::Duration;

use ogeom::core::FastMap;
use ogeom::core::clock::Instant;

use ogeom::algo::{Form, Material, Region, make_box, make_cylinder, volume_difference};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::{Deflection, triangulate};
use ogeom::topo::{Model, Shape, Triangulation};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn block(model: &mut Model, height: f64) -> Shape {
    make_box(model, Frame::WORLD, (10.0, 10.0, height), T)
        .unwrap()
        .shape
}

fn mesh_of(model: &Model, shape: &Shape) -> Triangulation {
    let deflection = Deflection {
        chord: 1e-3,
        ..Deflection::default()
    };
    triangulate(model, shape, deflection, T).unwrap()
}

/// The volume a closed mesh bounds, by the divergence theorem.
fn enclosed(mesh: &Triangulation) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            let cross = [
                b.y.mul_add(c.z, -(b.z * c.y)),
                b.z.mul_add(c.x, -(b.x * c.z)),
                b.x.mul_add(c.y, -(b.y * c.x)),
            ];
            a.z.mul_add(cross[2], a.x.mul_add(cross[0], a.y * cross[1])) / 6.0
        })
        .sum()
}

/// A block 10 wide and 10 high, and a scan of it with a cylinder of
/// radius 2 and height 3 standing on its top, at its centre.
fn bossed() -> (Model, Shape, Triangulation) {
    let mut model = Model::new();
    let solid = block(&mut model, 10.0);
    let boss = make_cylinder(&mut model, at(5.0, 5.0, 10.0), 2.0, 3.0, T)
        .unwrap()
        .shape;
    let fused = ogeom::boolean::fuse(&mut model, &solid, &boss, T)
        .unwrap()
        .shape;
    let scan = mesh_of(&model, &fused);
    (model, solid, scan)
}

fn only(regions: &[Region]) -> &Region {
    assert_eq!(regions.len(), 1, "{regions:#?}");
    &regions[0]
}

/// A boss the scan has and the solid lacks is one excess lump, as deep as
/// the boss is tall, of the boss's volume within the grid's error. The
/// error, averaged over seven cells about each width so that no one
/// alignment of the grid with the boss decides it, about halves with the
/// cell.
#[test]
fn a_boss_on_the_scan_is_one_excess_lump_its_height_deep() {
    let (model, solid, scan) = bossed();
    let boss = enclosed(&scan) - 1000.0;
    // The cells the boss's boundary passes through: its area, top, bottom
    // and side, times √3 cells.
    let area = 2.0 * std::f64::consts::PI * 2.0 * (2.0 + 3.0);
    let mut means = Vec::new();
    for base in [0.4, 0.2, 0.1, 0.05] {
        let mut sum = 0.0;
        for k in 0..7 {
            let cell = base * (1.0 + f64::from(k) / 7.0);
            let regions = volume_difference(&scan, &model, &solid, cell, T).unwrap();
            let r = only(&regions);
            assert_eq!(r.side, Material::Excess);
            assert_eq!(r.kind, Form::Lump);
            assert!(!r.unsure, "{r:?}");
            assert!((r.depth - 3.0).abs() < 1e-9, "depth {}", r.depth);
            let error = (r.volume - boss).abs();
            assert!(error <= area * 3.0f64.sqrt() * cell, "cell {cell}: {error}");
            sum += error;
        }
        means.push((base, sum / 7.0));
    }
    eprintln!("boss of {boss}: (cell, mean error) {means:?}");
    for pair in means.windows(2) {
        assert!(pair[1].1 < 0.7 * pair[0].1, "{means:?}");
    }
    assert!(
        means[3].1 / boss < 0.05,
        "within 5 % at a fortieth of its width"
    );
}

/// A pocket the scan has and the solid lacks is one missing region, as
/// deep as the pocket.
#[test]
fn a_pocket_in_the_scan_is_one_missing_lump_its_depth_deep() {
    let mut model = Model::new();
    let solid = block(&mut model, 10.0);
    let tool = make_cylinder(&mut model, at(5.0, 5.0, 7.0), 1.5, 4.0, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &solid, &tool, T)
        .unwrap()
        .shape;
    let scan = mesh_of(&model, &cut);
    let pocket = 1000.0 - enclosed(&scan);
    let regions = volume_difference(&scan, &model, &solid, 0.1, T).unwrap();
    let r = only(&regions);
    assert_eq!(r.side, Material::Missing);
    assert_eq!(r.kind, Form::Lump);
    assert!(!r.unsure);
    assert!((r.depth - 3.0).abs() < 1e-9, "depth {}", r.depth);
    assert!(
        (r.volume - pocket).abs() / pocket < 0.05,
        "{} of {pocket}",
        r.volume
    );
}

/// The solid's boss missing from the scan is one missing lump, read from
/// the solid's side.
#[test]
fn a_boss_the_scan_lacks_is_one_missing_lump_its_height_deep() {
    let (mut model, _, _) = bossed();
    let plain = block(&mut model, 10.0);
    let scan = mesh_of(&model, &plain);
    let boss = make_cylinder(&mut model, at(5.0, 5.0, 10.0), 2.0, 3.0, T)
        .unwrap()
        .shape;
    let solid = ogeom::boolean::fuse(&mut model, &plain, &boss, T)
        .unwrap()
        .shape;
    let regions = volume_difference(&scan, &model, &solid, 0.1, T).unwrap();
    let r = only(&regions);
    assert_eq!(r.side, Material::Missing);
    assert_eq!(r.kind, Form::Lump);
    assert!((r.depth - 3.0).abs() < 1e-3, "depth {}", r.depth);
    let exact = std::f64::consts::PI * 4.0 * 3.0;
    assert!(
        (r.volume - exact).abs() / exact < 0.05,
        "{} of {exact}",
        r.volume
    );
}

/// A face of the scan standing 0.3 off the solid's is one layer 0.3 deep,
/// of the slab's volume: the slab holds three whole layers of nodes.
#[test]
fn a_face_standing_off_is_a_layer_as_deep_as_the_offset() {
    let mut model = Model::new();
    let solid = block(&mut model, 10.0);
    let thicker = block(&mut model, 10.3);
    let scan = mesh_of(&model, &thicker);
    let regions = volume_difference(&scan, &model, &solid, 0.1, T).unwrap();
    let r = only(&regions);
    assert_eq!(r.side, Material::Excess);
    assert_eq!(r.kind, Form::Layer);
    assert!(!r.unsure);
    assert!((r.depth - 0.3).abs() < 1e-9, "depth {}", r.depth);
    assert!((r.volume - 30.0).abs() < 1e-9, "volume {}", r.volume);

    // At a cell wider than the offset the slab is below the resolution.
    let coarse = volume_difference(&scan, &model, &solid, 0.5, T).unwrap();
    assert!(coarse.is_empty(), "{coarse:#?}");
}

/// The scan with triangles taken out of its boss's top and of a side far
/// from it: the boss reads unsure, a hole over plain face reads nothing.
#[test]
fn a_region_by_a_hole_of_the_scan_is_unsure() {
    let (model, solid, closed) = bossed();
    let sure = only(&volume_difference(&closed, &model, &solid, 0.1, T).unwrap()).clone();

    let mut open = closed.clone();
    // The widest triangle of the boss's top.
    let area = |t: &[u32; 3]| {
        let [a, b, c] = t.map(|v| open.positions[v as usize]);
        ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)).abs()
    };
    let top = (0..open.triangles.len())
        .filter(|&t| {
            open.triangles[t]
                .iter()
                .all(|&v| open.positions[v as usize].z == 13.0)
        })
        .max_by(|&a, &b| area(&open.triangles[a]).total_cmp(&area(&open.triangles[b])))
        .unwrap();
    let side = open
        .triangles
        .iter()
        .position(|t| t.iter().all(|&v| open.positions[v as usize].x == 0.0))
        .unwrap();
    let (first, second) = (top.max(side), top.min(side));
    open.triangles.remove(first);
    open.triangles.remove(second);
    assert!(!open.is_closed());
    let regions = volume_difference(&open, &model, &solid, 0.1, T).unwrap();
    let r = only(&regions);
    assert_eq!(r.side, Material::Excess);
    assert!(r.unsure, "{r:?}");
    assert!(r.volume <= sure.volume, "{} over {}", r.volume, sure.volume);

    // Only the side's hole: nothing differs, nothing is reported.
    let mut model = model;
    let plain = block(&mut model, 10.0);
    let mut scan = mesh_of(&model, &plain);
    let side = scan
        .triangles
        .iter()
        .position(|t| t.iter().all(|&v| scan.positions[v as usize].x == 0.0))
        .unwrap();
    scan.triangles.remove(side);
    let regions = volume_difference(&scan, &model, &plain, 0.1, T).unwrap();
    assert!(regions.is_empty(), "{regions:#?}");
}

/// A scan stored as loose triangles, each with corners of its own, reads
/// as the indexed scan does.
#[test]
fn loose_triangles_read_as_the_surface_they_make() {
    let (model, solid, scan) = bossed();
    let mut loose = Triangulation::new();
    for t in &scan.triangles {
        let first = u32::try_from(loose.positions.len()).unwrap();
        loose
            .positions
            .extend(t.iter().map(|&v| scan.positions[v as usize]));
        loose.triangles.push([first, first + 1, first + 2]);
    }
    assert!(!loose.is_closed());
    let indexed = volume_difference(&scan, &model, &solid, 0.1, T).unwrap();
    let read = volume_difference(&loose, &model, &solid, 0.1, T).unwrap();
    assert_eq!(read, indexed);
}

/// A cell that is not a distance, a scan with no triangles and a shape
/// with no inside are refused by name.
#[test]
fn what_cannot_be_read_is_refused() {
    let (model, solid, scan) = bossed();
    for cell in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let e = volume_difference(&scan, &model, &solid, cell, T).unwrap_err();
        assert!(matches!(e, ogeom::core::OgeomError::Construction(_)), "{e}");
    }
    let e = volume_difference(&Triangulation::new(), &model, &solid, 0.1, T).unwrap_err();
    assert!(matches!(e, ogeom::core::OgeomError::NullObject(_)), "{e}");
    let face = ogeom::topo::explore_unique(&model, &solid, ogeom::topo::ShapeType::Face)
        .unwrap()
        .remove(0);
    let e = volume_difference(&scan, &model, &face, 0.1, T).unwrap_err();
    assert!(matches!(e, ogeom::core::OgeomError::Construction(_)), "{e}");
}

/// The closed surface of the block `[0, 10]³` on a grid of `n` by `n`
/// quads a face, each split in two, wound outward, each vertex moved by
/// `shift`.
fn gridded_block(n: u32, mut shift: impl FnMut(Point) -> Point) -> Triangulation {
    let mut mesh = Triangulation::new();
    let mut index = FastMap::default();
    let step = 10.0 / f64::from(n);
    let mut vertex = |mesh: &mut Triangulation, c: [u32; 3]| -> u32 {
        *index.entry(c).or_insert_with(|| {
            mesh.positions.push(shift(Point::new(
                f64::from(c[0]) * step,
                f64::from(c[1]) * step,
                f64::from(c[2]) * step,
            )));
            u32::try_from(mesh.positions.len() - 1).unwrap()
        })
    };
    // Each face: the fixed axis, its side, and the two axes it spans in
    // the order that winds it outward.
    for (axis, high, u, v) in [
        (0, false, 2, 1),
        (0, true, 1, 2),
        (1, false, 0, 2),
        (1, true, 2, 0),
        (2, false, 1, 0),
        (2, true, 0, 1),
    ] {
        for j in 0..n {
            for i in 0..n {
                let corner = |di: u32, dj: u32| {
                    let mut c = [0; 3];
                    c[axis] = if high { n } else { 0 };
                    c[u] = i + di;
                    c[v] = j + dj;
                    c
                };
                let a = vertex(&mut mesh, corner(0, 0));
                let b = vertex(&mut mesh, corner(1, 0));
                let c = vertex(&mut mesh, corner(1, 1));
                let d = vertex(&mut mesh, corner(0, 1));
                mesh.triangles.push([a, b, c]);
                mesh.triangles.push([a, c, d]);
            }
        }
    }
    mesh
}

/// A deterministic noise in `[-amplitude, amplitude]`.
fn noise(seed: &mut u64, amplitude: f64) -> f64 {
    *seed = seed
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    #[allow(clippy::cast_precision_loss, reason = "test noise")]
    let unit = (*seed >> 11) as f64 / (1u64 << 53) as f64;
    amplitude * 2.0f64.mul_add(unit, -1.0)
}

/// A scan of the block with noise of three quarters of a cell, which lifts
/// nodes off its faces here and there, reports nothing: each is shallower
/// than a cell. The answers are the same on one thread as on many.
#[test]
fn noise_under_a_cell_reports_nothing_and_threads_change_nothing() {
    let mut seed = 7_u64;
    let mut model = Model::new();
    let solid = block(&mut model, 10.0);
    let scan = gridded_block(40, |p| {
        let mut q = p;
        // Corners and edges stay put so the faces still meet.
        let on = |c: f64| c == 0.0 || c == 10.0;
        let shared = [p.x, p.y, p.z].iter().filter(|&&c| on(c)).count();
        if shared == 1 {
            let d = noise(&mut seed, 0.15);
            if on(p.x) {
                q.x += d;
            } else if on(p.y) {
                q.y += d;
            } else {
                q.z += d;
            }
        }
        q
    });
    assert!(scan.is_closed());
    let regions = volume_difference(&scan, &model, &solid, 0.2, T).unwrap();
    assert!(regions.is_empty(), "{regions:#?}");

    let (model, solid, scan) = bossed();
    ogeom::core::parallel::set_threads(1);
    let one = volume_difference(&scan, &model, &solid, 0.1, T).unwrap();
    ogeom::core::parallel::set_threads(8);
    let many = volume_difference(&scan, &model, &solid, 0.1, T).unwrap();
    ogeom::core::parallel::set_threads(0);
    assert_eq!(one, many);
}

/// A scan of about a million triangles, a bump of height 1 raised on its
/// top and noise over all of it, against the block.
#[test]
#[ignore = "heavy"]
fn a_million_triangle_scan_in_bounded_time() {
    let mut seed = 11_u64;
    let mut model = Model::new();
    let solid = block(&mut model, 10.0);
    let scan = gridded_block(290, |p| {
        let mut q = p;
        let r = (p.x - 5.0).hypot(p.y - 5.0);
        if p.z == 10.0 && r < 3.0 && p.x > 0.0 && p.x < 10.0 && p.y > 0.0 && p.y < 10.0 {
            q.z += 0.5 * (1.0 + (std::f64::consts::PI * r / 3.0).cos()) + noise(&mut seed, 0.01);
        }
        q
    });
    assert!(scan.triangles.len() >= 1_000_000);
    let start = Instant::now();
    let regions = volume_difference(&scan, &model, &solid, 0.05, T).unwrap();
    let took = start.elapsed();
    eprintln!("{} triangles: {took:?}, {regions:#?}", scan.triangles.len());
    let r = only(&regions);
    assert_eq!(r.side, Material::Excess);
    assert!((r.depth - 1.0).abs() < 0.02, "depth {}", r.depth);
    // The bump's volume, 2π ∫₀³ ½ (1 + cos(πr/3)) r dr.
    let bump = std::f64::consts::PI * (4.5 - 18.0 / std::f64::consts::PI.powi(2));
    assert!(
        (r.volume - bump).abs() / bump < 0.05,
        "{} of {bump}",
        r.volume
    );
    assert!(took < Duration::from_secs(40), "took {took:?}");
}

/// A sphere's scan of about a million triangles, its own tessellation,
/// against the sphere: every node near the curved face is read against the
/// true surface, and the chords' sliver is under the resolution.
#[test]
#[ignore = "heavy"]
fn a_million_triangle_scan_of_a_curved_solid_in_bounded_time() {
    let mut model = Model::new();
    let solid = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 5.0, T)
        .unwrap()
        .shape;
    let fine = Deflection {
        chord: 1e-4,
        angular: 0.01,
        ..Deflection::default()
    };
    let scan = triangulate(&model, &solid, fine, T).unwrap();
    assert!(scan.triangles.len() >= 1_000_000);
    let start = Instant::now();
    let regions = volume_difference(&scan, &model, &solid, 0.05, T).unwrap();
    let took = start.elapsed();
    eprintln!("{} triangles: {took:?}", scan.triangles.len());
    assert!(regions.is_empty(), "{regions:#?}");
    assert!(took < Duration::from_secs(60), "took {took:?}");
}
