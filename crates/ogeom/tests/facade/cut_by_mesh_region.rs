//! A solid cut by a break surface marked on a mesh: the break extended past
//! the solid, thickened into a slab and cut away, one side kept. The kept
//! piece's volume against a closed form or a quadrature, and the solid's
//! faces the break does not reach on their own surfaces.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{Built, check, make_box, make_cylinder, volume_properties};
use ogeom::boolean::cut_by_mesh_region;
use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::{Frame, Point};
use ogeom::mesh::{Deflection, Side};
use ogeom::topo::{Model, Shape, ShapeType, Triangulation, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A grid of `nu` by `nv` quads over `[0, 1]²` mapped by `at`, each split
/// in two, wound counter-clockwise about `∂u × ∂v`.
fn grid(nu: u32, nv: u32, at: impl Fn(f64, f64) -> Point) -> Triangulation {
    let mut mesh = Triangulation::new();
    for j in 0..=nv {
        for i in 0..=nu {
            mesh.positions.push(at(
                f64::from(i) / f64::from(nu),
                f64::from(j) / f64::from(nv),
            ));
        }
    }
    let id = |i: u32, j: u32| j * (nu + 1) + i;
    for j in 0..nv {
        for i in 0..nu {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
    }
    mesh
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

fn surface_of(model: &Model, face: &Shape) -> SurfaceGeometry {
    let data = model.node(face).unwrap().data().as_face().unwrap().clone();
    model.geometry().surface(data.surface).unwrap().clone()
}

fn assert_valid(model: &Model, shape: &Shape) {
    let diagnosis = check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

fn area(model: &Model, face: &Shape) -> f64 {
    ogeom::algo::surface_properties(model, face, Deflection::default(), T)
        .unwrap()
        .mass
}

/// Each face of `solid` the result holds a part of lies on the face's own
/// surface, placed as it was, and the faces `untouched` names by index come
/// through whole: one face each, of the area they had.
fn faces_keep_their_surfaces(model: &Model, solid: &Shape, built: &Built, untouched: &[usize]) {
    let result = explore_unique(model, &built.shape, ShapeType::Face).unwrap();
    for (k, face) in explore_unique(model, solid, ShapeType::Face)
        .unwrap()
        .iter()
        .enumerate()
    {
        let surface = surface_of(model, face);
        let location = model
            .node(face)
            .unwrap()
            .data()
            .as_face()
            .unwrap()
            .location
            .clone();
        let images: Vec<Shape> = built
            .history
            .copy_of(face)
            .map_or_else(
                || built.history.modified(face).to_vec(),
                |copy| vec![copy.clone()],
            )
            .into_iter()
            .filter(|f| result.iter().any(|g| g.is_same(f)))
            .collect();
        if untouched.contains(&k) {
            assert_eq!(images.len(), 1, "face {k}");
            let (was, is) = (area(model, face), area(model, &images[0]));
            assert!((was - is).abs() < 1e-9 * was, "face {k}: {was} became {is}");
        }
        for image in &images {
            assert_eq!(surface_of(model, image), surface, "face {k}");
            assert_eq!(
                &model
                    .node(image)
                    .unwrap()
                    .data()
                    .as_face()
                    .unwrap()
                    .location,
                &location,
                "face {k}"
            );
            assert_eq!(image.location(), face.location(), "face {k}");
        }
    }
}

fn block(model: &mut Model) -> Shape {
    make_box(model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape
}

/// A square patch of the plane `z = 4` inside a block of 10, facing up,
/// extended to the block's sides: the block above `z = 4.25` or below
/// `z = 3.75`, each side losing the offset of 0.25.
#[test]
fn a_block_cut_by_a_flat_patch_inside_it_is_the_block_cut_by_its_plane() {
    let patch = grid(4, 4, |u, v| Point::new(3.0 + 4.0 * u, 3.0 + 4.0 * v, 4.0));
    let region = vec![true; patch.triangles.len()];
    for (keep, low, high) in [(Side::Above, 4.25, 10.0), (Side::Below, 0.0, 3.75)] {
        let mut model = Model::new();
        let solid = block(&mut model);
        let built = cut_by_mesh_region(&mut model, &solid, &patch, &region, 0.25, keep, T).unwrap();
        assert_valid(&model, &built.shape);
        assert_eq!(model.node(&built.shape).unwrap().kind(), ShapeType::Solid);
        let v = volume(&model, &built.shape);
        let expected = 100.0 * (high - low);
        assert!(
            (v - expected).abs() < 1e-9 * expected,
            "{keep:?}: {v} against {expected}"
        );
        let bounds = ogeom::algo::tight_bounds(&model, &built.shape, T).unwrap();
        assert!((bounds.low().unwrap().z - low).abs() < T.confusion());
        assert!((bounds.high().unwrap().z - high).abs() < T.confusion());
        // Six planes, as the block cut by the plane has.
        let faces = explore_unique(&model, &built.shape, ShapeType::Face).unwrap();
        assert_eq!(faces.len(), 6, "{keep:?}");
        assert!(
            faces
                .iter()
                .all(|f| matches!(surface_of(&model, f), SurfaceGeometry::Plane(_)))
        );
        // The far cap comes through whole.
        let far = explore_unique(&model, &solid, ShapeType::Face)
            .unwrap()
            .iter()
            .position(|f| {
                let b = ogeom::algo::tight_bounds(&model, f, T).unwrap();
                let z = if keep == Side::Above { 10.0 } else { 0.0 };
                (b.low().unwrap().z - z).abs() < T.confusion()
                    && (b.high().unwrap().z - z).abs() < T.confusion()
            })
            .unwrap();
        faces_keep_their_surfaces(&model, &solid, &built, &[far]);
    }
}

/// A patch of a cylinder of radius 12 about the `x` axis, through
/// `(0, 0, 5)`, crossing a solid cylinder of radius 5 and height 10 stood
/// on the `z` axis: the part above the patch offset outward by 0.2, its
/// volume against a quadrature of the height left over the disk.
#[test]
fn a_cylinder_cut_by_a_curved_patch_keeps_what_lies_above_it() {
    let (radius, height, bend, offset) = (5.0, 10.0, 12.0, 0.2);
    let (centre, spread) = (5.0 - bend, 0.6_f64);
    let patch = grid(32, 96, |u, v| {
        let angle = spread * (2.0 * v - 1.0);
        Point::new(
            8.0 * (1.0 - 2.0 * u),
            -bend * angle.sin(),
            centre + bend * angle.cos(),
        )
    });
    // `u` runs along -x and `v` along -y at the top, so the patch faces up.
    let region = vec![true; patch.triangles.len()];
    let mut model = Model::new();
    let solid = make_cylinder(&mut model, Frame::WORLD, radius, height, T)
        .unwrap()
        .shape;
    let built =
        cut_by_mesh_region(&mut model, &solid, &patch, &region, offset, Side::Above, T).unwrap();
    assert_valid(&model, &built.shape);
    assert_eq!(model.node(&built.shape).unwrap().kind(), ShapeType::Solid);

    // Over the disk at y = radius sin φ, a chord 2 radius cos φ long,
    // under the cap and over the offset cylinder.
    let outer = bend + offset;
    let steps = 2000;
    let h = core::f64::consts::PI / f64::from(steps);
    let f = |phi: f64| {
        let y = radius * phi.sin();
        let chord = 2.0 * radius * phi.cos();
        chord * (height - centre - (outer * outer - y * y).sqrt()) * radius * phi.cos()
    };
    let mut expected = f(-core::f64::consts::FRAC_PI_2) + f(core::f64::consts::FRAC_PI_2);
    for k in 1..steps {
        let phi = f64::from(k).mul_add(h, -core::f64::consts::FRAC_PI_2);
        expected += if k % 2 == 1 { 4.0 } else { 2.0 } * f(phi);
    }
    expected *= h / 3.0;
    let v = volume(&model, &built.shape);
    // The patch's chords stand off the cylinder by at most
    // bend (1 - cos(step / 2)) over the disk's area.
    let step = 2.0 * spread / 96.0;
    let sag = bend * (1.0 - (step / 2.0).cos());
    let disk = core::f64::consts::PI * radius * radius;
    assert!(
        (v - expected).abs() <= sag * disk,
        "{v} against {expected}, apart by {} where the facets allow {}",
        (v - expected).abs(),
        sag * disk
    );
    // The top cap comes through whole, and the side is cut on its own
    // surface.
    let faces = explore_unique(&model, &solid, ShapeType::Face).unwrap();
    let top = faces
        .iter()
        .position(|f| {
            let b = ogeom::algo::tight_bounds(&model, f, T).unwrap();
            (b.low().unwrap().z - height).abs() < T.confusion()
        })
        .unwrap();
    faces_keep_their_surfaces(&model, &solid, &built, &[top]);
}

#[test]
fn an_empty_region_is_refused() {
    let patch = grid(2, 2, |u, v| Point::new(10.0 * u, 10.0 * v, 4.0));
    let mut model = Model::new();
    let solid = block(&mut model);
    let region = vec![false; patch.triangles.len()];
    let error = cut_by_mesh_region(&mut model, &solid, &patch, &region, 0.1, Side::Above, T)
        .unwrap_err()
        .to_string();
    assert!(error.contains("the region is empty"), "{error}");
}

/// A patch above the block: the slab passes over it, and the block lies
/// whole on one side.
#[test]
fn a_region_that_does_not_separate_the_shape_is_refused() {
    let patch = grid(2, 2, |u, v| Point::new(3.0 + 4.0 * u, 3.0 + 4.0 * v, 14.0));
    let region = vec![true; patch.triangles.len()];
    for keep in [Side::Above, Side::Below] {
        let mut model = Model::new();
        let solid = block(&mut model);
        let error = cut_by_mesh_region(&mut model, &solid, &patch, &region, 0.1, keep, T)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("the region does not separate the shape"),
            "{keep:?}: {error}"
        );
    }
}

/// A region folded over itself has no height over any plane to extend.
#[test]
fn a_region_folding_over_itself_is_refused() {
    // Three quarters of a tube about the x axis.
    let patch = grid(4, 24, |u, v| {
        let angle = 1.5 * core::f64::consts::PI * v;
        Point::new(10.0 * u, 5.0 + 3.0 * angle.cos(), 5.0 + 3.0 * angle.sin())
    });
    let region = vec![true; patch.triangles.len()];
    let mut model = Model::new();
    let solid = block(&mut model);
    let error = cut_by_mesh_region(&mut model, &solid, &patch, &region, 0.1, Side::Above, T)
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("faces away from its average normal"),
        "{error}"
    );
}

/// A region with a hole in it is closed over the hole: the cut is the
/// same as by the whole patch.
#[test]
fn a_hole_in_the_region_is_closed() {
    let patch = grid(6, 6, |u, v| Point::new(2.0 + 6.0 * u, 2.0 + 6.0 * v, 6.0));
    // The middle two by two quads left out.
    let region: Vec<bool> = (0..patch.triangles.len())
        .map(|t| {
            let quad = t / 2;
            let (i, j) = (quad % 6, quad / 6);
            !((2..4).contains(&i) && (2..4).contains(&j))
        })
        .collect();
    let mut model = Model::new();
    let solid = block(&mut model);
    let built =
        cut_by_mesh_region(&mut model, &solid, &patch, &region, 0.5, Side::Below, T).unwrap();
    assert_valid(&model, &built.shape);
    let v = volume(&model, &built.shape);
    assert!((v - 550.0).abs() < 1e-9 * 550.0, "{v}");
}

/// A synthetic fracture: a wavy surface over a square of 20, with noise
/// of up to 0.05 at each vertex, and a block standing on the square's
/// middle and through it, cut with an offset of 0.005 and its top kept.
/// The piece holds the volume over the fracture less the offset times
/// the fracture's area, and its mating face lies within 0.01 of it.
#[test]
fn a_block_through_a_synthetic_fracture_gives_the_piece_above_it() {
    let n = 40;
    let mut seed = 0x9e37_79b9_u32;
    let mut noise = || {
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        0.05 * (f64::from(seed) / f64::from(u32::MAX) * 2.0 - 1.0)
    };
    let mut fracture = grid(n, n, |u, v| {
        let (x, y) = (20.0 * u, 20.0 * v);
        Point::new(x, y, 0.3 * (1.3 * x).sin().mul_add((0.9 * y).cos(), 5.0))
    });
    for p in &mut fracture.positions {
        p.z += noise();
    }
    let region = vec![true; fracture.triangles.len()];
    let offset = 0.005;
    let mut model = Model::new();
    let at = Frame::new(
        Point::new(1.0, 1.0, 0.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let solid = make_box(&mut model, at, (18.0, 18.0, 10.0), T)
        .unwrap()
        .shape;
    let built = cut_by_mesh_region(
        &mut model,
        &solid,
        &fracture,
        &region,
        offset,
        Side::Above,
        T,
    )
    .unwrap();
    assert_valid(&model, &built.shape);

    // Over the block's square, which the grid's lines bound: the prisms
    // between each triangle and the top, and the triangles' areas.
    let (mut under, mut surface) = (0.0, 0.0);
    for &[a, b, c] in &fracture.triangles {
        let [pa, pb, pc] = [a, b, c].map(|k| fracture.positions[k as usize]);
        let inside = [pa, pb, pc]
            .iter()
            .all(|p| (1.0..=19.0).contains(&p.x) && (1.0..=19.0).contains(&p.y));
        if !inside {
            continue;
        }
        let normal = (pb - pa).cross(pc - pa);
        under += normal.z / 2.0 * (10.0 - (pa.z + pb.z + pc.z) / 3.0);
        surface += normal.magnitude() / 2.0;
    }
    let expected = offset.mul_add(-surface, under);
    let v = volume(&model, &built.shape);
    assert!(
        (v - expected).abs() < 5e-3 * expected,
        "{v} against {expected}"
    );

    // The mating face's vertices: those below the top.
    let near = ogeom::mesh::Bvh::build(&fracture).unwrap();
    let mut mating = 0;
    for vertex in explore_unique(&model, &built.shape, ShapeType::Vertex).unwrap() {
        let data = model
            .node(&vertex)
            .unwrap()
            .data()
            .as_vertex()
            .unwrap()
            .point;
        let p = vertex.transform(model.datums()).unwrap().apply(data);
        if p.z > 9.0 {
            continue;
        }
        let d = near.nearest(p, f64::INFINITY).unwrap().distance;
        assert!(
            d <= 0.01,
            "a vertex of the mating face stands {d} off the fracture"
        );
        mating += 1;
    }
    assert!(mating > 0);
}
