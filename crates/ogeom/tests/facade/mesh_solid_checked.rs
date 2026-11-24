//! A mesh converted back to a solid comes back passing `check`, or is
//! refused by name: a seam shorter than the scatter its curved face is
//! fitted to, and a mesh passing through itself.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{MeshSolidOptions, check, solid_from_mesh};
use ogeom::core::Tolerances;
use ogeom::math::Point;
use ogeom::topo::{Model, Triangulation};

const T: Tolerances = Tolerances::millimetres();

/// A plate `THICK` thick along `y` whose outline in `(x, z)` is a 2 by 1
/// rectangle with its far top corner rounded, radius 0.3, in five
/// fifteen-degree steps from the side up to 75 degrees. From there one
/// step, a kink, runs to the top 0.06 short of where the arc would have
/// met it: its corners on the top stand 0.006 off the arc's cylinder, and
/// the kink's edge with the top is as long as the plate is thick.
fn kinked_plate() -> Triangulation {
    const THICK: f64 = 2.17e-3;
    let (radius, short) = (0.3, 0.0603);
    let (cx, cz) = (2.0 - radius, 1.0 - radius);
    let mut outline = vec![(0.0, 0.0), (2.0, 0.0)];
    for k in 0..=5 {
        let a = f64::from(k * 15).to_radians();
        outline.push((cx + radius * a.cos(), cz + radius * a.sin()));
    }
    outline.push((cx - short, 1.0));
    outline.push((0.0, 1.0));
    let n = u32::try_from(outline.len()).unwrap();
    let mut mesh = Triangulation::new();
    for y in [0.0, THICK] {
        for &(x, z) in &outline {
            mesh.positions.push(Point::new(x, y, z));
        }
    }
    let (mx, mz) = outline
        .iter()
        .fold((0.0, 0.0), |(sx, sz), &(x, z)| (sx + x, sz + z));
    let count = f64::from(n);
    mesh.positions.push(Point::new(mx / count, 0.0, mz / count));
    mesh.positions
        .push(Point::new(mx / count, THICK, mz / count));
    let (bottom, top) = (2 * n, 2 * n + 1);
    for i in 0..n {
        let j = (i + 1) % n;
        mesh.triangles.push([bottom, j, i]);
        mesh.triangles.push([top, n + i, n + j]);
        mesh.triangles.push([i, j, n + j]);
        mesh.triangles.push([i, n + j, n + i]);
    }
    mesh
}

/// The kink, two triangles flat to each other, is given to the fillet's
/// cylinder, which then meets the top along the kink's edge there: a
/// ruling of the cylinder, 0.006 below the top, whose image on the top
/// stands that far off an edge 0.00217 long. `check` reports the edge
/// suspect. The kink is given back, a face of its own, and the cylinder
/// stays.
#[test]
fn a_kink_beside_a_thin_fillet_keeps_its_pcurves_on_their_edges() {
    let mesh = kinked_plate();
    let mut model = Model::new();
    let options = MeshSolidOptions {
        coplanar_distance: Some(4e-4),
        ..MeshSolidOptions::default()
    };
    let solid = solid_from_mesh(&mut model, &mesh, &options, T).unwrap();
    let found = check(&model, &solid.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert_eq!(
        (solid.report.curved_faces, solid.report.curved_faceted),
        (1, 0)
    );
}

/// A box whose top is drawn as four triangles fanned to a point half its
/// height below its bottom: the fan passes through the bottom, and every
/// face built on those triangles crosses it. Built, the solid is not one
/// region, and the conversion says where.
#[test]
fn a_mesh_passing_through_itself_is_refused_where_it_crosses() {
    let mut mesh = Triangulation::new();
    for (x, y, z) in [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
        (0.5, 0.5, -0.5),
    ] {
        mesh.positions.push(Point::new(x, y, z));
    }
    mesh.triangles.extend([
        [0, 2, 1],
        [0, 3, 2],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
        [4, 5, 8],
        [5, 6, 8],
        [6, 7, 8],
        [7, 4, 8],
    ]);
    let mut model = Model::new();
    let refused = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T)
        .expect_err("a mesh passing through itself bounds no solid");
    let text = refused.to_string();
    assert!(
        text.contains("the mesh passes through itself at (") && text.contains(", 0.000000)"),
        "{text}"
    );
}
