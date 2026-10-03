//! The mesh conversion's steps taken one by one: the regions found, merged,
//! split and fitted, then built. Each corrected conversion is measured
//! against the exact solid its mesh was drawn from, and against what the
//! automatic conversion builds where it finds the same regions.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    FallbackReason, FitConstraints, MeshRegions, MeshSolid, MeshSolidOptions, RegionFallback,
    RegionId, RegionRefusal, SurfaceKind, check, solid_from_mesh, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::math::{Axis, Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, Triangulation, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The kinds of surface a shape's faces are built on, counted: planes,
/// cylinders, cones, spheres, tori and anything else.
fn kinds(model: &Model, shape: &Shape) -> [usize; 6] {
    use ogeom::geom::SurfaceGeometry as S;
    let mut out = [0; 6];
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        out[match model.geometry().surface(data.surface).unwrap() {
            S::Plane(_) => 0,
            S::Cylinder(_) => 1,
            S::Cone(_) => 2,
            S::Sphere(_) => 3,
            S::Torus(_) => 4,
            _ => 5,
        }] += 1;
    }
    out
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

fn valid(model: &Model, out: &MeshSolid) {
    assert!(out.closed, "{:?}", out.report);
    let diagnosis = check(model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

/// A cylinder of radius `r` and height `h` about the z axis, drawn with
/// `n` facets round it and `rows` rows up it, every vertex on the
/// cylinder, its ends fanned from their centres.
fn faceted_cylinder(r: f64, h: f64, n: u32, rows: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    let at = |i: u32, j: u32| {
        let a = core::f64::consts::TAU * f64::from(i % n) / f64::from(n);
        Point::new(r * a.cos(), r * a.sin(), h * f64::from(j) / f64::from(rows))
    };
    let index = |i: u32, j: u32| j * n + i % n;
    for j in 0..=rows {
        for i in 0..n {
            mesh.positions.push(at(i, j));
        }
    }
    for j in 0..rows {
        for i in 0..n {
            let (a, b, c, d) = (
                index(i, j),
                index(i + 1, j),
                index(i + 1, j + 1),
                index(i, j + 1),
            );
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
    }
    let bottom = u32::try_from(mesh.positions.len()).unwrap();
    mesh.positions.push(Point::new(0.0, 0.0, 0.0));
    mesh.positions.push(Point::new(0.0, 0.0, h));
    for i in 0..n {
        mesh.triangles.push([bottom, index(i + 1, 0), index(i, 0)]);
        mesh.triangles
            .push([bottom + 1, index(i, rows), index(i + 1, rows)]);
    }
    mesh
}

/// A prism along y, ten long, whose section is a ten by two rectangle under
/// a roof rising `rise` to a ridge at its middle. Its volume is the
/// section's area times its length, exactly.
fn ridged_prism(rise: f64) -> (Triangulation, f64) {
    let section = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 2.0),
        (5.0, 2.0 + rise),
        (0.0, 2.0),
    ];
    let mut mesh = Triangulation::new();
    for y in [0.0, 10.0] {
        for &(x, z) in &section {
            mesh.positions.push(Point::new(x, y, z));
        }
    }
    for k in 1..4 {
        mesh.triangles.push([0, k + 1, k]);
        mesh.triangles.push([5, 5 + k, 5 + k + 1]);
    }
    for k in 0..5 {
        let next = (k + 1) % 5;
        mesh.triangles.push([k, next, 5 + next]);
        mesh.triangles.push([k, 5 + next, 5 + k]);
    }
    (mesh, 10.0 * (20.0 + 5.0 * rise))
}

/// The vertex of the regions' mesh at a point.
fn vertex_at(regions: &MeshRegions, p: Point) -> u32 {
    let at = regions
        .points()
        .iter()
        .position(|q| q.distance(p) < 1e-9)
        .expect("a vertex there");
    u32::try_from(at).unwrap()
}

/// The regions whose surface is a plane facing sideways: the facets round
/// a faceted cylinder.
fn side_facets(regions: &MeshRegions) -> Vec<RegionId> {
    regions
        .regions()
        .into_iter()
        .filter(|r| match &r.surface {
            Some(ogeom::algo::Canonical::Plane(p)) => p.frame().z().vector().z.abs() < 0.5,
            _ => false,
        })
        .map(|r| r.id)
        .collect()
}

/// Merge every side facet of a faceted cylinder into one region, each
/// merge taking a neighbour of what is merged so far.
fn merge_sides(regions: &mut MeshRegions) -> RegionId {
    let sides = side_facets(regions);
    merge_all(regions, sides)
}

/// Merge the regions into one, each merge taking a neighbour of what is
/// merged so far.
fn merge_all(regions: &mut MeshRegions, mut ids: Vec<RegionId>) -> RegionId {
    let mut merged = ids.remove(0);
    while !ids.is_empty() {
        let beside = regions.region(merged).unwrap().neighbours;
        let k = ids
            .iter()
            .position(|s| beside.contains(s))
            .expect("a region beside the merged ones");
        merged = regions.merge(merged, ids.remove(k)).unwrap();
    }
    merged
}

/// A mesh from its vertices and the faces round them, each face a convex
/// polygon fanned from its first vertex.
fn from_polygons(points: Vec<Point>, polygons: &[Vec<u32>]) -> Triangulation {
    let mut mesh = Triangulation::new();
    mesh.positions = points;
    for polygon in polygons {
        for k in 1..polygon.len() - 1 {
            mesh.triangles
                .push([polygon[0], polygon[k], polygon[k + 1]]);
        }
    }
    mesh
}

/// The vertex `i` of `n` round the rim at height `z` of a faceted cylinder
/// of radius `r`.
fn rim_vertex(regions: &MeshRegions, r: f64, n: u32, i: u32, z: f64) -> u32 {
    let a = core::f64::consts::TAU * f64::from(i) / f64::from(n);
    vertex_at(regions, Point::new(r * a.cos(), r * a.sin(), z))
}

/// The path that cuts the facets from angle step `from` to `to` off a
/// one-row band: up the first ruling, along the top rim, down the last.
fn across_the_band(regions: &MeshRegions, r: f64, h: f64, n: u32, from: u32, to: u32) -> Vec<u32> {
    let mut path = vec![rim_vertex(regions, r, n, from, 0.0)];
    path.extend((from..=to).map(|i| rim_vertex(regions, r, n, i, h)));
    path.push(rim_vertex(regions, r, n, to, 0.0));
    path
}

/// Found and built with no step between, the regions build exactly what
/// the automatic conversion does: the same faces on the same surfaces, the
/// same volume to the last bit, and the same report.
#[test]
fn finding_then_building_is_the_automatic_conversion() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, 4.0, 12.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &drilled, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &drilled, &edges[..4], 1.5, T)
        .map_or(drilled, |r| r.shape);
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 7.0, T)
        .unwrap()
        .shape;
    let mut meshes = vec![faceted_cylinder(10.0, 20.0, 10, 2), ridged_prism(0.02).0];
    for shape in [&rounded, &ball] {
        meshes.push(
            ogeom::mesh::triangulate(&model, shape, Deflection::with_chord(0.05).unwrap(), T)
                .unwrap(),
        );
    }
    for mesh in &meshes {
        let options = MeshSolidOptions::default();
        let mut auto = Model::new();
        let direct = solid_from_mesh(&mut auto, mesh, &options, T).unwrap();
        let mut stepped = Model::new();
        let built = MeshRegions::find(mesh, &options, T)
            .unwrap()
            .build(&mut stepped)
            .unwrap();
        assert_eq!(built.report, direct.report);
        assert_eq!(
            built.coplanar_distance.to_bits(),
            direct.coplanar_distance.to_bits()
        );
        assert_eq!(kinds(&stepped, &built.shape), kinds(&auto, &direct.shape));
        assert_eq!(
            volume(&stepped, &built.shape).to_bits(),
            volume(&auto, &direct.shape).to_bits()
        );
    }
}

/// A cylinder drawn with ten facets round it turns 36 degrees from facet
/// to facet, past the crease angle: the automatic conversion builds it as
/// ten flat faces, and its volume is the decagonal prism's. The facets
/// merged into one region are recognized as the cylinder their vertices
/// lie on, and build what the automatic conversion builds where the crease
/// angle lets it see the cylinder: three faces, the cylinder's volume.
#[test]
fn merged_facets_build_the_cylinder_they_lie_on() {
    let (r, h) = (10.0, 20.0);
    let exact = core::f64::consts::PI * r * r * h;
    let mesh = faceted_cylinder(r, h, 10, 2);
    let options = MeshSolidOptions::default();

    let mut auto = Model::new();
    let faceted = solid_from_mesh(&mut auto, &mesh, &options, T).unwrap();
    valid(&auto, &faceted);
    assert_eq!(kinds(&auto, &faceted.shape), [12, 0, 0, 0, 0, 0]);
    let decagon = 5.0 * r * r * (core::f64::consts::TAU / 10.0).sin() * h;
    assert!((volume(&auto, &faceted.shape) - decagon).abs() < 1e-9 * decagon);

    let mut regions = MeshRegions::find(&mesh, &options, T).unwrap();
    assert_eq!(side_facets(&regions).len(), 10);
    let merged = merge_sides(&mut regions);
    let region = regions.region(merged).unwrap();
    let Some(ogeom::algo::Canonical::Cylinder(c)) = region.surface else {
        panic!("a cylinder: {region:?}");
    };
    assert!((c.radius() - r).abs() < 1e-9, "{c:?}");
    assert!(region.deviation.unwrap() <= regions.distance());
    assert_eq!(regions.regions().len(), 3);
    let mut back = Model::new();
    let out = regions.build(&mut back).unwrap();
    valid(&back, &out);
    assert_eq!(kinds(&back, &out.shape), [2, 1, 0, 0, 0, 0]);
    let corrected = volume(&back, &out.shape);
    assert!(
        (corrected - exact).abs() < 1e-6 * exact,
        "{corrected} against {exact}"
    );

    // Allowed the 36 degree turn, the automatic conversion finds the same.
    let lenient = MeshSolidOptions {
        crease: 0.7,
        ..options
    };
    let mut seen = Model::new();
    let automatic = solid_from_mesh(&mut seen, &mesh, &lenient, T).unwrap();
    valid(&seen, &automatic);
    assert_eq!(kinds(&seen, &automatic.shape), kinds(&back, &out.shape));
    let found = volume(&seen, &automatic.shape);
    assert!(
        (found - corrected).abs() < 1e-9 * exact,
        "{found} found, {corrected} corrected"
    );
}

/// A roof rising two hundredths over a half width of five is two planes
/// meeting at a ridge, leaning 0.004 apart. Found with a coplanar angle of
/// a hundredth and a distance of five hundredths, the two halves are one
/// plane region. Split along the ridge, each half is its own plane, and
/// the prism builds with its exact volume, as the automatic conversion
/// with the default distance builds it.
#[test]
fn a_region_split_along_a_ridge_becomes_two_planes() {
    let rise = 0.02;
    let (mesh, exact) = ridged_prism(rise);
    let loose = MeshSolidOptions {
        coplanar_angle: 0.01,
        coplanar_distance: Some(0.05),
        ..MeshSolidOptions::default()
    };
    let mut regions = MeshRegions::find(&mesh, &loose, T).unwrap();
    let ridge = [
        vertex_at(&regions, Point::new(5.0, 0.0, 2.0 + rise)),
        vertex_at(&regions, Point::new(5.0, 10.0, 2.0 + rise)),
    ];
    let roof = regions
        .regions()
        .into_iter()
        .find(|r| {
            r.triangles
                .iter()
                .any(|&t| regions.triangles()[t].contains(&ridge[0]))
                && r.triangles.len() == 4
        })
        .expect("the roof as one region");
    assert_eq!(regions.regions().len(), 6);

    let (left, right) = regions.split(roof.id, &ridge).unwrap();
    for id in [left, right] {
        let half = regions.region(id).unwrap();
        assert_eq!(half.triangles.len(), 2);
        assert!(matches!(
            half.surface,
            Some(ogeom::algo::Canonical::Plane(_))
        ));
        assert!(half.deviation.unwrap() < 1e-12, "{half:?}");
    }
    let mut back = Model::new();
    let out = regions.build(&mut back).unwrap();
    valid(&back, &out);
    assert_eq!(kinds(&back, &out.shape), [7, 0, 0, 0, 0, 0]);
    let split = volume(&back, &out.shape);
    assert!(
        (split - exact).abs() < 1e-9 * exact,
        "{split} against {exact}"
    );

    let mut auto = Model::new();
    let automatic = solid_from_mesh(&mut auto, &mesh, &MeshSolidOptions::default(), T).unwrap();
    valid(&auto, &automatic);
    assert_eq!(automatic.report.faces, out.report.faces);
    let found = volume(&auto, &automatic.shape);
    assert!(
        (found - split).abs() < 1e-9 * exact,
        "{found} against {split}"
    );
}

/// A cylinder fitted with its radius fixed, its axis fixed, or both, to
/// the facets of a faceted cylinder: each fit verifies at the vertices,
/// the solid builds with the cylinder's volume, and a fit about the true
/// axis lands on the true radius.
#[test]
fn a_cylinder_fitted_with_a_fixed_radius_or_axis() {
    let (r, h) = (10.0, 20.0);
    let exact = core::f64::consts::PI * r * r * h;
    let mesh = faceted_cylinder(r, h, 10, 2);
    let options = MeshSolidOptions::default();
    let base = MeshRegions::find(&mesh, &options, T).unwrap();
    let z = Axis {
        location: Point::new(0.0, 0.0, 3.0),
        direction: Direction::Z,
    };
    for constraints in [
        FitConstraints {
            radius: Some(r),
            ..FitConstraints::default()
        },
        FitConstraints {
            axis: Some(z),
            ..FitConstraints::default()
        },
        FitConstraints {
            axis: Some(z),
            radius: Some(r),
        },
    ] {
        let mut regions = base.clone();
        let sides = merge_sides(&mut regions);
        let deviation = regions
            .fit(sides, SurfaceKind::Cylinder, &constraints)
            .unwrap();
        assert!(deviation < 1e-9, "{constraints:?}: {deviation}");
        let Some(ogeom::algo::Canonical::Cylinder(c)) = regions.region(sides).unwrap().surface
        else {
            panic!("a cylinder");
        };
        assert!((c.radius() - r).abs() < 1e-9, "{constraints:?}: {c:?}");
        let axis = c.frame().z().vector();
        assert!(axis.cross(Vector::Z).magnitude() < 1e-9, "{axis:?}");
        let mut back = Model::new();
        let out = regions.build(&mut back).unwrap();
        valid(&back, &out);
        assert_eq!(kinds(&back, &out.shape), [2, 1, 0, 0, 0, 0]);
        let v = volume(&back, &out.shape);
        assert!((v - exact).abs() < 1e-6 * exact, "{constraints:?}: {v}");
    }
}

/// The curved region of a meshed primitive, refitted holding the
/// constraints, lands on the primitive's own surface, and the solid builds
/// with the primitive's volume.
fn refits_onto(
    model: &mut Model,
    shape: &Shape,
    kind: SurfaceKind,
    constraints: FitConstraints,
    lands: impl Fn(&ogeom::algo::Canonical) -> bool,
) {
    let mesh =
        ogeom::mesh::triangulate(model, shape, Deflection::with_chord(0.01).unwrap(), T).unwrap();
    let mut regions = MeshRegions::find(&mesh, &MeshSolidOptions::default(), T).unwrap();
    let curved = regions
        .regions()
        .into_iter()
        .find(|r| !matches!(r.surface, Some(ogeom::algo::Canonical::Plane(_))))
        .expect("a curved region");
    let deviation = regions.fit(curved.id, kind, &constraints).unwrap();
    assert!(deviation <= regions.distance(), "{kind:?}: {deviation}");
    let surface = regions.region(curved.id).unwrap().surface.unwrap();
    assert!(lands(&surface), "{kind:?}: {surface:?}");
    let mut back = Model::new();
    let out = regions.build(&mut back).unwrap();
    valid(&back, &out);
    let (want, got) = (volume(model, shape), volume(&back, &out.shape));
    assert!(
        (got - want).abs() < 1e-6 * want,
        "{kind:?}: {got} against {want}"
    );
}

/// A sphere with its radius fixed, a cone about its axis, and a torus
/// about its axis or with its tube's radius fixed: each fit lands on the
/// primitive the mesh was drawn from.
#[test]
fn spheres_cones_and_tori_fitted_with_what_they_fix() {
    let mut model = Model::new();
    let near = |a: f64, b: f64| (a - b).abs() < 1e-6;
    let z = Axis {
        location: Point::new(0.0, 0.0, 0.0),
        direction: Direction::Z,
    };
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 7.0, T)
        .unwrap()
        .shape;
    refits_onto(
        &mut model,
        &ball,
        SurfaceKind::Sphere,
        FitConstraints {
            radius: Some(7.0),
            ..FitConstraints::default()
        },
        |s| {
            matches!(s, ogeom::algo::Canonical::Sphere(s)
                if s.radius() == 7.0 && s.centre().distance(Point::new(0.0, 0.0, 0.0)) < 1e-6)
        },
    );
    let cone = ogeom::algo::make_cone(&mut model, Frame::WORLD, 6.0, 3.0, 10.0, T)
        .unwrap()
        .shape;
    refits_onto(
        &mut model,
        &cone,
        SurfaceKind::Cone,
        FitConstraints {
            axis: Some(z),
            ..FitConstraints::default()
        },
        |s| {
            matches!(s, ogeom::algo::Canonical::Cone(c)
                if near(c.half_angle(), (3.0_f64 / 10.0).atan())
                    && c.frame().z().vector().cross(Vector::Z).magnitude() < 1e-12)
        },
    );
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 20.0, 5.0, T)
        .unwrap()
        .shape;
    for constraints in [
        FitConstraints {
            axis: Some(z),
            ..FitConstraints::default()
        },
        FitConstraints {
            radius: Some(5.0),
            ..FitConstraints::default()
        },
    ] {
        refits_onto(&mut model, &ring, SurfaceKind::Torus, constraints, |s| {
            matches!(s, ogeom::algo::Canonical::Torus(t)
                    if near(t.major_radius(), 20.0) && near(t.minor_radius(), 5.0))
        });
    }
}

/// The steps that cannot be taken are refused by name, and leave the
/// regions as they were.
#[test]
fn steps_that_cannot_be_taken_are_refused_by_name() {
    let (r, h) = (10.0, 20.0);
    let mesh = faceted_cylinder(r, h, 10, 2);
    let mut regions = MeshRegions::find(&mesh, &MeshSolidOptions::default(), T).unwrap();
    let ends: Vec<RegionId> = regions
        .regions()
        .into_iter()
        .filter(|g| !side_facets(&regions).contains(&g.id))
        .map(|g| g.id)
        .collect();
    assert_eq!(ends.len(), 2);
    assert_eq!(
        regions.merge(ends[0], ends[1]),
        Err(RegionRefusal::NotAdjacent(ends[0], ends[1]))
    );
    assert_eq!(
        regions.merge(ends[0], ends[0]),
        Err(RegionRefusal::SameRegion(ends[0]))
    );
    let sides = merge_sides(&mut regions);
    let before = regions.regions();

    // A path up one ruling from rim to rim cuts the band open but leaves
    // it one piece.
    let ruling: Vec<u32> = [0.0, 10.0, 20.0]
        .map(|z| vertex_at(&regions, Point::new(r, 0.0, z)))
        .to_vec();
    assert_eq!(
        regions.split(sides, &ruling),
        Err(RegionRefusal::DoesNotCut { pieces: 1 })
    );
    // Two vertices on opposite sides share no edge.
    let across = [
        vertex_at(&regions, Point::new(r, 0.0, 0.0)),
        vertex_at(&regions, Point::new(-r, 0.0, 0.0)),
    ];
    assert_eq!(
        regions.split(sides, &across),
        Err(RegionRefusal::NotAnEdgeOfTheRegion(across[0], across[1]))
    );
    // The band is not a sphere, nor a cylinder a twentieth wider, nor one
    // about an axis a thousandth off.
    let free = FitConstraints::default();
    assert!(matches!(
        regions.fit(sides, SurfaceKind::Sphere, &free),
        Err(RegionRefusal::DoesNotVerify { .. })
    ));
    let wide = FitConstraints {
        radius: Some(r * 1.05),
        ..free
    };
    let Err(RegionRefusal::DoesNotVerify {
        deviation,
        distance,
        kind: SurfaceKind::Cylinder,
    }) = regions.fit(sides, SurfaceKind::Cylinder, &wide)
    else {
        panic!("a wider cylinder does not verify");
    };
    assert!(deviation > distance);
    let off = FitConstraints {
        axis: Some(Axis {
            location: Point::new(1e-3, 0.0, 0.0),
            direction: Direction::Z,
        }),
        ..free
    };
    assert!(matches!(
        regions.fit(sides, SurfaceKind::Cylinder, &off),
        Err(RegionRefusal::DoesNotVerify { .. })
    ));
    // A cone's radius changes along it, and a plane has no axis.
    assert_eq!(
        regions.fit(sides, SurfaceKind::Cone, &wide),
        Err(RegionRefusal::ConstraintDoesNotApply {
            kind: SurfaceKind::Cone,
            constraint: "radius",
        })
    );
    assert_eq!(
        regions.fit(sides, SurfaceKind::Plane, &off),
        Err(RegionRefusal::ConstraintDoesNotApply {
            kind: SurfaceKind::Plane,
            constraint: "axis",
        })
    );
    // A band round the axis is no disk, and a patch holds no radius.
    assert_eq!(
        regions.fit(sides, SurfaceKind::Patch, &free),
        Err(RegionRefusal::NotADisk)
    );
    assert_eq!(
        regions.fit(sides, SurfaceKind::Patch, &wide),
        Err(RegionRefusal::ConstraintDoesNotApply {
            kind: SurfaceKind::Patch,
            constraint: "radius",
        })
    );
    assert_eq!(regions.regions(), before);
}

/// The vertices of a band drawn one row high lie on its rim circles, and
/// so on a sphere through both as well as on the cylinder. A piece of the
/// band split off and fitted as that sphere verifies at every vertex, but
/// the sphere meets the cylinder beside it only on the rims, so no seam
/// between them can be placed along the rulings: the build facets both,
/// and names them. The solid is the polygonal prism the facets bound.
#[test]
fn a_region_an_edit_leaves_unbuildable_is_named_with_why() {
    let (r, h, n) = (10.0, 20.0, 64);
    let mesh = faceted_cylinder(r, h, n, 1);
    let options = MeshSolidOptions::default();
    let mut regions = MeshRegions::find(&mesh, &options, T).unwrap();
    let mut back = Model::new();
    let whole = regions.build(&mut back).unwrap();
    assert_eq!(kinds(&back, &whole.shape), [2, 1, 0, 0, 0, 0]);
    assert!(whole.report.fallbacks.is_empty(), "{:?}", whole.report);

    let band = regions.region_of(0).unwrap();
    let path = across_the_band(&regions, r, h, n, 0, 8);
    let (piece, rest) = regions.split(band, &path).unwrap();
    let deviation = regions
        .fit(piece, SurfaceKind::Sphere, &FitConstraints::default())
        .unwrap();
    assert!(deviation < 1e-12, "{deviation}");
    let mut back = Model::new();
    let out = regions.build(&mut back).unwrap();
    valid(&back, &out);
    let report = &out.report;
    assert_eq!(report.curved_faceted, report.fallbacks.len());
    for id in [piece, rest] {
        assert!(
            report.fallbacks.contains(&RegionFallback {
                region: id,
                reason: FallbackReason::BoundaryNotPlaced,
            }),
            "{report:?}"
        );
    }
    assert_eq!(kinds(&back, &out.shape), [66, 0, 0, 0, 0, 0]);
    let prism = f64::from(n) / 2.0 * r * r * (core::f64::consts::TAU / f64::from(n)).sin() * h;
    let v = volume(&back, &out.shape);
    assert!((v - prism).abs() < 1e-9 * prism, "{v} against {prism}");
}

/// A block whose edge is rounded with a radius of three, the round drawn
/// as two facets that turn 45 degrees: each facet is a plane of its own.
/// One row of the round's vertices stands a hundred-thousandth outward, so
/// a cylinder fitted freely through the round's three rows leans off the
/// faces beside it. Merged, the round is put on the cylinder tangent to
/// both faces, as recognition puts its own rounds, at the radius most of
/// its vertices give: the block's own round, and the solid's exact volume.
#[test]
fn a_merged_round_between_two_planes_is_tangent_to_both() {
    let (radius, length, rows) = (3.0, 10.0, 4_u32);
    let lift = 1e-5;
    let diagonal = core::f64::consts::FRAC_1_SQRT_2;
    let section = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 2.0),
        (
            7.0 + (radius + lift) * diagonal,
            2.0 + (radius + lift) * diagonal,
        ),
        (7.0, 5.0),
        (0.0, 5.0),
    ];
    let m = u32::try_from(section.len()).unwrap();
    let mut points = Vec::new();
    for j in 0..=rows {
        let y = length * f64::from(j) / f64::from(rows);
        points.extend(section.iter().map(|&(x, z)| Point::new(x, y, z)));
    }
    let at = |j: u32, k: u32| j * m + k % m;
    let mut polygons: Vec<Vec<u32>> = Vec::new();
    for j in 0..rows {
        for k in 0..m {
            polygons.push(vec![at(j, k), at(j + 1, k), at(j + 1, k + 1), at(j, k + 1)]);
        }
    }
    polygons.push((0..m).map(|k| at(0, k)).collect());
    polygons.push((0..m).rev().map(|k| at(rows, k)).collect());
    let mesh = from_polygons(points, &polygons);
    let exact = length * radius.mul_add(-radius * (1.0 - core::f64::consts::FRAC_PI_4), 50.0);

    let mut regions = MeshRegions::find(&mesh, &MeshSolidOptions::default(), T).unwrap();
    let facets: Vec<RegionId> = regions
        .regions()
        .into_iter()
        .filter(|g| match &g.surface {
            Some(ogeom::algo::Canonical::Plane(p)) => {
                let n = p.frame().z().vector();
                n.x.abs() > 0.1 && n.z.abs() > 0.1
            }
            _ => false,
        })
        .map(|g| g.id)
        .collect();
    assert_eq!(facets.len(), 2);
    let round = merge_all(&mut regions, facets);
    let region = regions.region(round).unwrap();
    let Some(ogeom::algo::Canonical::Cylinder(c)) = region.surface else {
        panic!("a cylinder: {region:?}");
    };
    assert!((c.radius() - radius).abs() < 1e-12, "{c:?}");
    let o = c.frame().origin();
    assert!(
        (o.x - 7.0).abs() < 1e-12 && (o.z - 2.0).abs() < 1e-12,
        "{c:?}"
    );
    assert!(region.deviation.unwrap() <= regions.distance());

    let mut back = Model::new();
    let out = regions.build(&mut back).unwrap();
    valid(&back, &out);
    assert!(out.report.fallbacks.is_empty(), "{:?}", out.report);
    assert_eq!(kinds(&back, &out.shape), [6, 1, 0, 0, 0, 0]);
    let v = volume(&back, &out.shape);
    assert!((v - exact).abs() < 1e-9 * exact, "{v} against {exact}");
}

/// A cylinder drawn with 64 facets round it is recognized whole. One facet
/// split off is a plane of two triangles with its corners on the
/// cylinder, which the build gives to the curved face beside it where
/// recognition left it. Split off by the caller it stays a face, and the
/// solid is the cylinder with that one flat cut along it.
#[test]
fn a_facet_the_caller_split_off_is_not_absorbed() {
    let (r, h, n) = (10.0, 20.0, 64);
    let mesh = faceted_cylinder(r, h, n, 1);
    let mut regions = MeshRegions::find(&mesh, &MeshSolidOptions::default(), T).unwrap();
    let band = regions.region_of(0).unwrap();
    let path = across_the_band(&regions, r, h, n, 0, 1);
    let (flat, rest) = regions.split(band, &path).unwrap();
    let facet = regions.region(flat).unwrap();
    assert_eq!(facet.triangles.len(), 2);
    assert!(matches!(
        facet.surface,
        Some(ogeom::algo::Canonical::Plane(_))
    ));
    assert!(matches!(
        regions.region(rest).unwrap().surface,
        Some(ogeom::algo::Canonical::Cylinder(_))
    ));
    let mut back = Model::new();
    let out = regions.build(&mut back).unwrap();
    valid(&back, &out);
    assert_eq!(kinds(&back, &out.shape), [3, 1, 0, 0, 0, 0]);
    let theta = core::f64::consts::TAU / f64::from(n);
    let segment = r * r / 2.0 * (theta - theta.sin());
    let exact = (core::f64::consts::PI * r * r - segment) * h;
    let v = volume(&back, &out.shape);
    assert!((v - exact).abs() < 1e-9 * exact, "{v} against {exact}");
}

/// A block ten wide whose top is the saddle `z = 5 + x y / 20`, drawn on a
/// grid of seven cells across. Its volume is the base's area times five
/// exactly: the saddle's rise over the square cancels.
fn saddle_block() -> (Triangulation, f64) {
    let (half, height, c, cells) = (5.0_f64, 5.0, 0.05_f64, 7_u32);
    let top = |x: f64, y: f64| c.mul_add(x * y, height);
    let step = 2.0 * half / f64::from(cells);
    let coordinate = |i: u32| step.mul_add(f64::from(i), -half);
    let mut points = Vec::new();
    for j in 0..=cells {
        for i in 0..=cells {
            let (x, y) = (coordinate(i), coordinate(j));
            points.push(Point::new(x, y, top(x, y)));
        }
    }
    let grid = |i: u32, j: u32| j * (cells + 1) + i;
    // The rim round the top, counterclockwise from above, each with the
    // vertex below it on the base.
    let mut rim: Vec<(u32, u32)> = Vec::new();
    let mut ring = Vec::new();
    ring.extend((0..cells).map(|i| (i, 0)));
    ring.extend((0..cells).map(|j| (cells, j)));
    ring.extend((0..cells).map(|i| (cells - i, cells)));
    ring.extend((0..cells).map(|j| (0, cells - j)));
    for (i, j) in ring {
        let below = u32::try_from(points.len()).unwrap();
        points.push(Point::new(coordinate(i), coordinate(j), 0.0));
        rim.push((grid(i, j), below));
    }
    let centre = u32::try_from(points.len()).unwrap();
    points.push(Point::new(0.0, 0.0, 0.0));
    let mut polygons: Vec<Vec<u32>> = Vec::new();
    for j in 0..cells {
        for i in 0..cells {
            polygons.push(vec![grid(i, j), grid(i + 1, j), grid(i + 1, j + 1)]);
            polygons.push(vec![grid(i, j), grid(i + 1, j + 1), grid(i, j + 1)]);
        }
    }
    for k in 0..rim.len() {
        let (a, a_low) = rim[k];
        let (b, b_low) = rim[(k + 1) % rim.len()];
        polygons.push(vec![a_low, b_low, b, a]);
        polygons.push(vec![centre, b_low, a_low]);
    }
    let base = 2.0 * half;
    (from_polygons(points, &polygons), base * base * height)
}

/// The facets of the saddle block's top: its regions on planes facing up.
fn top_facets(regions: &MeshRegions) -> Vec<RegionId> {
    regions
        .regions()
        .into_iter()
        .filter(|g| match &g.surface {
            Some(ogeom::algo::Canonical::Plane(p)) => p.frame().z().vector().z > 0.5,
            _ => false,
        })
        .map(|g| g.id)
        .collect()
}

/// Found with a crease angle so small that every turn of the saddle is a
/// crease, the block's top is a facet per triangle. Merged, no plane or
/// canonical surface holds the facets, and the merge falls back to a
/// patch over the disk they make; the block builds with its exact volume,
/// within the coplanar distance over the base's area. With the patches off
/// the merged region has no surface, and a patch fitted to it on request
/// builds the same. A patch is refused by name where the region folds over
/// an edge of the block or is a strip with no vertex inside.
#[test]
fn a_merged_smooth_region_falls_back_to_a_patch() {
    let (mesh, exact) = saddle_block();
    let creased = MeshSolidOptions {
        crease: 0.01,
        ..MeshSolidOptions::default()
    };
    let check = |regions: &MeshRegions| {
        let mut back = Model::new();
        let out = regions.build(&mut back).unwrap();
        valid(&back, &out);
        assert_eq!(out.report.patch_faces, 1, "{:?}", out.report);
        assert_eq!(kinds(&back, &out.shape), [5, 0, 0, 0, 0, 1]);
        let v = volume(&back, &out.shape);
        let allowance = 100.0 * regions.distance();
        assert!((v - exact).abs() < allowance, "{v} against {exact}");
    };

    let mut regions = MeshRegions::find(&mesh, &creased, T).unwrap();
    let facets = top_facets(&regions);
    assert_eq!(facets.len(), 98);
    let top = merge_all(&mut regions, facets);
    let region = regions.region(top).unwrap();
    assert!(
        matches!(region.surface, Some(ogeom::algo::Canonical::Swept(_))),
        "{region:?}"
    );
    assert!(region.deviation.unwrap() <= regions.distance());
    check(&regions);

    // A side of the block is a strip with every vertex on its rim, and
    // the top merged with it folds over the block's edge.
    let side = regions
        .regions()
        .into_iter()
        .find(|g| match &g.surface {
            Some(ogeom::algo::Canonical::Plane(p)) => p.frame().z().vector().x > 0.5,
            _ => false,
        })
        .unwrap()
        .id;
    let free = FitConstraints::default();
    assert_eq!(
        regions.fit(side, SurfaceKind::Patch, &free),
        Err(RegionRefusal::TooNarrowForAPatch)
    );
    let mut folded = regions.clone();
    let both = folded.merge(top, side).unwrap();
    assert_eq!(folded.region(both).unwrap().surface, None);
    assert_eq!(
        folded.fit(both, SurfaceKind::Patch, &free),
        Err(RegionRefusal::PatchDoesNotVerify)
    );

    let unpatched = MeshSolidOptions {
        patches: false,
        ..creased
    };
    let mut regions = MeshRegions::find(&mesh, &unpatched, T).unwrap();
    let facets = top_facets(&regions);
    let top = merge_all(&mut regions, facets);
    assert_eq!(regions.region(top).unwrap().surface, None);
    let deviation = regions.fit(top, SurfaceKind::Patch, &free).unwrap();
    assert!(deviation <= regions.distance(), "{deviation}");
    check(&regions);
}
