//! §E1 of `docs/PLAN.md`: one pick structure over several deflections.
//!
//! The claim under test is not that the hierarchy is faster — that is what it
//! is *for*, and a timing is not a proof of anything. The claim is that it
//! changes nothing: what the descent returns is what the finest level alone
//! returns, ray for ray, hit for hit. A structure that skipped work and got a
//! different answer would be worse than useless.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::select::{PickHierarchy, Pickable, Ray};
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn at(p: Point) -> Frame {
    Frame::new(p, Direction::Z, Direction::X, T).unwrap()
}

/// A part with curvature in it, so the deflections genuinely differ: a block
/// with a bore, a boss on top and a ball taken out of one corner.
fn part(model: &mut Model) -> Shape {
    let block = ogeom::algo::make_box(model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let bore = ogeom::algo::make_cylinder(model, at(Point::new(6.0, 6.0, -5.0)), 2.5, 20.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(model, &block, &bore, T).unwrap().shape;
    let boss = ogeom::algo::make_cylinder(model, at(Point::new(14.0, 14.0, 10.0)), 3.0, 4.0, T)
        .unwrap()
        .shape;
    ogeom::boolean::fuse(model, &drilled, &boss, T)
        .unwrap()
        .shape
}

/// A fan of rays that samples the part from several directions, including
/// ones that graze it.
fn rays() -> Vec<Ray> {
    let mut out = Vec::new();
    for i in 0..9 {
        for j in 0..9 {
            let (x, y) = (f64::from(i) * 2.4, f64::from(j) * 2.4);
            out.push(Ray {
                origin: Point::new(x, y, 40.0),
                direction: Vector::new(0.0, 0.0, -1.0),
            });
            out.push(Ray {
                origin: Point::new(-30.0, x, y * 0.6),
                direction: Vector::new(1.0, 0.1, 0.05),
            });
        }
    }
    out
}

/// The descent returns the finest level's own answer.
#[test]
fn a_hierarchy_answers_exactly_as_its_finest_level_does() {
    let mut model = Model::new();
    let shape = part(&mut model);

    let chords = [1.0, 0.25, 0.05];
    let deflections: Vec<Deflection> = chords
        .iter()
        .map(|c| Deflection::with_chord(*c).unwrap())
        .collect();
    let hierarchy = PickHierarchy::build(&model, &shape, &deflections, T).unwrap();
    assert_eq!(hierarchy.level_count(), 3);
    // Coarsest first.
    assert!(hierarchy.chord(0).unwrap() > hierarchy.chord(2).unwrap());

    let finest = Pickable::build(&model, &shape, Deflection::with_chord(0.05).unwrap(), T).unwrap();

    let mut struck = 0;
    for ray in rays() {
        let theirs = finest.pick(ray, 0.2);
        let ours = hierarchy.pick(ray, 0.2);
        assert_eq!(
            ours.len(),
            theirs.len(),
            "the same hits at {:?}: {} against {}",
            ray.origin,
            ours.len(),
            theirs.len()
        );
        for (a, b) in ours.iter().zip(&theirs) {
            assert!(a.shape.is_same(&b.shape), "and the same sub-shape");
            assert_eq!(a.kind, b.kind);
            assert!(
                (a.distance - b.distance).abs() < 1e-12,
                "and the same depth: {} against {}",
                a.distance,
                b.distance
            );
        }
        struck += ours.len();
    }
    assert!(struck > 100, "the fan actually hits the part: {struck}");
}

/// The same for the refined answers, which are the exact ones.
#[test]
fn refining_through_the_hierarchy_lands_on_the_same_surface_points() {
    let mut model = Model::new();
    let shape = part(&mut model);
    let deflections = [
        Deflection::with_chord(0.8).unwrap(),
        Deflection::with_chord(0.05).unwrap(),
    ];
    let hierarchy = PickHierarchy::build(&model, &shape, &deflections, T).unwrap();
    let finest = Pickable::build(&model, &shape, deflections[1], T).unwrap();

    let mut refined = 0;
    for ray in rays() {
        let ours = hierarchy.pick_refined(ray, 0.0, T);
        let theirs = finest.pick_refined(ray, 0.0, T);
        assert_eq!(ours.len(), theirs.len());
        for (a, b) in ours.iter().zip(&theirs) {
            assert_eq!(a.refined, b.refined);
            assert!(
                a.position.distance(b.position) < 1e-12,
                "the same exact point: {:?} against {:?}",
                a.position,
                b.position
            );
            if a.refined {
                refined += 1;
            }
        }
    }
    assert!(refined > 50, "the exact refinement did run: {refined}");
}

/// A level is chosen by the detail a view wants, without building anything.
#[test]
fn a_view_asks_for_the_detail_it_needs_and_gets_a_level_that_exists() {
    let mut model = Model::new();
    let shape = part(&mut model);
    let deflections = [
        Deflection::with_chord(1.0).unwrap(),
        Deflection::with_chord(0.1).unwrap(),
        Deflection::with_chord(0.01).unwrap(),
    ];
    let hierarchy = PickHierarchy::build(&model, &shape, &deflections, T).unwrap();

    // A distant view asks for little and is given the coarsest that serves.
    let far = hierarchy.for_chord(2.0);
    assert_eq!(far.triangle_count(), hierarchy.coarsest().triangle_count());
    // A close one asks for more than any level has and gets the finest.
    let close = hierarchy.for_chord(1e-6);
    assert_eq!(close.triangle_count(), hierarchy.finest().triangle_count());
    // And the levels genuinely differ, or there was nothing to choose.
    assert!(
        hierarchy.finest().triangle_count() > hierarchy.coarsest().triangle_count() * 2,
        "the levels are different scenes: {} against {}",
        hierarchy.finest().triangle_count(),
        hierarchy.coarsest().triangle_count()
    );

    // Every level names the same faces in the same order — which is what lets
    // an answer found at one level be carried to another.
    let count = hierarchy.coarsest().face_count();
    for index in 0..hierarchy.level_count() {
        assert_eq!(hierarchy.level(index).unwrap().face_count(), count);
    }
}

/// And the descent actually rules things out. Not a timing — a count: for a
/// ray down the part, the coarse level admits a small fraction of the faces,
/// which is the whole reason the structure exists.
#[test]
fn the_coarse_level_rules_out_most_of_the_scene() {
    let mut model = Model::new();
    let shape = part(&mut model);
    let coarse = Pickable::build(&model, &shape, Deflection::with_chord(1.0).unwrap(), T).unwrap();
    let fine = Pickable::build(&model, &shape, Deflection::with_chord(0.05).unwrap(), T).unwrap();
    let total = coarse.face_count();
    assert!(total >= 8, "the part has faces to rule out: {total}");

    let mut admitted = 0;
    let mut asked = 0;
    for ray in rays() {
        let near = coarse.faces_near(ray, 1.0 + 0.05);
        admitted += near.iter().filter(|x| **x).count();
        asked += total;
        // And nothing the fine level hits is ever ruled out — the property
        // the descent stands on, checked ray by ray rather than argued.
        for hit in fine.pick(ray, 0.0) {
            let owner = fine.face_index(hit.triangle).unwrap();
            assert!(
                near[owner],
                "face {owner} is hit but was ruled out at {:?}",
                ray.origin
            );
        }
    }
    #[allow(clippy::cast_precision_loss, reason = "counts, far below the mantissa")]
    let fraction = admitted as f64 / asked as f64;
    assert!(
        fraction < 0.5,
        "the coarse level rules most of the scene out: {fraction} admitted"
    );
}

/// A hierarchy of no levels answers nothing, and says so.
#[test]
fn a_hierarchy_with_no_levels_is_refused() {
    let mut model = Model::new();
    let shape = part(&mut model);
    assert!(PickHierarchy::build(&model, &shape, &[], T).is_err());
}
