//! The five configurations §A of `docs/PLAN.md` names, each measured against
//! a closed form.
//!
//! They are one family: every one of them is the classifier being asked where
//! a piece stands at a place where the question has no answer — on a
//! coincidence, at a chart degeneracy, in a cusp, on a tangency. They are
//! written here before the interference table exists, so that what the table
//! buys is a measurement and not an impression.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> Deflection {
    Deflection::with_chord(1e-3).unwrap()
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

fn at(p: Point) -> Frame {
    Frame::new(p, Direction::Z, Direction::X, T).unwrap()
}

/// **A1 — a tool flush with the part.** A block with a bore through it,
/// refilled by the very cylinder that cut it. The tool's wall coincides with
/// the bore's wall over its whole length and its caps are flush with the two
/// faces the bore broke, so every piece of the tool lies *on* the part rather
/// than in or out of it.
///
/// The measured claim is that the refilled block is the block again: same
/// volume, and one closed shell.
#[test]
fn a1_a_tool_flush_with_the_part_refills_the_bore() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    // A bore straight through, its axis at the block's centre. Long enough
    // that its walls break both the bottom and the top face.
    let drill =
        ogeom::algo::make_cylinder(&mut model, at(Point::new(10.0, 10.0, -5.0)), 3.0, 20.0, T)
            .unwrap()
            .shape;
    let drilled = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;
    let pi = core::f64::consts::PI;
    let bore = pi * 9.0 * 10.0;
    assert!(
        (volume(&model, &drilled) - (4000.0 - bore)).abs() < bore * 1e-3,
        "the bore came out: {}",
        volume(&model, &drilled)
    );

    // The plug: exactly the bore, ends flush with the faces it broke.
    let plug =
        ogeom::algo::make_cylinder(&mut model, at(Point::new(10.0, 10.0, 0.0)), 3.0, 10.0, T)
            .unwrap()
            .shape;
    let refilled = ogeom::boolean::fuse(&mut model, &drilled, &plug, T)
        .unwrap()
        .shape;
    assert!(
        (volume(&model, &refilled) - 4000.0).abs() < 4000.0 * 1e-3,
        "the block is whole again: {}",
        volume(&model, &refilled)
    );
    assert!(
        ogeom::algo::check(&model, &refilled, T).unwrap().is_valid(),
        "and it is a valid solid"
    );
}

/// **A2 — a section through a chart pole.** A plane through a ball's own axis
/// ends its section exactly at both poles, where the sphere's chart
/// degenerates: a whole line of the chart is one point in space, and the
/// pieces either side of the section meet there.
///
/// Measured against the closed form for a half ball.
#[test]
fn a2_a_section_through_both_poles_halves_the_ball() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 4.0, T)
        .unwrap()
        .shape;
    // A block filling x ≥ 0 around the ball. Its face x = 0 contains the
    // ball's axis, so the section runs pole to pole.
    let half = ogeom::algo::make_box(
        &mut model,
        at(Point::new(0.0, -10.0, -10.0)),
        (10.0, 20.0, 20.0),
        T,
    )
    .unwrap()
    .shape;
    let pi = core::f64::consts::PI;
    let ball_volume = 4.0 / 3.0 * pi * 64.0;

    let kept = ogeom::boolean::common(&mut model, &ball, &half, T)
        .unwrap()
        .shape;
    assert!(
        (volume(&model, &kept) - ball_volume / 2.0).abs() < ball_volume * 1e-3,
        "half the ball: {}",
        volume(&model, &kept)
    );

    let rest = ogeom::boolean::cut(&mut model, &ball, &half, T)
        .unwrap()
        .shape;
    assert!(
        (volume(&model, &rest) - ball_volume / 2.0).abs() < ball_volume * 1e-3,
        "and the other half: {}",
        volume(&model, &rest)
    );
}

/// **A3 — a ball on a box's corner vertex.** A vertex-against-face
/// interference is a first-class thing; it has no business surfacing as an
/// anomaly during face splitting.
///
/// Three placements of the same pair, each a different way for the corner to
/// meet the ball. Touching the vertex from outside bounds no material, so the
/// fuse is both volumes and the cut takes nothing. Centred *on* the vertex,
/// the three faces meeting there cut the ball through its own centre and the
/// closed form is an octant. And with the vertex exactly on the sphere, the
/// ball straddles the corner and the two halves it is cut into must still
/// account for all of it.
#[test]
fn a3_a_ball_on_a_corner_vertex_bounds_nothing() {
    let pi = core::f64::consts::PI;
    let corner = Point::new(10.0, 10.0, 10.0);
    let block_of = |model: &mut Model| {
        ogeom::algo::make_box(model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape
    };

    // Tangent at the vertex from outside: nothing is bounded.
    {
        let mut model = Model::new();
        let block = block_of(&mut model);
        let d = 3.0 / 3.0_f64.sqrt();
        let ball = ogeom::algo::make_sphere(
            &mut model,
            at(Point::new(10.0 + d, 10.0 + d, 10.0 + d)),
            3.0,
            T,
        )
        .unwrap()
        .shape;
        let ball_volume = 4.0 / 3.0 * pi * 27.0;
        let fused = ogeom::boolean::fuse(&mut model, &block, &ball, T)
            .unwrap()
            .shape;
        assert!(
            (volume(&model, &fused) - (1000.0 + ball_volume)).abs() < 1.0,
            "both volumes, joined at the corner they share: {}",
            volume(&model, &fused)
        );
        let cut = ogeom::boolean::cut(&mut model, &block, &ball, T)
            .unwrap()
            .shape;
        assert!(
            (volume(&model, &cut) - 1000.0).abs() < 1e-6,
            "a vertex removes nothing: {}",
            volume(&model, &cut)
        );
    }

    // Centred on the vertex: the box's three faces there cut the ball through
    // its centre, so exactly one octant of it is inside.
    {
        let mut model = Model::new();
        let block = block_of(&mut model);
        let ball = ogeom::algo::make_sphere(&mut model, at(corner), 3.0, T)
            .unwrap()
            .shape;
        let octant = 4.0 / 3.0 * pi * 27.0 / 8.0;
        let shared = ogeom::boolean::common(&mut model, &block, &ball, T)
            .unwrap()
            .shape;
        assert!(
            (volume(&model, &shared) - octant).abs() < octant * 5e-3,
            "one octant of the ball: {} vs {octant}",
            volume(&model, &shared)
        );
        let cut = ogeom::boolean::cut(&mut model, &block, &ball, T)
            .unwrap()
            .shape;
        assert!(
            (volume(&model, &cut) - (1000.0 - octant)).abs() < octant * 5e-3,
            "and the block less that octant: {}",
            volume(&model, &cut)
        );
    }

    // The vertex exactly on the sphere, the ball straddling the corner. The
    // closed form here is awkward — three overlapping caps meeting at the
    // point — so the claim is an independent one: the two pieces the box cuts
    // the ball into are the whole ball, and neither is empty.
    {
        let mut model = Model::new();
        let block = block_of(&mut model);
        let d = 3.0 / 3.0_f64.sqrt();
        let ball = ogeom::algo::make_sphere(
            &mut model,
            at(Point::new(10.0 - d, 10.0 - d, 10.0 - d)),
            3.0,
            T,
        )
        .unwrap()
        .shape;
        let ball_volume = 4.0 / 3.0 * pi * 27.0;
        let inside = ogeom::boolean::common(&mut model, &block, &ball, T)
            .unwrap()
            .shape;
        let outside = ogeom::boolean::cut(&mut model, &ball, &block, T)
            .unwrap()
            .shape;
        let (vi, vo) = (volume(&model, &inside), volume(&model, &outside));
        assert!(vi > 0.0 && vo > 0.0, "both pieces are real: {vi} and {vo}");
        assert!(
            (vi + vo - ball_volume).abs() < ball_volume * 5e-3,
            "and together they are the ball: {vi} + {vo} against {ball_volume}"
        );
    }
}

/// **A4 — a section through tangential contact.** The textbook half-section:
/// a plane through a bore's axis. It meets the bore's wall along the wall's
/// own rulings, and it meets the bore's rim circles at exactly the points
/// those rulings end — a place where the crossing question is degenerate from
/// every side at once.
///
/// The measured claim is the half of the drilled block, and — because this is
/// what the drawing feature D2 needs — a section whose curves are there.
#[test]
fn a4_an_on_axis_half_section_of_a_bore() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let drill =
        ogeom::algo::make_cylinder(&mut model, at(Point::new(10.0, 10.0, -5.0)), 3.0, 20.0, T)
            .unwrap()
            .shape;
    let drilled = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;

    // The cutting solid's face y = 10 passes straight down the bore's axis.
    let knife = ogeom::algo::make_box(
        &mut model,
        at(Point::new(-5.0, 10.0, -5.0)),
        (30.0, 20.0, 20.0),
        T,
    )
    .unwrap()
    .shape;
    let pi = core::f64::consts::PI;
    let whole = 4000.0 - pi * 9.0 * 10.0;
    let half = ogeom::boolean::cut(&mut model, &drilled, &knife, T)
        .unwrap()
        .shape;
    assert!(
        (volume(&model, &half) - whole / 2.0).abs() < whole * 1e-3,
        "half the drilled block: {}",
        volume(&model, &half)
    );

    // And the section itself carries the two rulings the plane cuts from the
    // bore's wall, each the bore's full depth.
    let section = ogeom::boolean::section(&mut model, &drilled, &knife, T)
        .unwrap()
        .shape;
    let edges = ogeom::topo::explore_unique(&model, &section, ShapeType::Edge).unwrap();
    let mut rulings = 0;
    for e in &edges {
        let length = ogeom::algo::linear_properties(&model, e, fine(), T)
            .unwrap()
            .mass;
        let Some(mid) = ogeom::algo::shape_bounds(&model, e, T).unwrap().centre() else {
            continue;
        };
        // A ruling of the bore stands at radius 3 from the axis and runs the
        // full depth.
        let radius = ((mid.x - 10.0) * (mid.x - 10.0) + (mid.y - 10.0) * (mid.y - 10.0)).sqrt();
        if (radius - 3.0).abs() < 1e-6 && (length - 10.0).abs() < 1e-6 {
            rulings += 1;
        }
    }
    assert_eq!(rulings, 2, "the plane cuts the wall along two rulings");
}

/// **A5 — a shell around a three-cylinder tip.** Three mutually perpendicular
/// cylinders of one radius meet at a box's corner: the three fillet surfaces a
/// corner blend would leave. The tool that clears the leftover spike is the
/// corner block minus the ball, and removing it is what fails.
///
/// The measured claim is the closed form for the block's corner rounded to a
/// sphere octant: the corner block of side `r` less the spike, which is that
/// block less an eighth of the ball.
#[test]
fn a5_a_corner_block_less_a_ball_clears_the_spike() {
    let mut model = Model::new();
    let r = 3.0;
    // The corner block: the cube of side r seated at the box's corner,
    // running from the ball's centre out past the corner.
    let block = ogeom::algo::make_box(
        &mut model,
        at(Point::new(10.0 - r, 10.0 - r, 10.0 - r)),
        (r, r, r),
        T,
    )
    .unwrap()
    .shape;
    let ball = ogeom::algo::make_sphere(
        &mut model,
        at(Point::new(10.0 - r, 10.0 - r, 10.0 - r)),
        r,
        T,
    )
    .unwrap()
    .shape;
    let pi = core::f64::consts::PI;
    let octant = 4.0 / 3.0 * pi * r * r * r / 8.0;
    let spike = ogeom::boolean::cut(&mut model, &block, &ball, T)
        .unwrap()
        .shape;
    let want = r * r * r - octant;
    // The volume is measured off a mesh, and a mesh of a sphere is inscribed
    // in it: at chord `δ` on radius `r` the deficit runs to about `3δ/r`, so
    // the instrument is refined until its own error is well under the claim.
    // Measured at two deflections a decade apart, which is what says the
    // remaining gap is the mesh and not the boolean.
    let mut previous = f64::INFINITY;
    for chord in [1e-3, 1e-4] {
        let fine = Deflection::with_chord(chord).unwrap();
        let measured = ogeom::algo::volume_properties(&model, &spike, fine, T)
            .unwrap()
            .mass;
        let error = (measured - want).abs() / want;
        assert!(
            error < previous,
            "refining the mesh brings the measurement closer: {measured} vs {want}"
        );
        assert!(
            error < chord * 4.0 / r,
            "the corner block less the ball's octant, within the mesh's own \
             deficit at chord {chord}: {measured} vs {want}"
        );
        previous = error;
    }
    assert!(
        ogeom::algo::check(&model, &spike, T).unwrap().is_valid(),
        "and the spike is a valid solid"
    );
}
