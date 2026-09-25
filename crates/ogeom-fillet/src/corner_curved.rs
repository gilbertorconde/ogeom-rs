//! Rounding a vertex where a curved face meets: the one-ball corner on any
//! three hosts.
//!
//! Three faces meet at the vertex, any of them curved (a drum's wall where
//! a box's side and top cut it). The ball that touches all three from
//! inside is found by walking its centre: from a guess, each host's nearest
//! point and tangent plane there give three planes, the point a radius in
//! from those planes is the next guess, and the walk settles where every
//! host's foot point stands a radius from the centre along that host's
//! normal. Three planes through the centre, each square to one of the
//! corner's edges at the edge's point nearest the centre, bound the
//! corner's compartment: the cone they make toward the vertex, truncated
//! past it. The compartment clipped to the solid is the corner's block,
//! the block less the ball is the spike the rounded corner sheds, and the
//! cut sheds it. The solid itself supplies the host sides, so a curved
//! host is spoken exactly as a planar one.
//!
//! Each band that follows ends on the plane square to its edge through the
//! ball's centre, where its ball is the corner's ball.

use crate::support::planar_face;
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Orientation, Shape, ShapeType};

/// A host of the corner: its placed surface and the sign that turns the
/// surface's own normal outward.
struct Host {
    faces: Vec<Shape>,
    surface: SurfaceGeometry,
    outward: f64,
}

/// Round the convex vertex `vertex` of `solid`, at `corner`, whose three
/// faces are not all planes.
pub(crate) fn curved_corner(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    corner: Point,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    // The hosts are the faces whose surfaces pass through the corner, read
    // from geometry as the planar tool reads its planes: the bands of an
    // earlier call trim the hosts back from the vertex and consume it, but
    // every host's surface still holds the corner. Pieces of one surface
    // (trims a boolean split apart) are one host, told by their normals.
    let mut hosts: Vec<Host> = Vec::new();
    let mut normals_at: Vec<Vector> = Vec::new();
    for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face().cloned()) else {
            continue;
        };
        let Some(stored) = model.geometry().surface(data.surface).cloned() else {
            continue;
        };
        use ogeom_geom::Transformable as _;
        let surface = stored.transformed(&face.transform(model.datums())?, tol)?;
        let Ok(projection) = ogeom_algo::project_on_surface(&surface, corner, 32, tol) else {
            continue;
        };
        if projection.point.distance(corner) > tol.confusion() * 100.0 {
            continue;
        }
        let (u, v) = projection.parameters;
        let (du, dv) = surface.d1_at(u, v, tol)?;
        let n = du.cross(dv);
        let m = n.magnitude();
        if m <= tol.angular() {
            continue;
        }
        let n = n / m;
        if normals_at
            .iter()
            .any(|held| held.cross(n).magnitude() < tol.angular() * 10.0)
        {
            let index = normals_at
                .iter()
                .position(|held| held.cross(n).magnitude() < tol.angular() * 10.0)
                .unwrap_or(0);
            hosts[index].faces.push(face);
            continue;
        }
        let outward = if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        normals_at.push(n);
        hosts.push(Host {
            faces: vec![face],
            surface,
            outward,
        });
    }
    if hosts.len() != 3 {
        ogeom_bail!(
            Construction,
            "a curved corner rounds where three faces meet; {} surfaces pass \
             through this vertex; see docs/PARITY.md, fillet.edge-blends",
            hosts.len()
        );
    }

    let centre = ball_centre(&hosts, corner, radius, tol)?;
    // The ball must touch each host inside its face, and sit in the
    // material: a concave corner's ball would add material, which a cut
    // cannot say.
    let deflection = ogeom_mesh::Deflection {
        chord: (radius * 1e-3).max(tol.confusion() * 1e2),
        ..ogeom_mesh::Deflection::default()
    };
    let mut touches: Vec<Point> = Vec::with_capacity(3);
    for host in &hosts {
        let foot = ogeom_algo::project_on_surface(&host.surface, centre, 32, tol)?.point;
        let mut seated = false;
        for face in &host.faces {
            seated |= ogeom_algo::classify_on_face(model, face, foot, deflection, tol)?
                != ogeom_algo::Containment::Out;
        }
        if !seated {
            ogeom_bail!(
                Construction,
                "a ball of radius {radius} does not seat inside every face of \
                 the corner; the corner is smaller than the ball"
            );
        }
        touches.push(foot);
    }
    if ogeom_algo::classify_in_solid(model, solid, centre, deflection, tol)?
        != ogeom_algo::Containment::In
    {
        ogeom_bail!(
            Construction,
            "the corner is concave; a ball rounding it adds material, which \
             the corner tool's cut cannot say"
        );
    }

    // The compartment's three planes, one per pair of hosts: through the
    // centre and the ball's touch points on the two, which is the plane of
    // the end section of the band between them, where its ball is this
    // ball. Two of them meet along the line from the centre through the
    // touch point they share, so each edge of the compartment's cone
    // passes exactly through a touch point.
    let toward = corner - centre;
    let mut normals: Vec<Vector> = Vec::with_capacity(3);
    for (i, j) in [(0, 1), (1, 2), (2, 0)] {
        let n = (touches[i] - centre).cross(touches[j] - centre);
        let m = n.magnitude();
        if m <= tol.angular() * radius * radius {
            ogeom_bail!(
                Construction,
                "two of the corner's faces meet tangentially; there is no \
                 crease between them to end a band on"
            );
        }
        let n = n / m;
        let facing = n.dot(toward);
        if facing.abs() <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "an edge of the corner runs square to the way to its vertex"
            );
        }
        normals.push(if facing > 0.0 { n } else { -n });
    }
    let axis = toward / toward.magnitude();
    // The cone's edges, each where two of the planes meet, turned into the
    // third plane's side.
    let mut rays: Vec<Vector> = Vec::with_capacity(3);
    for (j, k, l) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
        let e = normals[j].cross(normals[k]);
        let m = e.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "two of the corner's edges run parallel");
        }
        let e = e / m;
        let e = if e.dot(normals[l]) >= 0.0 { e } else { -e };
        if e.dot(axis) <= tol.angular() {
            ogeom_bail!(
                Construction,
                "the corner's compartment opens wider than a half space"
            );
        }
        rays.push(e);
    }
    // Truncated a fifth past the vertex along the axis, which holds the
    // corner's spike and reaches no further into the part.
    let reach = toward.dot(axis) * 1.2 + radius * 0.1;
    let tips: Vec<Point> = rays
        .iter()
        .map(|e| centre + *e * (reach / e.dot(axis)))
        .collect();
    // Each triangle wound counter-clockwise about its outward normal: a
    // side plane faces away from the cone's inside, against its own
    // normal, and the cap faces along the axis.
    let triangle = |model: &mut Model, corners: [Point; 3], outward: Vector| {
        let turn = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
        let wound = if turn.dot(outward) >= 0.0 {
            corners
        } else {
            [corners[0], corners[2], corners[1]]
        };
        planar_face(model, &wound, outward, tol)
    };
    // Ray 0 lies in planes 0 and 1, ray 1 in 1 and 2, ray 2 in 2 and 0.
    let faces = vec![
        triangle(model, [centre, tips[2], tips[0]], -normals[0])?,
        triangle(model, [centre, tips[0], tips[1]], -normals[1])?,
        triangle(model, [centre, tips[1], tips[2]], -normals[2])?,
        triangle(model, [tips[0], tips[1], tips[2]], axis)?,
    ];
    let sewn = ogeom_algo::sew(model, &faces, tol)?;
    let [shell] = sewn.shells.as_slice() else {
        ogeom_bail!(Construction, "the corner's compartment did not close");
    };
    if !ogeom_algo::is_shell_closed(model, shell)? {
        ogeom_bail!(Construction, "the corner's compartment did not close");
    }
    let compartment = ogeom_algo::make_solid(model, std::slice::from_ref(shell))?.shape;

    // The ball, its poles square to the axis and its seam turned from the
    // corner, so neither stands in the patch the compartment keeps.
    let pole = {
        let p = axis.cross(rays[0]);
        let m = p.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "the corner's compartment has no width");
        }
        p / m
    };
    let frame = Frame::new(
        centre,
        Direction::new(pole, tol)?,
        Direction::new(-axis, tol)?,
        tol,
    )?;
    let ball = ogeom_algo::make_sphere(model, frame, radius, tol)?.shape;
    // The block first: the compartment clipped to the solid. Each of the
    // cone's edges crosses a host square to it at that host's touch point,
    // so the touch points are the block's own vertices, and the ball's
    // tangency to each host lands on a vertex rather than being solved
    // where a host would cut the ball's rim in a tangent line.
    let block = ogeom_bool::common(model, solid, &compartment, tol)?;
    let spike = ogeom_bool::cut(model, &block.shape, &ball, tol)?;
    let rounded = ogeom_bool::cut(model, solid, &spike.shape, tol)?;
    let mut built = Built {
        shape: rounded.shape,
        history: block.history.then(&spike.history).then(&rounded.history),
    };
    built.history.modify(vertex, built.shape.clone());
    Ok(built)
}

/// The centre of the ball a radius in from all three hosts, walked from
/// the plane corner's own answer.
fn ball_centre(hosts: &[Host], corner: Point, radius: f64, tol: Tolerances) -> OgeomResult<Point> {
    let inward_at = |host: &Host, p: Point| -> OgeomResult<(Point, Vector)> {
        let projection = ogeom_algo::project_on_surface(&host.surface, p, 32, tol)?;
        let (u, v) = projection.parameters;
        let (du, dv) = host.surface.d1_at(u, v, tol)?;
        let n = du.cross(dv);
        let m = n.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "a corner host has no normal at the ball");
        }
        Ok((projection.point, n / m * (-host.outward)))
    };
    // m_k · (c − p_k) = radius, for each host's foot p_k and inward normal
    // m_k: three planes, one point.
    let solve = |planes: &[(Point, Vector)]| -> OgeomResult<Point> {
        let a = [
            planes[0].1.to_array(),
            planes[1].1.to_array(),
            planes[2].1.to_array(),
        ];
        let b = [
            radius + planes[0].1.dot(planes[0].0.to_vector()),
            radius + planes[1].1.dot(planes[1].0.to_vector()),
            radius + planes[2].1.dot(planes[2].0.to_vector()),
        ];
        // Cramer's rule: three planes, one point.
        let det = |m: [[f64; 3]; 3]| {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let d = det(a);
        if d.abs() <= tol.angular() {
            ogeom_bail!(
                Construction,
                "the corner's tangent planes do not meet in one point"
            );
        }
        let column = |c: usize| {
            let mut m = a;
            for (row, value) in m.iter_mut().zip(b) {
                row[c] = value;
            }
            det(m) / d
        };
        Ok(Point::new(column(0), column(1), column(2)))
    };
    let mut planes = Vec::with_capacity(3);
    for host in hosts {
        planes.push(inward_at(host, corner)?);
    }
    let mut centre = solve(&planes)?;
    for _ in 0..100 {
        let mut next_planes = Vec::with_capacity(3);
        for host in hosts {
            next_planes.push(inward_at(host, centre)?);
        }
        let next = solve(&next_planes)?;
        let moved = next.distance(centre);
        centre = next;
        if moved <= tol.confusion() * 1e-3 {
            return Ok(centre);
        }
    }
    ogeom_bail!(
        NotDone,
        "the corner ball's centre did not settle a radius in from its hosts"
    )
}
