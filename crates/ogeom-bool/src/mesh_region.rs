//! Cutting a solid by a region of a mesh: a break surface, extended past
//! the solid and thickened into a slab, cut away so the pieces on one side
//! of it are kept.
//!
//! The region is seen along its average normal, the sum of its triangles'
//! area vectors. Where every triangle faces that way and the rims meet
//! nothing but themselves, the region is the graph of a height over its
//! average plane: the triangles covering a point of the plane number the
//! times the rims wind about it, each triangle's winding about the point
//! one inside it and nothing outside, and with every winding one or nothing
//! each point is covered once at most. The plane is then triangulated
//! between the rims and a rectangle standing clear of the shape's box, with
//! the rims' edges as constraints, and the triangles the rims wind about
//! nothing (outside the region, and in its holes) are lifted to heights
//! taken linearly between the rims and the average plane, where the
//! rectangle lies. The sheet so extended is one graph over the rectangle,
//! so it cannot cross itself, and its only rim is the rectangle's.
//!
//! The sheet is thickened by the offset to each side with
//! [`thicken`](ogeom_mesh::thicken()), which refuses layers that run into
//! each other, built into a solid with
//! [`solid_from_mesh`](ogeom_algo::solid_from_mesh) and cut from the shape.
//! The offset layers' rim vertices stand within
//! [`OFFSET_FARTHEST`](ogeom_mesh::OFFSET_FARTHEST) times the offset of the
//! rectangle, which stands off the shape's box by the box's diagonal and
//! that reach, so the slab's walls meet nothing of the shape: the shape's
//! faces are cut only by the break and its extension, and each face the
//! slab misses comes through on its own surface.
//!
//! Every piece left lies wholly to one side of the slab. A piece is read by
//! its vertex farthest from the sheet along the normal: the sheet lies
//! over every point of the rectangle once, so a line along the normal
//! through the vertex meets it on one side only.

use ogeom_algo::{Built, History, MeshSolidOptions, solid_from_mesh, tight_bounds};
use ogeom_core::{Exact, OgeomResult, Predicates, Sign, Tolerances, ogeom_bail, ogeom_err};
use ogeom_math::{Point, Vector};
use ogeom_mesh::{Bvh, OFFSET_FARTHEST, Side, Thicken, thicken};
use ogeom_topo::{Model, Shape, ShapeType, Triangulation, explore_unique};
use spade::handles::{FixedDirectedEdgeHandle, FixedVertexHandle};
use spade::{ConstrainedDelaunayTriangulation, Triangulation as _};

/// Cut `shape` by the break surface `region` marks on `mesh`, thickened by
/// `offset` to each side, and keep the pieces on the `keep` side.
///
/// `region` holds one flag per triangle of `mesh`; the flagged triangles
/// are the break. They must be wound alike (each edge two of them share
/// run opposite ways in them) and face, as their corners wind
/// counter-clockwise, one way: seen along their average normal none may
/// face away or overlap another. [`Side::Above`] keeps the side the break's
/// faces face, [`Side::Below`] the other.
///
/// The break is extended past its rims along its average plane until it is
/// clear of the shape's box, its holes are closed, and the sheet is
/// thickened by `offset` to each side into a slab `2 * offset` thick, which
/// is cut from the shape: the face the kept pieces take along the break
/// stands `offset` off it, on the kept side. The slab's faces are the
/// thickened mesh's, merged where coplanar and recognized where they lie on
/// a cylinder, cone, sphere or torus. The shape's faces keep their
/// surfaces, those the break does not reach whole. The result is the kept solid, or a compound
/// of the kept solids where the break leaves several on that side. Its
/// history is the cut's, with the pieces on the other side deleted.
///
/// # Errors
///
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if
/// `region` does not hold one flag per triangle.
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `offset` is not a positive distance, `shape` holds no solid, the region
/// is empty, a triangle of it names a vertex twice or one that is not
/// there, two of its triangles run along an edge the same way, its
/// triangles' areas cancel so it has no average normal, a triangle faces
/// away from that normal, its rims cross or touch seen along it, or the
/// region does not separate the shape (no piece lies on one side of it,
/// the side kept included); and as [`thicken`](ogeom_mesh::thicken()),
/// [`solid_from_mesh`](ogeom_algo::solid_from_mesh) and [`cut`](crate::cut)
/// fail.
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the
/// thickened sheet does not build into one closed solid.
pub fn cut_by_mesh_region(
    model: &mut Model,
    shape: &Shape,
    mesh: &Triangulation,
    region: &[bool],
    offset: f64,
    keep: Side,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if region.len() != mesh.triangles.len() {
        ogeom_bail!(
            Dimension,
            "the region holds {} flags for {} triangles",
            region.len(),
            mesh.triangles.len()
        );
    }
    if !(offset.is_finite() && offset > 0.0) {
        ogeom_bail!(Construction, "an offset of {offset} is not a distance");
    }
    if explore_unique(model, shape, ShapeType::Solid)?.is_empty() {
        ogeom_bail!(Construction, "the shape to cut holds no solid");
    }
    let break_mesh = flagged(mesh, region)?;
    let bounds = tight_bounds(model, shape, tol)?;
    let (Some(low), Some(high)) = (bounds.low(), bounds.high()) else {
        ogeom_bail!(Construction, "the shape to cut has no extent");
    };
    let corners = [
        Point::new(low.x, low.y, low.z),
        Point::new(high.x, low.y, low.z),
        Point::new(low.x, high.y, low.z),
        Point::new(high.x, high.y, low.z),
        Point::new(low.x, low.y, high.z),
        Point::new(high.x, low.y, high.z),
        Point::new(low.x, high.y, high.z),
        Point::new(high.x, high.y, high.z),
    ];
    let margin = low.distance(high) + OFFSET_FARTHEST * offset;
    let sheet = extended(&break_mesh, &corners, margin, tol)?;

    let slab = thicken(&sheet.mesh, 2.0 * offset, Thicken::Both)?;
    let built = solid_from_mesh(model, &slab.mesh, &MeshSolidOptions::default(), tol)?;
    if !built.closed || explore_unique(model, &built.shape, ShapeType::Solid)?.len() != 1 {
        ogeom_bail!(
            NotDone,
            "the thickened break did not build into one closed solid"
        );
    }
    let cut = crate::cut(model, shape, &built.shape, tol)?;

    let bvh = Bvh::build(&sheet.mesh)?;
    let mut kept = Vec::new();
    let mut dropped = Vec::new();
    for piece in explore_unique(model, &cut.shape, ShapeType::Solid)? {
        if side_of(model, &piece, &bvh, sheet.normal)? == keep {
            kept.push(piece);
        } else {
            dropped.push(piece);
        }
    }
    if kept.is_empty() || dropped.is_empty() {
        ogeom_bail!(
            Construction,
            "the region does not separate the shape: {} pieces lie on the side kept \
             and {} on the other",
            kept.len(),
            dropped.len()
        );
    }

    let result = if let [one] = kept.as_slice() {
        one.clone()
    } else {
        model.add_compound(&kept)?
    };
    let mut step = History::new();
    for piece in &dropped {
        for kind in [
            ShapeType::Solid,
            ShapeType::Shell,
            ShapeType::Face,
            ShapeType::Wire,
            ShapeType::Edge,
            ShapeType::Vertex,
        ] {
            for sub in explore_unique(model, piece, kind)? {
                step.delete(&sub);
            }
        }
    }
    step.modify(&cut.shape, result.clone());
    Ok(Built::new(result, cut.history.followed_by(&step)))
}

/// The flagged triangles of `mesh`, their vertices renumbered from nought.
fn flagged(mesh: &Triangulation, region: &[bool]) -> OgeomResult<Triangulation> {
    let mut slot = vec![u32::MAX; mesh.positions.len()];
    let mut out = Triangulation::new();
    for (t, _) in mesh.triangles.iter().zip(region).filter(|&(_, &f)| f) {
        let mut corners = [0_u32; 3];
        for (k, &v) in t.iter().enumerate() {
            let Some(p) = mesh.positions.get(v as usize) else {
                ogeom_bail!(
                    Construction,
                    "a triangle of the region names vertex {v}, which is not there"
                );
            };
            if !p.is_finite() {
                ogeom_bail!(Construction, "vertex {v} of the region is not finite");
            }
            if slot[v as usize] == u32::MAX {
                #[allow(
                    clippy::cast_possible_truncation,
                    reason = "no more vertices than the mesh's indices name"
                )]
                let next = out.positions.len() as u32;
                slot[v as usize] = next;
                out.positions.push(*p);
            }
            corners[k] = slot[v as usize];
        }
        if corners[0] == corners[1] || corners[1] == corners[2] || corners[2] == corners[0] {
            ogeom_bail!(
                Construction,
                "a triangle of the region names a vertex twice"
            );
        }
        out.triangles.push(corners);
    }
    if out.triangles.is_empty() {
        ogeom_bail!(Construction, "the region is empty: no triangle is flagged");
    }
    Ok(out)
}

/// The break extended to a rectangle over its average plane.
struct Extended {
    /// The sheet: the break's triangles, then the ones closing its holes and
    /// reaching to the rectangle, all facing along `normal`.
    mesh: Triangulation,
    /// The break's average normal.
    normal: Vector,
}

/// The break extended over its average plane to a rectangle holding the
/// points `clear` and the break, grown by `margin` on each side.
fn extended(
    sheet: &Triangulation,
    clear: &[Point],
    margin: f64,
    tol: Tolerances,
) -> OgeomResult<Extended> {
    let at = |v: u32| sheet.positions[v as usize];

    // Each directed edge once, or two triangles run along it the same way.
    let mut edges: ogeom_core::FastMap<(u32, u32), ()> = ogeom_core::FastMap::default();
    for &[a, b, c] in &sheet.triangles {
        for (x, y) in [(a, b), (b, c), (c, a)] {
            if edges.insert((x, y), ()).is_some() {
                ogeom_bail!(
                    Construction,
                    "two triangles of the region run along one edge the same way: \
                     it is not wound alike or not 2-manifold"
                );
            }
        }
    }
    let mut rim: Vec<(u32, u32)> = edges
        .keys()
        .filter(|&&(x, y)| !edges.contains_key(&(y, x)))
        .copied()
        .collect();
    rim.sort_unstable();

    // The average normal and plane: the area vectors summed, through the
    // area-weighted centroid.
    let mut area = Vector::new(0.0, 0.0, 0.0);
    let mut moment = Vector::new(0.0, 0.0, 0.0);
    let mut total = 0.0;
    for &[a, b, c] in &sheet.triangles {
        let (pa, pb, pc) = (at(a), at(b), at(c));
        let twice = (pb - pa).cross(pc - pa);
        area += twice;
        let weight = twice.magnitude();
        total += weight;
        moment += (pa.to_vector() + pb.to_vector() + pc.to_vector()) * (weight / 3.0);
    }
    // A sum of areas that normalizes leaves `total` positive.
    let Ok(normal) = area.normalized(tol) else {
        ogeom_bail!(
            Construction,
            "the region's triangles' areas cancel: it has no average normal to extend along"
        );
    };
    let origin = Point::from_vector(moment / total);
    let seed = if normal.x.abs() <= normal.y.abs() && normal.x.abs() <= normal.z.abs() {
        Vector::new(1.0, 0.0, 0.0)
    } else if normal.y.abs() <= normal.z.abs() {
        Vector::new(0.0, 1.0, 0.0)
    } else {
        Vector::new(0.0, 0.0, 1.0)
    };
    let u = normal.cross(seed).normalized(tol)?;
    let v = normal.cross(u);
    let flat = |p: Point| {
        let d = p - origin;
        [d.dot(u), d.dot(v)]
    };
    let plane: Vec<[f64; 2]> = sheet.positions.iter().map(|&p| flat(p)).collect();

    for &[a, b, c] in &sheet.triangles {
        if Exact::orient2d(plane[a as usize], plane[b as usize], plane[c as usize])
            != Sign::Positive
        {
            ogeom_bail!(
                Construction,
                "a triangle of the region faces away from its average normal: \
                 the region folds over and cannot be extended without crossing itself"
            );
        }
    }

    // The rectangle, counter-clockwise.
    let mut lo = [f64::INFINITY; 2];
    let mut hi = [f64::NEG_INFINITY; 2];
    for q in plane.iter().copied().chain(clear.iter().map(|&p| flat(p))) {
        for k in 0..2 {
            lo[k] = lo[k].min(q[k]);
            hi[k] = hi[k].max(q[k]);
        }
    }
    let (x0, y0, x1, y1) = (
        lo[0] - margin,
        lo[1] - margin,
        hi[0] + margin,
        hi[1] + margin,
    );
    let frame = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];

    // The rims' vertices, then the rectangle's corners.
    let mut rim_vertices: Vec<u32> = rim.iter().flat_map(|&(x, y)| [x, y]).collect();
    rim_vertices.sort_unstable();
    rim_vertices.dedup();
    let mut coordinates: Vec<[f64; 2]> = rim_vertices.iter().map(|&w| plane[w as usize]).collect();
    coordinates.extend(frame);
    if coordinates
        .iter()
        .flatten()
        .any(|&x| spade::validate_coordinate(x).is_err())
    {
        ogeom_bail!(
            Construction,
            "the region's rim cannot be triangulated over its plane"
        );
    }
    let mut order: Vec<usize> = (0..coordinates.len()).collect();
    order.sort_by(|&a, &b| {
        let (p, q) = (coordinates[a], coordinates[b]);
        p[0].total_cmp(&q[0]).then(p[1].total_cmp(&q[1]))
    });
    if order
        .windows(2)
        .any(|pair| coordinates[pair[0]] == coordinates[pair[1]])
    {
        ogeom_bail!(
            Construction,
            "two vertices of the region's rims meet seen along its average normal"
        );
    }
    let slot = |w: u32| rim_vertices.binary_search(&w).unwrap_or(0);
    let points: Vec<spade::Point2<f64>> = coordinates
        .iter()
        .map(|c| spade::Point2::new(c[0], c[1]))
        .collect();
    let Ok(mut cdt) =
        ConstrainedDelaunayTriangulation::<spade::Point2<f64>>::bulk_load_cdt(points, Vec::new())
    else {
        ogeom_bail!(
            Construction,
            "the region's rim cannot be triangulated over its plane"
        );
    };
    let crossed = || {
        ogeom_err!(
            Construction,
            "the region's rims cross or touch seen along its average normal"
        )
    };
    let mut constrained: Vec<FixedDirectedEdgeHandle> = Vec::with_capacity(rim.len());
    for &(x, y) in &rim {
        let (from, to) = (
            FixedVertexHandle::from_index(slot(x)),
            FixedVertexHandle::from_index(slot(y)),
        );
        let made = cdt.try_add_constraint(from, to);
        let [e] = made[..] else {
            return Err(crossed());
        };
        let forward = cdt.directed_edge(e).from().fix() == from;
        constrained.push(if forward { e } else { e.rev() });
    }

    // How many times the rims wind about each face, walked from the
    // unbounded face: the face on a rim edge's left counts one more than
    // the face on its right.
    let mut weight = vec![0_i64; cdt.num_directed_edges()];
    for &e in &constrained {
        weight[e.index()] += 1;
        weight[e.rev().index()] -= 1;
    }
    let mut count: Vec<Option<i64>> = vec![None; cdt.num_all_faces()];
    let outer = cdt.outer_face().fix().index();
    count[outer] = Some(0);
    let mut sides = vec![[None; 3]; cdt.num_all_faces()];
    for face in cdt.inner_faces() {
        sides[face.fix().index()] = face.adjacent_edges().map(|e| Some(e.fix()));
    }
    let mut stack = vec![outer];
    let mut around: Vec<FixedDirectedEdgeHandle> = Vec::with_capacity(4);
    while let Some(f) = stack.pop() {
        let Some(here) = count[f] else { continue };
        around.clear();
        if f == outer {
            for e in cdt.convex_hull() {
                around.push(if e.face().is_outer() {
                    e.fix()
                } else {
                    e.rev().fix()
                });
            }
        } else {
            around.extend(sides[f].iter().flatten());
        }
        for &e in &around {
            let g = cdt.directed_edge(e).rev().face().fix().index();
            if count[g].is_none() {
                count[g] = Some(here - weight[e.index()]);
                stack.push(g);
            }
        }
    }
    if count.iter().any(|c| !matches!(c, Some(0 | 1))) {
        return Err(crossed());
    }

    // The sheet: the break, then the rectangle's corners on the average
    // plane, then the faces the rims wind about nothing.
    let mut out = sheet.clone();
    out.normals.clear();
    out.parameters.clear();
    let first_corner = u32::try_from(out.positions.len())
        .ok()
        .filter(|&n| n.checked_add(4).is_some())
        .ok_or_else(|| ogeom_err!(Construction, "the region has too many vertices to extend"))?;
    for [x, y] in frame {
        out.positions.push(origin + u * x + v * y);
    }
    let vertex = |k: usize| -> u32 {
        if k < rim_vertices.len() {
            rim_vertices[k]
        } else {
            #[allow(clippy::cast_possible_truncation, reason = "one of four corners")]
            let corner = (k - rim_vertices.len()) as u32;
            first_corner + corner
        }
    };
    for face in cdt.inner_faces() {
        if count[face.fix().index()] != Some(0) {
            continue;
        }
        let [a, b, c] = face.vertices().map(|h| h.fix().index());
        if Exact::orient2d(coordinates[a], coordinates[b], coordinates[c]) != Sign::Positive {
            ogeom_bail!(
                Construction,
                "the region's rim cannot be triangulated over its plane"
            );
        }
        out.triangles.push([vertex(a), vertex(b), vertex(c)]);
    }
    Ok(Extended { mesh: out, normal })
}

/// The side of the sheet `piece` lies on: read at its vertex farthest from
/// the sheet along the normal.
fn side_of(model: &Model, piece: &Shape, sheet: &Bvh, normal: Vector) -> OgeomResult<Side> {
    let mut best: Option<(f64, Side)> = None;
    for vertex in explore_unique(model, piece, ShapeType::Vertex)? {
        let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
            ogeom_bail!(Dangling, "a vertex of the cut is not in this model");
        };
        let p = vertex.transform(model.datums())?.apply(data.point);
        let up = sheet.ray(p, normal, f64::INFINITY, true);
        let down = sheet.ray(p, -normal, f64::INFINITY, true);
        let read = match (up, down) {
            (Some(hit), None) => (hit.distance, Side::Below),
            (None, Some(hit)) => (hit.distance, Side::Above),
            _ => continue,
        };
        if best.is_none_or(|(d, _)| read.0 > d) {
            best = Some(read);
        }
    }
    match best {
        Some((d, side)) if d > 0.0 => Ok(side),
        _ => ogeom_bail!(
            NotDone,
            "a piece of the cut has no vertex off the break to tell its side by"
        ),
    }
}
