//! Seams: edges between two faces on one surface facing one way, which a
//! boolean leaves where it joined two flats into one.
//!
//! A seam is no edge of the part's form: the surface runs on across it
//! unbroken and the material stays on the one side (mainstream CAD merges
//! such faces after a boolean, so the edge cannot even be picked). A
//! request naming one names nothing, and the engine reads a solid holding
//! seams with them dissolved, by the unification
//! [`ogeom_heal::unify_same_domain`] applies, its history carried back to
//! the solid as given.

use ogeom_algo::History;
use ogeom_core::{FastMap, OgeomError, OgeomResult, Tolerances};
use ogeom_math::Point;
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of};

use super::classify::{Planar, read_planar};
use super::{Id, id};
use crate::chamfer::Chamfer;

/// Whether `edge` is a seam of `solid`: bounded by exactly two faces of
/// one shell lying on one carrier facing one way
/// ([`ogeom_heal::same_domain`]). Two coplanar faces of different shells,
/// or facing opposite ways (the material changing side across the edge),
/// meet at a real edge.
pub(crate) fn is_seam(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let faces = ancestors_of(model, solid, edge, ShapeType::Face)?;
    let [a, b] = faces.as_slice() else {
        return Ok(false);
    };
    if a.node() == b.node() {
        return Ok(false);
    }
    let (around_a, around_b) = (
        ancestors_of(model, solid, a, ShapeType::Shell)?,
        ancestors_of(model, solid, b, ShapeType::Shell)?,
    );
    if !around_a
        .iter()
        .any(|s| around_b.iter().any(|t| t.node() == s.node()))
    {
        return Ok(false);
    }
    ogeom_heal::same_domain(model, a, b, tol)
}

/// A solid read for the engine with respect to its seams.
pub(crate) enum Dissolved {
    /// Nothing to dissolve, or not a solid bounded by planes: the engine
    /// reads the solid as given.
    Plain,
    /// Seams dissolved into a solid the engine cannot read: the engine
    /// declines.
    Outside,
    /// The solid with its seams dissolved.
    Into(Box<Seamless>),
}

/// A solid bounded by planes with its seams dissolved, and what each of
/// its elements became.
pub(crate) struct Seamless {
    /// The solid with its seams dissolved: each set of faces on one plane
    /// joined across seams one face, each run of edges on one line joined
    /// one edge.
    pub solid: Shape,
    input: Planar,
    unified: Planar,
    /// Each face of the input: the face of `unified` holding it.
    faces: Vec<usize>,
    /// Each edge of the input: the edge of `unified` running along it, or
    /// none for a seam.
    edges: FastMap<Id, Option<Id>>,
    /// Each vertex of the input: the vertex of `unified` there, or none
    /// where it dissolved into an edge or a face.
    vertices: FastMap<Id, Option<Id>>,
}

/// `solid` with its seams dissolved, read as the engine reads a solid.
///
/// # Errors
///
/// [`OgeomError::Cancelled`] when cancelled; any other failure to dissolve
/// declines.
pub(crate) fn dissolve(
    model: &mut Model,
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<Dissolved> {
    let Some(input) = read_planar(model, solid, tol)? else {
        return Ok(Dissolved::Plain);
    };
    let unified = match ogeom_heal::unify_same_domain(model, solid, tol) {
        Ok((_, 0)) => return Ok(Dissolved::Plain),
        Ok((unified, _)) => unified,
        Err(e @ OgeomError::Cancelled) => return Err(e),
        Err(_) => return Ok(Dissolved::Outside),
    };
    let joined = match ogeom_heal::merge_edges(model, &unified.shape, tol) {
        Ok((joined, _)) => joined.shape,
        Err(e @ OgeomError::Cancelled) => return Err(e),
        Err(_) => return Ok(Dissolved::Outside),
    };
    let Some(read) = read_planar(model, &joined, tol)? else {
        return Ok(Dissolved::Outside);
    };
    Ok(match Seamless::new(joined, input, read, tol) {
        Some(seamless) => Dissolved::Into(Box::new(seamless)),
        None => Dissolved::Outside,
    })
}

impl Seamless {
    /// Each element of `input` matched to the one of `unified` holding it,
    /// by where it lies: both solids are planes bounded by straight edges,
    /// so a vertex is matched by its point, an edge by the segment its ends
    /// lie on, a face by the face across one of its edges facing its way.
    /// `None` where a face of the input has no edge left to be matched by.
    fn new(solid: Shape, input: Planar, unified: Planar, tol: Tolerances) -> Option<Self> {
        let near = tol.confusion();
        let mut vertices: FastMap<Id, Option<Id>> = FastMap::default();
        for (v, p) in &input.points {
            let at = unified
                .points
                .iter()
                .find(|(_, q)| q.distance(*p) <= near)
                .map(|(w, _)| *w);
            vertices.insert(*v, at);
        }
        let on = |p: Point, (a, b): (Point, Point)| {
            let span = b - a;
            let length = span.magnitude();
            let t = (p - a).dot(span) / length;
            (-near..=length + near).contains(&t)
                && (p - (a + span * (t / length))).magnitude() <= near
        };
        let mut edges: FastMap<Id, Option<Id>> = FastMap::default();
        for ((a, b), e) in &input.edges {
            let ends = (input.points[a], input.points[b]);
            let along = unified.edges.iter().find(|((c, d), _)| {
                let segment = (unified.points[c], unified.points[d]);
                on(ends.0, segment) && on(ends.1, segment)
            });
            edges.insert(*e, along.map(|(_, f)| *f));
        }
        // Each unified edge by its vertex pair, with the faces walking it.
        let mut walking: FastMap<(Id, Id), Vec<usize>> = FastMap::default();
        for (gi, face) in unified.faces.iter().enumerate() {
            let n = face.corners.len();
            for k in 0..n {
                let (a, b) = (face.corners[k], face.corners[(k + 1) % n]);
                walking.entry((a.min(b), a.max(b))).or_default().push(gi);
            }
        }
        let pair_of: FastMap<Id, (Id, Id)> = unified.edges.iter().map(|(k, e)| (*e, *k)).collect();
        let mut faces = Vec::with_capacity(input.faces.len());
        for face in &input.faces {
            let n = face.corners.len();
            let image = (0..n).find_map(|k| {
                let (a, b) = (face.corners[k], face.corners[(k + 1) % n]);
                let e = input.edges.get(&(a.min(b), a.max(b)))?;
                let f = (*edges.get(e)?)?;
                walking
                    .get(pair_of.get(&f)?)?
                    .iter()
                    .copied()
                    .max_by(|x, y| {
                        let dx = unified.faces[*x].normal.dot(face.normal);
                        let dy = unified.faces[*y].normal.dot(face.normal);
                        dx.total_cmp(&dy)
                    })
            })?;
            faces.push(image);
        }
        Some(Self {
            solid,
            input,
            unified,
            faces,
            edges,
            vertices,
        })
    }

    /// `specs`, named on the solid as given, named on [`Seamless::solid`]:
    /// the pieces of one edge the seams split, asked alike, ask once.
    /// `None` where a spec names a seam or no element of the solid, or two
    /// pieces of one edge are asked differently.
    pub fn specs(&self, specs: &[(Shape, Chamfer)]) -> Option<Vec<(Shape, Chamfer)>> {
        let face = |named: &Shape| {
            let fi = self
                .input
                .faces
                .iter()
                .position(|f| f.shape.is_same(named))?;
            Some(self.unified.faces[self.faces[fi]].shape.clone())
        };
        let mut out: Vec<(Shape, Chamfer)> = Vec::with_capacity(specs.len());
        for (edge, spec) in specs {
            let f = (*self.edges.get(&id(edge))?)?;
            let edge = self.unified.shapes[&f].clone();
            let spec = match spec {
                Chamfer::Symmetric(d) => Chamfer::Symmetric(*d),
                Chamfer::Distances {
                    face: named,
                    on_face,
                    on_other,
                } => Chamfer::Distances {
                    face: face(named)?,
                    on_face: *on_face,
                    on_other: *on_other,
                },
                Chamfer::Angle {
                    face: named,
                    distance,
                    angle,
                } => Chamfer::Angle {
                    face: face(named)?,
                    distance: *distance,
                    angle: *angle,
                },
            };
            match out.iter().find(|(e, _)| e.is_same(&edge)) {
                Some((_, asked)) if alike(asked, &spec) => {}
                Some(_) => return None,
                None => out.push((edge, spec)),
            }
        }
        Some(out)
    }

    /// The history of a blend of [`Seamless::solid`], `made`, told of the
    /// solid as given: each face as the face holding it fared, each edge as
    /// the edge running along it fared, a seam deleted, each vertex as the
    /// vertex there fared, one dissolved into an edge or a face deleted.
    pub fn carry(&self, made: &History) -> History {
        let mut history = History::new();
        let mut pass = |from: &Shape, via: &Shape| {
            if made.is_deleted(via) {
                history.delete(from);
            }
            for shape in made.generated(via) {
                history.generate(from, shape.clone());
            }
            for shape in made.modified(via) {
                history.modify(from, shape.clone());
            }
        };
        for (fi, face) in self.input.faces.iter().enumerate() {
            pass(&face.shape, &self.unified.faces[self.faces[fi]].shape);
        }
        let mut ids: Vec<(Id, Option<Id>)> = self
            .edges
            .iter()
            .chain(&self.vertices)
            .map(|(k, v)| (*k, *v))
            .collect();
        ids.sort_unstable();
        let mut gone = Vec::new();
        for (from, via) in ids {
            match via {
                Some(via) => pass(&self.input.shapes[&from], &self.unified.shapes[&via]),
                None => gone.push(from),
            }
        }
        for from in gone {
            history.delete(&self.input.shapes[&from]);
        }
        history
    }
}

/// Whether two specs ask one chamfer: one spelling, the same numbers, the
/// same face named.
fn alike(a: &Chamfer, b: &Chamfer) -> bool {
    match (a, b) {
        (Chamfer::Symmetric(x), Chamfer::Symmetric(y)) => x == y,
        (
            Chamfer::Distances {
                face: f,
                on_face: a1,
                on_other: a2,
            },
            Chamfer::Distances {
                face: g,
                on_face: b1,
                on_other: b2,
            },
        ) => f.is_same(g) && a1 == b1 && a2 == b2,
        (
            Chamfer::Angle {
                face: f,
                distance: a1,
                angle: a2,
            },
            Chamfer::Angle {
                face: g,
                distance: b1,
                angle: b2,
            },
        ) => f.is_same(g) && a1 == b1 && a2 == b2,
        _ => false,
    }
}
