//! The mesh conversion's steps, between finding the regions and building
//! the solid.

use core::fmt;

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::{Axis, Cylinder, Frame, Plane, Point, Sphere, Torus, Vector};
use ogeom_topo::{Model, Triangulation};

use super::{
    Carrier, Curved, Found, MeshSolid, MeshSolidOptions, MeshSolidReport, align_one, axis_frame,
    from_to, hole_frame, plane_normals, plane_through, region_patch, sags_as_the_surface,
    slit_band, sphere_axis, swept_claim, tangent_blends, tangent_rounds, unit_normal,
};
use crate::recognize::{self, Canonical, recognize_curved, worst_deviation};
use crate::recognize_patch::Refused;

/// The regions of a mesh, found and not yet built.
///
/// [`MeshRegions::find`] welds, cleans and orients the mesh and gathers its
/// triangles into regions, each on the surface recognition found for it;
/// [`MeshRegions::build`] places the seams and corners and builds the
/// faces and the solid. [`solid_from_mesh`](super::solid_from_mesh) is the
/// one followed by the other. Between them a region can be merged with a
/// neighbour ([`MeshRegions::merge`]), split along a path of mesh vertices
/// ([`MeshRegions::split`]), or given a surface of a chosen kind
/// ([`MeshRegions::fit`]), and the build that follows goes through the
/// same seams, corners, faces and checks as an automatic one: a region
/// corrected to what recognition would have found builds what recognition
/// would have built.
///
/// Every surface a step puts on a region is verified as recognition's are:
/// every vertex of the region within the coplanar distance of it, and
/// every curved surface's triangles no farther off it in their middles
/// than their own sag. A step that cannot be taken is refused by name
/// ([`RegionRefusal`]), the regions unchanged.
///
/// A region a step left on a surface goes through the passes recognition
/// runs over its own regions: a round meeting two planes across smooth
/// edges is put on the cylinder tangent to both, and a round, ball or
/// torus among fillets on the surface its supports fix (see
/// [`MeshRegions::fit`] for the fits that keep what they were given). A
/// planar region of one or two triangles a step made is built as its own
/// face, never given to a curved neighbour as the build gives the facets
/// recognition leaves. A curved region the build cannot place falls back
/// to facets and is named in [`MeshSolidReport::fallbacks`], with why.
///
/// Indices name the mesh as the regions read it: [`MeshRegions::points`]
/// and [`MeshRegions::triangles`], welded, cleaned and oriented, not the
/// mesh handed to [`MeshRegions::find`] ([`MeshRegions::vertex_of`] maps
/// its vertices).
#[derive(Clone)]
pub struct MeshRegions {
    found: Found,
    options: MeshSolidOptions,
    tol: Tolerances,
    /// The regions a step made or changed, by index.
    touched: Vec<usize>,
}

/// Names one region of a [`MeshRegions`]. A merge keeps the first region's
/// name and retires the second's; a split keeps the name for one piece and
/// gives the other a new one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RegionId(pub(super) usize);

/// A curved region the build faceted, and why: see
/// [`MeshSolidReport::fallbacks`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RegionFallback {
    /// The region.
    pub region: RegionId,
    /// Why its face was not built on its surface.
    pub reason: FallbackReason,
}

/// Why a curved region was built as facets instead of on its surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum FallbackReason {
    /// Its seams or corners could not be placed on its surface and its
    /// neighbours' within their tolerance.
    BoundaryNotPlaced,
    /// Its face, built, reached past the triangles it replaces: its
    /// boundary closed on the wrong part of the surface.
    ReachesPast,
    /// Its face, built, faced into the material, against its triangles.
    TurnedIn,
    /// The body it is part of did not hold its mesh's volume, within what
    /// its curved faces may add over their facets.
    VolumeOff,
    /// One of its seams folded back in its chart, deeper than the seam's
    /// tolerance.
    FoldedSeam,
    /// Its face overlapped a face beside it.
    Overlaps,
    /// Its face's own mesh did not meet its neighbours': the solid would
    /// not tessellate closed.
    MeshesOpen,
    /// Building its face, or meshing or measuring it once built, failed.
    BuildFailed,
}

impl fmt::Display for FallbackReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::BoundaryNotPlaced => "its boundary could not be placed on its surface",
            Self::ReachesPast => "its face reached past its triangles",
            Self::TurnedIn => "its face faced into the material",
            Self::VolumeOff => "its body did not hold the mesh's volume",
            Self::FoldedSeam => "a seam of it folded in its chart",
            Self::Overlaps => "its face overlapped a neighbour",
            Self::MeshesOpen => "its face did not mesh closed with its neighbours",
            Self::BuildFailed => "building or meshing its face failed",
        })
    }
}

impl RegionId {
    /// The region's index among every region the mesh has had, retired
    /// ones included.
    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

/// One region as found: its triangles, the surface they are built on, and
/// the regions beside it.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshRegion {
    /// Its name.
    pub id: RegionId,
    /// Its triangles, indices into [`MeshRegions::triangles`], ascending.
    pub triangles: Vec<usize>,
    /// The surface its face is built on: a plane, a recognized curved
    /// surface, or a fitted sweep or patch ([`Canonical::Swept`]). `None`
    /// for a region with no surface, whose triangles are built as planar
    /// facets.
    pub surface: Option<Canonical>,
    /// The worst distance from any of its vertices to the surface,
    /// measured; `None` with no surface.
    pub deviation: Option<f64>,
    /// The regions sharing a mesh edge with it, ascending.
    pub neighbours: Vec<RegionId>,
}

/// A kind of surface to fit to a region with [`MeshRegions::fit`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SurfaceKind {
    /// A plane.
    Plane,
    /// A cylinder.
    Cylinder,
    /// A cone.
    Cone,
    /// A sphere.
    Sphere,
    /// A torus.
    Torus,
    /// A fitted B-spline patch, as the automatic conversion fits a smooth
    /// region nothing else fits ([`MeshSolidOptions::patches`]): over a
    /// region that is one disk bounded by one loop, verified at every
    /// vertex and across every triangle. It holds no axis or radius.
    Patch,
}

/// What a fit holds fixed.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FitConstraints {
    /// The surface's axis, as a line (its sense is free): a cylinder's, a
    /// cone's or a torus's. A plane and a sphere have none to fix.
    pub axis: Option<Axis>,
    /// A cylinder's or a sphere's radius, or a torus's tube (minor)
    /// radius. A plane has none, and a cone's changes along it.
    pub radius: Option<f64>,
}

/// Why a step on the regions was refused. The regions are unchanged.
#[derive(Debug, Clone, PartialEq)]
pub enum RegionRefusal {
    /// No region has this name, or it holds no triangle (a merge retired
    /// it).
    NoSuchRegion(RegionId),
    /// A region cannot be merged with itself.
    SameRegion(RegionId),
    /// The two regions share no mesh edge.
    NotAdjacent(RegionId, RegionId),
    /// A path needs at least two vertices.
    PathTooShort,
    /// Two consecutive vertices of the path are not the ends of an edge of
    /// the region's triangles.
    NotAnEdgeOfTheRegion(u32, u32),
    /// Cut along the path, the region falls into this many pieces, not two.
    DoesNotCut {
        /// How many pieces the path leaves.
        pieces: usize,
    },
    /// The kind of surface has no such thing to hold fixed.
    ConstraintDoesNotApply {
        /// The kind asked for.
        kind: SurfaceKind,
        /// The constraint, `"axis"` or `"radius"`.
        constraint: &'static str,
    },
    /// A fixed radius that is not a finite positive distance.
    InvalidRadius(f64),
    /// The region has too few vertices to determine the kind.
    TooFewVertices {
        /// The kind asked for.
        kind: SurfaceKind,
        /// How many it needs.
        needed: usize,
        /// How many the region has.
        found: usize,
    },
    /// No surface of the kind could be fitted to the vertices at all.
    NoFit(SurfaceKind),
    /// The fitted surface misses a vertex by more than the coplanar
    /// distance.
    DoesNotVerify {
        /// The kind asked for.
        kind: SurfaceKind,
        /// The worst distance from a vertex to it.
        deviation: f64,
        /// The coplanar distance it had to be within.
        distance: f64,
    },
    /// A triangle's middle stands off the fitted surface by more than its
    /// sag allows: its corners are on the surface and it spans off it.
    TriangleOffTheSurface {
        /// The kind asked for.
        kind: SurfaceKind,
        /// The triangle, an index into [`MeshRegions::triangles`].
        triangle: usize,
    },
    /// A patch needs a region that is one disk bounded by one loop; this
    /// one has a hole or touches itself.
    NotADisk,
    /// Too few of the region's vertices lie inside its boundary to hold a
    /// patch across it.
    TooNarrowForAPatch,
    /// No patch fitted to the region passed the verification: within the
    /// coplanar distance of every vertex and across every triangle, its
    /// normal agreeing with theirs.
    PatchDoesNotVerify,
}

impl fmt::Display for RegionRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoSuchRegion(id) => write!(f, "no region {} with triangles", id.0),
            Self::SameRegion(id) => write!(f, "region {} cannot be merged with itself", id.0),
            Self::NotAdjacent(a, b) => {
                write!(f, "regions {} and {} share no mesh edge", a.0, b.0)
            }
            Self::PathTooShort => write!(f, "a path needs at least two vertices"),
            Self::NotAnEdgeOfTheRegion(a, b) => {
                write!(
                    f,
                    "vertices {a} and {b} are not joined by an edge of the region"
                )
            }
            Self::DoesNotCut { pieces } => write!(
                f,
                "the path leaves the region in {pieces} piece(s), not two"
            ),
            Self::ConstraintDoesNotApply { kind, constraint } => {
                write!(f, "a {kind:?} has no {constraint} to fix")
            }
            Self::InvalidRadius(r) => write!(f, "the radius {r} is not a positive distance"),
            Self::TooFewVertices {
                kind,
                needed,
                found,
            } => write!(
                f,
                "a {kind:?} needs {needed} vertices and the region has {found}"
            ),
            Self::NoFit(kind) => write!(f, "no {kind:?} fits the region's vertices"),
            Self::DoesNotVerify {
                kind,
                deviation,
                distance,
            } => write!(
                f,
                "the fitted {kind:?} misses a vertex by {deviation:e}, past the distance {distance:e}"
            ),
            Self::TriangleOffTheSurface { kind, triangle } => write!(
                f,
                "triangle {triangle} spans off the fitted {kind:?} by more than its sag"
            ),
            Self::NotADisk => write!(f, "the region is not one disk, as a patch needs"),
            Self::TooNarrowForAPatch => {
                write!(f, "too few vertices lie inside the region to hold a patch")
            }
            Self::PatchDoesNotVerify => write!(f, "no patch fitted to the region verifies"),
        }
    }
}

impl std::error::Error for RegionRefusal {}

impl fmt::Debug for MeshRegions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MeshRegions")
            .field("triangles", &self.found.triangles.len())
            .field("regions", &self.ids().len())
            .field("distance", &self.found.flat)
            .finish_non_exhaustive()
    }
}

impl MeshRegions {
    /// Weld, clean and orient the mesh, and gather its triangles into
    /// regions as [`solid_from_mesh`](super::solid_from_mesh) does before it builds.
    ///
    /// # Errors
    ///
    /// As [`solid_from_mesh`](super::solid_from_mesh).
    pub fn find(
        mesh: &Triangulation,
        options: &MeshSolidOptions,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        Ok(Self {
            found: super::find(mesh, options, tol)?,
            options: *options,
            tol,
            touched: Vec::new(),
        })
    }

    /// Build the solid from the regions as they stand, through the same
    /// seams, corners, faces and checks as
    /// [`solid_from_mesh`](super::solid_from_mesh). A region with no
    /// surface is built as planar facets; a curved region whose boundary
    /// cannot be placed on its surface falls back to facets, and the
    /// report names it and why ([`MeshSolidReport::fallbacks`]): a region
    /// a step put on a surface and named there did not hold. A planar
    /// region of one or two triangles a step made stays a face of its own.
    ///
    /// # Errors
    ///
    /// As [`solid_from_mesh`](super::solid_from_mesh).
    pub fn build(&self, model: &mut Model) -> OgeomResult<MeshSolid> {
        let protected: Vec<bool> = if self.touched.is_empty() {
            Vec::new()
        } else {
            let mut by_region = vec![false; self.found.groups.carriers.len()];
            for &g in &self.touched {
                by_region[g] = true;
            }
            self.found
                .groups
                .of
                .iter()
                .map(|&g| by_region.get(g).copied().unwrap_or(false))
                .collect()
        };
        super::build(model, &self.found, &self.options, &protected, self.tol)
    }

    /// The mesh's vertices as the regions read them, welded.
    #[must_use]
    pub fn points(&self) -> &[Point] {
        &self.found.points
    }

    /// The mesh's triangles as the regions read them: without the ones of
    /// no area or repeated, wound consistently and outward.
    #[must_use]
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.found.triangles
    }

    /// The vertex of [`MeshRegions::points`] the input mesh's vertex
    /// `input` was welded into; `None` past the input's vertices.
    #[must_use]
    pub fn vertex_of(&self, input: usize) -> Option<u32> {
        self.found.remap.get(input).copied()
    }

    /// The coplanar distance the regions were found to, and every surface
    /// put on a region is verified at.
    #[must_use]
    pub fn distance(&self) -> f64 {
        self.found.flat
    }

    /// What finding the regions measured of the mesh: its welds, dropped
    /// triangles, windings and closure. The build's report starts from it.
    #[must_use]
    pub fn report(&self) -> &MeshSolidReport {
        &self.found.report
    }

    /// The region holding a triangle.
    #[must_use]
    pub fn region_of(&self, triangle: usize) -> Option<RegionId> {
        self.found
            .groups
            .of
            .get(triangle)
            .filter(|&&g| g < self.found.groups.carriers.len())
            .map(|&g| RegionId(g))
    }

    /// The names of the regions holding triangles, ascending.
    #[must_use]
    pub fn ids(&self) -> Vec<RegionId> {
        let mut held = vec![false; self.found.groups.carriers.len()];
        for &g in &self.found.groups.of {
            if let Some(h) = held.get_mut(g) {
                *h = true;
            }
        }
        held.iter()
            .enumerate()
            .filter(|(_, h)| **h)
            .map(|(g, _)| RegionId(g))
            .collect()
    }

    /// Every region holding triangles, in the order of their names.
    #[must_use]
    pub fn regions(&self) -> Vec<MeshRegion> {
        let members = self.members_all();
        members
            .into_iter()
            .enumerate()
            .filter(|(_, m)| !m.is_empty())
            .map(|(g, m)| self.describe(g, m))
            .collect()
    }

    /// One region, or `None` where no region of that name holds triangles.
    #[must_use]
    pub fn region(&self, id: RegionId) -> Option<MeshRegion> {
        let members = self.members(id.0);
        (!members.is_empty()).then(|| self.describe(id.0, members))
    }

    /// Merge region `b` into region `a`, which must share a mesh edge.
    ///
    /// The merged region takes the surface recognition finds for its
    /// triangles as a whole (a plane where they all lie within the
    /// distance of one, otherwise a cylinder, cone, sphere, torus or sweep
    /// where one verifies); where none does, the surface of the larger of
    /// the two where it holds the other's vertices (unless that is a patch,
    /// fitted to the larger alone); then, with
    /// [`MeshSolidOptions::patches`], a patch fitted to the whole where it
    /// is one disk and the patch verifies; and otherwise none, which
    /// builds as facets. `b` is retired and the merged region keeps the
    /// name `a`.
    ///
    /// # Errors
    ///
    /// [`RegionRefusal::NoSuchRegion`], [`RegionRefusal::SameRegion`] or
    /// [`RegionRefusal::NotAdjacent`].
    pub fn merge(&mut self, a: RegionId, b: RegionId) -> Result<RegionId, RegionRefusal> {
        let first = self.held(a)?;
        let second = self.held(b)?;
        if a == b {
            return Err(RegionRefusal::SameRegion(a));
        }
        let groups = &self.found.groups;
        let adjacent = first.iter().any(|&t| {
            (3 * t..3 * t + 3)
                .any(|h| self.found.adjacency.twin[h].is_some_and(|o| groups.of[o / 3] == b.0))
        });
        if !adjacent {
            return Err(RegionRefusal::NotAdjacent(a, b));
        }
        let mut union: Vec<usize> = first.iter().chain(&second).copied().collect();
        union.sort_unstable();
        let carrier = self.refitted(&union, || {
            let (large, small) = if first.len() >= second.len() {
                (a.0, &second)
            } else {
                (b.0, &first)
            };
            match &self.found.groups.carriers[large] {
                Carrier::Curved(c) if c.patch.is_some() => Carrier::Gone,
                carrier => self.kept(carrier, small, &union),
            }
        });
        for &t in &union {
            self.found.groups.of[t] = a.0;
        }
        self.found.groups.carriers[b.0] = Carrier::Gone;
        self.found.groups.carriers[a.0] = carrier;
        self.touched.retain(|&g| g != b.0);
        self.lay_out(a.0, false, true);
        Ok(a)
    }

    /// Split a region in two along a path of mesh vertices, consecutive
    /// ones joined by edges of the region's triangles: a path across the
    /// region from boundary to boundary, or a closed one (its last vertex
    /// its first) round part of it.
    ///
    /// Each piece takes the surface recognition finds for it, as a merged
    /// region does, and otherwise keeps the region's own (which holds it,
    /// as it held the whole); a piece of a region with no surface takes a
    /// patch as a merged region does. The piece holding the region's first
    /// triangle keeps the name; the other is named anew.
    ///
    /// # Errors
    ///
    /// [`RegionRefusal::NoSuchRegion`], [`RegionRefusal::PathTooShort`],
    /// [`RegionRefusal::NotAnEdgeOfTheRegion`] or
    /// [`RegionRefusal::DoesNotCut`] where the path leaves the region in
    /// one piece or in more than two.
    pub fn split(
        &mut self,
        id: RegionId,
        path: &[u32],
    ) -> Result<(RegionId, RegionId), RegionRefusal> {
        let members = self.held(id)?;
        if path.len() < 2 {
            return Err(RegionRefusal::PathTooShort);
        }
        let triangles = &self.found.triangles;
        let mut cut: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
        for pair in path.windows(2) {
            let (u, v) = (pair[0], pair[1]);
            let on = members.iter().any(|&t| {
                (3 * t..3 * t + 3).any(|h| {
                    let (a, b) = from_to(triangles, h);
                    (a, b) == (u, v) || (a, b) == (v, u)
                })
            });
            if !on {
                return Err(RegionRefusal::NotAnEdgeOfTheRegion(u, v));
            }
            cut.insert((u.min(v), u.max(v)));
        }
        let mut piece = vec![usize::MAX; triangles.len()];
        let mut pieces = 0_usize;
        for &start in &members {
            if piece[start] != usize::MAX {
                continue;
            }
            piece[start] = pieces;
            let mut stack = vec![start];
            while let Some(t) = stack.pop() {
                for h in 3 * t..3 * t + 3 {
                    let (a, b) = from_to(triangles, h);
                    if cut.contains(&(a.min(b), a.max(b))) {
                        continue;
                    }
                    if let Some(twin) = self.found.adjacency.twin[h] {
                        let o = twin / 3;
                        if self.found.groups.of[o] == id.0 && piece[o] == usize::MAX {
                            piece[o] = pieces;
                            stack.push(o);
                        }
                    }
                }
            }
            pieces += 1;
        }
        if pieces != 2 {
            return Err(RegionRefusal::DoesNotCut { pieces });
        }
        let (kept, other): (Vec<usize>, Vec<usize>) = members.iter().partition(|&&t| piece[t] == 0);
        let parent = self.found.groups.carriers[id.0].clone();
        let carriers =
            [&kept, &other].map(|part| self.refitted(part, || self.kept(&parent, part, part)));
        let [first, second] = carriers;
        let new = self.found.groups.carriers.len();
        self.found.groups.carriers.push(second);
        self.found.groups.carriers[id.0] = first;
        for &t in &other {
            self.found.groups.of[t] = new;
        }
        self.lay_out(id.0, false, true);
        self.lay_out(new, false, true);
        Ok((id, RegionId(new)))
    }

    /// Fit a surface of the kind to the region, holding what `constraints`
    /// fix, and put the region on it where it verifies: every vertex of
    /// the region within the coplanar distance, and every triangle's middle
    /// no farther off it than its sag. Returns the worst distance from a
    /// vertex to the surface.
    ///
    /// The fit is least squares on the vertices' distances to the surface.
    /// With an axis fixed it is the surface about that line (a cylinder's
    /// radius the vertices' mean distance from it, a torus's spine from
    /// the circle its profile traces); with a radius fixed the rest moves
    /// about it. An axis fixed is kept as given; otherwise the axis is set
    /// square to a plane of the solid it all but is square to, and onto a
    /// coaxial region's, where the surface still verifies there, as
    /// recognition does. A fit holding nothing then goes through the
    /// tangent passes as recognition's regions do (a round between two
    /// planes put on the cylinder tangent to both, a fillet on the surface
    /// its supports fix), where the surface they give still verifies; a fit
    /// holding an axis or a radius keeps the surface it was given.
    ///
    /// A [`SurfaceKind::Patch`] is fitted over the region's triangles by
    /// the automatic conversion's own patch fit and verification, whatever
    /// [`MeshSolidOptions::patches`] says, and is kept as fitted.
    ///
    /// # Errors
    ///
    /// [`RegionRefusal::NoSuchRegion`],
    /// [`RegionRefusal::ConstraintDoesNotApply`],
    /// [`RegionRefusal::InvalidRadius`], [`RegionRefusal::TooFewVertices`],
    /// [`RegionRefusal::NoFit`], [`RegionRefusal::DoesNotVerify`] or
    /// [`RegionRefusal::TriangleOffTheSurface`]; for a patch,
    /// [`RegionRefusal::NotADisk`], [`RegionRefusal::TooNarrowForAPatch`]
    /// or [`RegionRefusal::PatchDoesNotVerify`].
    pub fn fit(
        &mut self,
        id: RegionId,
        kind: SurfaceKind,
        constraints: &FitConstraints,
    ) -> Result<f64, RegionRefusal> {
        let members = self.held(id)?;
        applies(kind, constraints, self.tol)?;
        if kind == SurfaceKind::Patch {
            let claim = self.patch(&members).map_err(|refused| match refused {
                Refused::NotDisk => RegionRefusal::NotADisk,
                Refused::Narrow => RegionRefusal::TooNarrowForAPatch,
                Refused::Unverified => RegionRefusal::PatchDoesNotVerify,
            })?;
            let deviation = claim.deviation;
            self.found.groups.carriers[id.0] = Carrier::Curved(claim);
            self.lay_out(id.0, false, false);
            return Ok(deviation);
        }
        let vertices = self.vertices(&members);
        let (pts, nrm) = self.samples(&vertices, &members);
        let shape = fitted(kind, constraints, &pts, &nrm, self.found.flat, self.tol)?;
        let flat = self.found.flat;
        let deviation = worst_deviation(&shape, &pts);
        if deviation.is_nan() || deviation > flat {
            return Err(RegionRefusal::DoesNotVerify {
                kind,
                deviation,
                distance: flat,
            });
        }
        if kind != SurfaceKind::Plane
            && let Some(&t) = members.iter().find(|&&t| !self.sags(&shape, t))
        {
            return Err(RegionRefusal::TriangleOffTheSurface { kind, triangle: t });
        }
        let carrier = match shape {
            Canonical::Plane(plane) => Carrier::Plane(self.outward(plane, &members)),
            shape => Carrier::Curved(curved(shape, deviation, vertices)),
        };
        self.found.groups.carriers[id.0] = carrier;
        let free = constraints.axis.is_none() && constraints.radius.is_none();
        self.lay_out(id.0, constraints.axis.is_some(), free);
        // Set onto a shared axis, the surface is measured again there.
        Ok(match &self.found.groups.carriers[id.0] {
            Carrier::Curved(c) => c.deviation,
            _ => deviation,
        })
    }

    /// Put `surface` on a region as given, unverified, keeping how its face
    /// is laid out: a plane on a planar region, or a curved surface of the
    /// kind the region is on. What the build makes of a surface its
    /// vertices do not lie on is what its fallbacks catch; this is how a
    /// test hands it one. Not part of the supported interface.
    ///
    /// # Errors
    ///
    /// [`RegionRefusal::NoSuchRegion`], or [`RegionRefusal::NoFit`] where
    /// the surface is not of the region's kind.
    #[doc(hidden)]
    pub fn put_surface_unverified(
        &mut self,
        id: RegionId,
        surface: Canonical,
    ) -> Result<(), RegionRefusal> {
        self.held(id)?;
        let carrier = &mut self.found.groups.carriers[id.0];
        match (carrier, surface) {
            (Carrier::Plane(plane), Canonical::Plane(given)) => *plane = given,
            (Carrier::Curved(curved), given)
                if core::mem::discriminant(&curved.shape) == core::mem::discriminant(&given) =>
            {
                curved.shape = given;
            }
            (_, given) => {
                return Err(RegionRefusal::NoFit(match given {
                    Canonical::Plane(_) => SurfaceKind::Plane,
                    Canonical::Cylinder(_) => SurfaceKind::Cylinder,
                    Canonical::Cone(_) => SurfaceKind::Cone,
                    Canonical::Sphere(_) => SurfaceKind::Sphere,
                    Canonical::Torus(_) => SurfaceKind::Torus,
                    Canonical::Swept(_) => SurfaceKind::Patch,
                }));
            }
        }
        Ok(())
    }

    /// Each region's triangles, by region.
    fn members_all(&self) -> Vec<Vec<usize>> {
        let mut members = vec![Vec::new(); self.found.groups.carriers.len()];
        for (t, &g) in self.found.groups.of.iter().enumerate() {
            if let Some(m) = members.get_mut(g) {
                m.push(t);
            }
        }
        members
    }

    fn members(&self, g: usize) -> Vec<usize> {
        (0..self.found.triangles.len())
            .filter(|&t| self.found.groups.of[t] == g)
            .collect()
    }

    /// A region's triangles, refused where it holds none.
    fn held(&self, id: RegionId) -> Result<Vec<usize>, RegionRefusal> {
        let members = self.members(id.0);
        if members.is_empty() {
            return Err(RegionRefusal::NoSuchRegion(id));
        }
        Ok(members)
    }

    fn vertices(&self, members: &[usize]) -> Vec<u32> {
        let mut vertices: Vec<u32> = members
            .iter()
            .flat_map(|&t| self.found.triangles[t])
            .collect();
        vertices.sort_unstable();
        vertices.dedup();
        vertices
    }

    fn points_of(&self, vertices: &[u32]) -> Vec<Point> {
        vertices
            .iter()
            .map(|&v| self.found.points[v as usize])
            .collect()
    }

    /// The vertices with normals averaged over the region's triangles.
    fn samples(&self, vertices: &[u32], members: &[usize]) -> (Vec<Point>, Vec<Vector>) {
        let mut sum: std::collections::HashMap<u32, Vector> = std::collections::HashMap::new();
        for &t in members {
            let n = unit_normal(&self.found.points, self.found.triangles[t]);
            for &v in &self.found.triangles[t] {
                *sum.entry(v).or_insert(Vector::ZERO) += n;
            }
        }
        let nrm = vertices
            .iter()
            .map(|v| {
                let s = sum.get(v).copied().unwrap_or(Vector::Z);
                let m = s.magnitude();
                if m > 0.0 { s / m } else { Vector::Z }
            })
            .collect();
        (self.points_of(vertices), nrm)
    }

    /// Whether a triangle sags as the surface does.
    fn sags(&self, shape: &Canonical, t: usize) -> bool {
        let corners = self.found.triangles[t].map(|v| self.found.points[v as usize]);
        sags_as_the_surface(shape, corners, self.found.flat)
    }

    /// A plane turned to face the way the region's triangles do.
    fn outward(&self, plane: Plane, members: &[usize]) -> Plane {
        let facing = members.iter().fold(Vector::ZERO, |s, &t| {
            s + unit_normal(&self.found.points, self.found.triangles[t])
        });
        let frame = plane.frame();
        if facing.dot(frame.z().vector()) >= 0.0 {
            plane
        } else {
            Plane::new(Frame::about(frame.origin(), -frame.z()))
        }
    }

    /// The surface recognition finds for a set of triangles taken whole: a
    /// plane, a canonical curved surface, or (with the sweeps on) a sweep.
    fn recognized(&self, members: &[usize]) -> Option<Carrier> {
        let flat = self.found.flat;
        let vertices = self.vertices(members);
        let (pts, nrm) = self.samples(&vertices, members);
        if pts.len() >= 3 && recognize::is_flat(&pts, flat, self.tol) {
            let (centre, normal) = plane_through(&pts, self.tol)?;
            let plane = Plane::new(Frame::about(centre, normal));
            return Some(Carrier::Plane(self.outward(plane, members)));
        }
        if !(self.options.recognize && self.options.merge_coplanar) {
            return None;
        }
        let stride = members.len().div_ceil(100).max(1);
        let chords: Vec<(Point, Point)> = members
            .iter()
            .step_by(stride)
            .flat_map(|&t| {
                let [a, b, c] = self.found.triangles[t].map(|v| self.found.points[v as usize]);
                [(a, b), (b, c), (c, a)]
            })
            .collect();
        if let Some(found) = recognize_curved(&pts, &nrm, &chords, flat, self.tol)
            && members.iter().all(|&t| self.sags(&found.surface, t))
        {
            return Some(Carrier::Curved(curved(
                found.surface,
                found.deviation,
                vertices,
            )));
        }
        if self.options.sweeps {
            let normals: Vec<Vector> = self
                .found
                .triangles
                .iter()
                .map(|t| unit_normal(&self.found.points, *t))
                .collect();
            return swept_claim(
                &self.found.points,
                &self.found.triangles,
                &normals,
                members,
                flat,
                self.tol,
            )
            .map(Carrier::Curved);
        }
        None
    }

    /// The surface recognition finds for a set of triangles taken whole;
    /// else the one `kept` gives; else, with the patches on, a patch over
    /// them; else none.
    fn refitted(&self, members: &[usize], kept: impl FnOnce() -> Carrier) -> Carrier {
        if let Some(carrier) = self.recognized(members) {
            return carrier;
        }
        let carrier = kept();
        if !matches!(carrier, Carrier::Gone) {
            return carrier;
        }
        if self.options.recognize && self.options.merge_coplanar && self.options.patches {
            return self.patch(members).map_or(Carrier::Gone, Carrier::Curved);
        }
        Carrier::Gone
    }

    /// A patch fitted to a set of triangles and verified, as the automatic
    /// conversion fits one to a smooth region, the faces across its
    /// boundary read from the regions as they stand.
    fn patch(&self, members: &[usize]) -> Result<Curved, Refused> {
        let f = &self.found;
        region_patch(
            &f.points,
            &f.triangles,
            &f.adjacency,
            &f.groups,
            members,
            f.flat,
            self.tol,
        )
    }

    /// A carrier kept for a region's new triangles where it holds the
    /// vertices of `added` within the distance, measured over `all`; gone
    /// otherwise.
    fn kept(&self, carrier: &Carrier, added: &[usize], all: &[usize]) -> Carrier {
        let flat = self.found.flat;
        let pts = self.points_of(&self.vertices(added));
        match carrier {
            Carrier::Plane(plane) => {
                if pts
                    .iter()
                    .all(|p| plane.signed_distance_to(*p).abs() <= flat)
                {
                    Carrier::Plane(*plane)
                } else {
                    Carrier::Gone
                }
            }
            Carrier::Curved(c) => {
                if worst_deviation(&c.shape, &pts) > flat {
                    return Carrier::Gone;
                }
                let vertices = self.vertices(all);
                let deviation = worst_deviation(&c.shape, &self.points_of(&vertices));
                Carrier::Curved(Curved {
                    deviation,
                    vertices,
                    ..c.clone()
                })
            }
            Carrier::Gone => Carrier::Gone,
        }
    }

    /// Lay a curved region out on its surface as recognition lays its own:
    /// a sphere turned about the normal its boundary circles share, with
    /// `tangent` a round between two planes put on the cylinder tangent to
    /// both, the axis squared to a plane and set onto a coaxial region's
    /// unless it was fixed, with `tangent` a fillet put on the surface its
    /// supports fix, the chart's branch and whether the region wraps, and a
    /// closed surface's seams set clear of its holes. The region is marked
    /// as one a step touched. A patch takes no tangent pass: it is what is
    /// left where no surface of theirs fits.
    fn lay_out(&mut self, g: usize, axis_fixed: bool, tangent: bool) {
        if !self.touched.contains(&g) {
            self.touched.push(g);
        }
        let cos_crease = self.options.crease.cos();
        let f = &mut self.found;
        let tol = self.tol;
        sphere_axis(
            &f.points,
            &f.triangles,
            &f.adjacency,
            &mut f.groups,
            g,
            f.flat,
            tol,
        );
        let tangent =
            tangent && matches!(&f.groups.carriers[g], Carrier::Curved(c) if c.patch.is_none());
        if tangent {
            let planes = f.groups.clone();
            tangent_rounds(
                &f.points,
                &f.triangles,
                &f.adjacency,
                &mut f.groups,
                &planes,
                cos_crease,
                f.flat,
                Some(g),
                tol,
            );
        }
        let Carrier::Curved(mut curved) = f.groups.carriers[g].clone() else {
            return;
        };
        let (mut leaders, normals) = if axis_fixed {
            (Vec::new(), Vec::new())
        } else {
            let mut others: Vec<(usize, usize, Frame)> = f
                .groups
                .carriers
                .iter()
                .enumerate()
                .filter(|&(o, _)| o != g)
                .filter_map(|(_, c)| match c {
                    Carrier::Curved(c) => {
                        let rank = match c.shape {
                            Canonical::Cylinder(_) => 0,
                            Canonical::Cone(_) => 1,
                            Canonical::Torus(_) => 2,
                            _ => return None,
                        };
                        Some((rank, usize::MAX - c.vertices.len(), axis_frame(&c.shape)?))
                    }
                    _ => None,
                })
                .collect();
            others.sort_by_key(|&(rank, size, _)| (rank, size));
            (
                others.into_iter().map(|(_, _, frame)| frame).collect(),
                plane_normals(&f.groups),
            )
        };
        align_one(&f.points, &mut curved, &mut leaders, &normals, f.flat, tol);
        f.groups.carriers[g] = Carrier::Curved(curved);
        if tangent {
            let planes = f.groups.clone();
            tangent_blends(
                &f.points,
                &f.triangles,
                &f.adjacency,
                &mut f.groups,
                &planes,
                cos_crease,
                f.flat,
                Some(g),
                tol,
            );
        }
        hole_frame(&f.points, &f.triangles, &f.adjacency, &mut f.groups, g, tol);
        slit_band(&f.points, &f.triangles, &f.adjacency, &mut f.groups, g, tol);
    }

    fn describe(&self, g: usize, triangles: Vec<usize>) -> MeshRegion {
        let groups = &self.found.groups;
        let mut neighbours: Vec<RegionId> = triangles
            .iter()
            .flat_map(|&t| {
                (3 * t..3 * t + 3).filter_map(|h| self.found.adjacency.twin[h].map(|o| o / 3))
            })
            .map(|o| groups.of[o])
            .filter(|&o| o != g && o < groups.carriers.len())
            .map(RegionId)
            .collect();
        neighbours.sort_unstable();
        neighbours.dedup();
        let (surface, deviation) = match &groups.carriers[g] {
            Carrier::Plane(plane) => {
                let pts = self.points_of(&self.vertices(&triangles));
                let deviation = pts
                    .iter()
                    .map(|p| plane.signed_distance_to(*p).abs())
                    .fold(0.0_f64, f64::max);
                (Some(Canonical::Plane(*plane)), Some(deviation))
            }
            Carrier::Curved(c) => (Some(c.shape.clone()), Some(c.deviation)),
            Carrier::Gone => (None, None),
        };
        MeshRegion {
            id: RegionId(g),
            triangles,
            surface,
            deviation,
            neighbours,
        }
    }
}

/// A region on a curved surface, not yet laid out.
fn curved(shape: Canonical, deviation: f64, vertices: Vec<u32>) -> Curved {
    Curved {
        shape,
        deviation,
        fitted: deviation,
        centre: (0.0, 0.0),
        wraps: false,
        wraps_v: false,
        fixed: false,
        vertices,
        patch: None,
    }
}

/// The fewest vertices each kind is fitted to: half as many again as its
/// unknowns, as recognition asks.
/// A patch counts the vertices inside its boundary instead.
const fn floor(kind: SurfaceKind) -> usize {
    match kind {
        SurfaceKind::Plane => 3,
        SurfaceKind::Sphere => 6,
        SurfaceKind::Cylinder => 8,
        SurfaceKind::Cone => 9,
        SurfaceKind::Torus => 11,
        SurfaceKind::Patch => 0,
    }
}

/// Whether the kind has what the constraints fix, and a fixed radius is a
/// distance.
fn applies(
    kind: SurfaceKind,
    constraints: &FitConstraints,
    tol: Tolerances,
) -> Result<(), RegionRefusal> {
    let refuse = |constraint| Err(RegionRefusal::ConstraintDoesNotApply { kind, constraint });
    match kind {
        SurfaceKind::Plane | SurfaceKind::Sphere | SurfaceKind::Patch
            if constraints.axis.is_some() =>
        {
            return refuse("axis");
        }
        SurfaceKind::Plane | SurfaceKind::Cone | SurfaceKind::Patch
            if constraints.radius.is_some() =>
        {
            return refuse("radius");
        }
        _ => {}
    }
    if let Some(r) = constraints.radius
        && !(r.is_finite() && r > tol.confusion())
    {
        return Err(RegionRefusal::InvalidRadius(r));
    }
    Ok(())
}

/// A canonical surface of the kind fitted to the points, holding the
/// constraints, which [`applies`] has checked.
fn fitted(
    kind: SurfaceKind,
    constraints: &FitConstraints,
    pts: &[Point],
    nrm: &[Vector],
    flat: f64,
    tol: Tolerances,
) -> Result<Canonical, RegionRefusal> {
    let needed = floor(kind);
    if pts.len() < needed {
        return Err(RegionRefusal::TooFewVertices {
            kind,
            needed,
            found: pts.len(),
        });
    }
    let none = RegionRefusal::NoFit(kind);
    let free = |seed: Option<Canonical>| -> Result<Canonical, RegionRefusal> {
        let seed = seed.ok_or_else(|| none.clone())?;
        Ok(recognize::refine(seed.clone(), pts, flat, tol).unwrap_or(seed))
    };
    let (axis, radius) = (constraints.axis, constraints.radius);
    let centroid = Point::from_vector(
        pts.iter().fold(Vector::ZERO, |s, p| s + p.to_vector()) / count(pts.len()),
    );
    match kind {
        // Fitted over the triangles, not the points alone.
        SurfaceKind::Patch => Err(none),
        SurfaceKind::Plane => {
            let (centre, normal) = plane_through(pts, tol).ok_or(none)?;
            Ok(Canonical::Plane(Plane::new(Frame::about(centre, normal))))
        }
        SurfaceKind::Sphere => {
            let seed = free(recognize::fit_sphere(pts, tol))?;
            let Some(r) = radius else {
                return Ok(seed);
            };
            let Canonical::Sphere(s) = seed else {
                return Err(none);
            };
            let c = s.centre();
            let y = solve(vec![c.x, c.y, c.z], pts, &|y, p| {
                p.distance(Point::new(y[0], y[1], y[2])) - r
            });
            Sphere::centred(Point::new(y[0], y[1], y[2]), r, tol)
                .map(Canonical::Sphere)
                .map_err(|_| none)
        }
        SurfaceKind::Cylinder => match (axis, radius) {
            (Some(axis), r) => {
                let origin = on_line(axis, centroid);
                let r = r.unwrap_or_else(|| {
                    pts.iter()
                        .map(|p| distance_from_line(axis, *p))
                        .sum::<f64>()
                        / count(pts.len())
                });
                Cylinder::new(Frame::about(origin, axis.direction), r, tol)
                    .map(Canonical::Cylinder)
                    .map_err(|_| none)
            }
            (None, Some(r)) => {
                let seed = free(recognize::fit_cylinder(pts, nrm, tol))?;
                let Canonical::Cylinder(c) = seed else {
                    return Err(none);
                };
                let (o, d) = (c.frame().origin(), c.frame().z().vector());
                let like = Canonical::Cylinder(c);
                let y = solve(vec![o.x, o.y, o.z, d.x, d.y, d.z], pts, &|y, p| {
                    recognize::residual(&like, &[y[0], y[1], y[2], y[3], y[4], y[5], r], p)
                });
                recognize::rebuild(&like, &[y[0], y[1], y[2], y[3], y[4], y[5], r], tol).ok_or(none)
            }
            (None, None) => free(recognize::fit_cylinder(pts, nrm, tol)),
        },
        SurfaceKind::Cone => {
            let Some(axis) = axis else {
                return free(recognize::fit_cone(pts, nrm, tol));
            };
            let Some(Canonical::Cone(seed)) =
                recognize::ruled_about(pts, axis.location, axis.direction, 0.0, tol)
            else {
                return Err(none);
            };
            let (o, d) = (seed.frame().origin(), seed.frame().z().vector());
            let like = Canonical::Cone(seed);
            let x = |y: &[f64]| [o.x, o.y, o.z, d.x, d.y, d.z, y[0], y[1]];
            let y = solve(
                vec![seed.reference_radius(), seed.half_angle()],
                pts,
                &|y, p| recognize::residual(&like, &x(y), p),
            );
            recognize::rebuild(&like, &x(&y), tol).ok_or(none)
        }
        SurfaceKind::Torus => match axis {
            Some(axis) => {
                let d = axis.direction.vector();
                let base = on_line(axis, centroid);
                let (major, height, minor) = profile_circle(pts, base, axis).ok_or(none.clone())?;
                let like = Canonical::Torus(
                    Torus::new(
                        Frame::about(base + d * height, axis.direction),
                        major,
                        radius.unwrap_or(minor),
                        tol,
                    )
                    .map_err(|_| none.clone())?,
                );
                let x = |y: &[f64]| {
                    let o = base + d * y[0];
                    [
                        o.x,
                        o.y,
                        o.z,
                        d.x,
                        d.y,
                        d.z,
                        y[1],
                        radius.unwrap_or_else(|| y[2]),
                    ]
                };
                let y0 = match radius {
                    Some(_) => vec![height, major],
                    None => vec![height, major, minor],
                };
                let y = solve(y0, pts, &|y, p| recognize::residual(&like, &x(y), p));
                recognize::rebuild(&like, &x(&y), tol).ok_or(none)
            }
            None => {
                let seed = free(recognize::fit_torus(pts, nrm, tol))?;
                let Some(r) = radius else {
                    return Ok(seed);
                };
                let Canonical::Torus(t) = seed else {
                    return Err(none);
                };
                let (o, d) = (t.frame().origin(), t.frame().z().vector());
                let like = Canonical::Torus(t);
                let x = |y: &[f64]| [y[0], y[1], y[2], y[3], y[4], y[5], y[6], r];
                let y = solve(
                    vec![o.x, o.y, o.z, d.x, d.y, d.z, t.major_radius()],
                    pts,
                    &|y, p| recognize::residual(&like, &x(y), p),
                );
                recognize::rebuild(&like, &x(&y), tol).ok_or(none)
            }
        },
    }
}

#[allow(
    clippy::cast_precision_loss,
    reason = "vertex counts are far below 2^52"
)]
fn count(n: usize) -> f64 {
    n.max(1) as f64
}

/// The point of the axis nearest `p`.
fn on_line(axis: Axis, p: Point) -> Point {
    let d = axis.direction.vector();
    axis.location + d * (p - axis.location).dot(d)
}

fn distance_from_line(axis: Axis, p: Point) -> f64 {
    let d = axis.direction.vector();
    let w = p - axis.location;
    (w - d * w.dot(d)).magnitude()
}

/// The circle the points trace in the half-plane through the axis (their
/// distance from it against their height along it from `base`): its
/// centre's distance and height, and its radius. A torus about the axis
/// is that circle swept round it.
fn profile_circle(pts: &[Point], base: Point, axis: Axis) -> Option<(f64, f64, f64)> {
    let d = axis.direction.vector();
    let mut m = nalgebra::Matrix3::<f64>::zeros();
    let mut b = nalgebra::Vector3::<f64>::zeros();
    for p in pts {
        let w = *p - base;
        let h = w.dot(d);
        let rho = (w - d * h).magnitude();
        let row = nalgebra::Vector3::new(2.0 * rho, 2.0 * h, 1.0);
        m += row * row.transpose();
        b += row * rho.mul_add(rho, h * h);
    }
    let x = m.lu().solve(&b)?;
    let r2 = x[2] + x[0].mul_add(x[0], x[1] * x[1]);
    (r2 > 0.0).then(|| (x[0], x[1], r2.sqrt()))
}

/// Least squares on a residual over the points, from `y`: Gauss-Newton
/// steps through the Jacobian's singular value decomposition, the
/// directions the points do not determine left where they are, each step
/// halved until it lowers the cost. Fitted on an even subsample of a large
/// set, as recognition fits.
fn solve(mut y: Vec<f64>, points: &[Point], residual: &dyn Fn(&[f64], Point) -> f64) -> Vec<f64> {
    const FIT_SAMPLES: usize = 1500;
    let stride = points.len().div_ceil(FIT_SAMPLES).max(1);
    let points: Vec<Point> = points.iter().step_by(stride).copied().collect();
    let n = y.len();
    let cost = |y: &[f64]| -> f64 { points.iter().map(|p| residual(y, *p).powi(2)).sum() };
    let scale = points
        .iter()
        .map(|p| p.distance(points[0]))
        .fold(f64::MIN_POSITIVE, f64::max);
    let mut current = cost(&y);
    for _ in 0..60 {
        let steps: Vec<f64> = y.iter().map(|v| 1e-7 * v.abs().max(scale)).collect();
        let mut jacobian = nalgebra::DMatrix::<f64>::zeros(points.len(), n);
        let mut r = nalgebra::DVector::<f64>::zeros(points.len());
        for (i, p) in points.iter().enumerate() {
            r[i] = residual(&y, *p);
            for k in 0..n {
                let mut up = y.clone();
                let mut down = y.clone();
                up[k] += steps[k];
                down[k] -= steps[k];
                jacobian[(i, k)] = (residual(&up, *p) - residual(&down, *p)) / (2.0 * steps[k]);
            }
        }
        let svd = jacobian.svd(true, true);
        let largest = svd.singular_values.max();
        if largest.is_nan() || largest <= 0.0 {
            break;
        }
        let Ok(step) = svd.solve(&(-r), largest * 1e-8) else {
            break;
        };
        let mut scale_step = 1.0;
        let mut gained = 0.0;
        for _ in 0..20 {
            let trial: Vec<f64> = y
                .iter()
                .zip(step.iter())
                .map(|(a, b)| b.mul_add(scale_step, *a))
                .collect();
            let c = cost(&trial);
            if c.is_finite() && c < current {
                gained = current - c;
                y = trial;
                current = c;
                break;
            }
            scale_step *= 0.5;
        }
        if gained <= current * 1e-12 {
            break;
        }
    }
    y
}
