//! Gathering: each operand's faces, boundary edges and the overlaps between
//! them, read once before the filler runs.

use super::*;

// --- gathering ---------------------------------------------------------------

/// An edge occurrence's identity: its node and where it is placed.
///
/// One node placed twice is two edges in space: a prism's far cap is its
/// near cap moved, the same nodes under a displacement. Paves shared by
/// node alone would split the near cap's edges wherever the far cap's were
/// crossed. The placement is folded in as its hash, which keeps the key a
/// plain value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct EdgeKey {
    pub(crate) node: ogeom_topo::TShapeId,
    pub(crate) placement: u64,
}

impl EdgeKey {
    pub(crate) fn of(edge: &Shape) -> Self {
        use core::hash::{Hash as _, Hasher as _};
        // An edge where its node puts it, the usual case, needs no hash.
        if edge.location().is_identity() {
            return Self {
                node: edge.node(),
                placement: 0,
            };
        }
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        edge.location().hash(&mut hasher);
        Self {
            node: edge.node(),
            placement: hasher.finish(),
        }
    }

    /// The node's index, for the debug traces.
    pub(crate) fn index(self) -> u32 {
        self.node.index()
    }
}

/// What a boolean reads of one edge occurrence, the same for the faces on
/// either side of it.
pub(crate) struct Occurrence {
    key: EdgeKey,
    /// The edge's curve in world space, placement applied.
    curve: std::sync::Arc<Curve>,
    /// Each end where the curve puts it, with the tolerance of the vertex
    /// there, matched to the curve's ends by position: the edge's vertex
    /// order and its curve's direction need not agree once the edge is
    /// reversed in its wire.
    ends: [(Point, f64); 2],
    /// Each end's vertex where the vertex stands, which a file may put
    /// between the ends of the curves it closes, each within its radius;
    /// the curve's end where it has none.
    vertices: [Point; 2],
    /// The loosest tolerance either of its vertices states.
    ends_tolerance: f64,
    /// The box of seventeen points along the curve, the box of the chord
    /// midpoints between them, and twice the farthest a midpoint stands off
    /// its chord: how far the true curve may sag from the sampled polyline.
    sampled: (ogeom_math::Aabb, ogeom_math::Aabb, f64),
}

impl Occurrence {
    fn of(
        model: &Model,
        edge: &Shape,
        curve: Curve,
        range: (f64, f64),
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        let vertices = model.children_of(edge)?;
        let ends_tolerance = vertices
            .iter()
            .filter_map(|v| model.node(v).and_then(|n| n.data().as_vertex()))
            .fold(0.0_f64, |acc, d| acc.max(d.tolerance.get()));
        let (a, b) = (curve.point_at(range.0, tol)?, curve.point_at(range.1, tol)?);
        let mut ends = [(a, tol.confusion()), (b, tol.confusion())];
        let mut at_vertex = [a, b];
        for v in &vertices {
            if let Some(d) = model.node(v).and_then(|n| n.data().as_vertex()) {
                let at = v.transform(model.datums())?.apply(d.point);
                let k = usize::from(at.distance(b) < at.distance(a));
                ends[k].1 = ends[k].1.max(d.tolerance.get());
                at_vertex[k] = at;
            }
        }
        let mut previous: Option<Point> = None;
        let (mut points, mut mids) = (ogeom_math::Aabb::EMPTY, ogeom_math::Aabb::EMPTY);
        let mut sag = 0.0_f64;
        for i in 0..=16 {
            #[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
            let t = range.0 + (range.1 - range.0) * f64::from(i) / 16.0;
            let p = curve.point_at(t, tol)?;
            if let Some(q) = previous {
                let step = (range.1 - range.0) / 32.0;
                let mid = curve.point_at(t - step, tol)?;
                sag = sag.max(mid.distance(Point::midpoint(q, p)) * 2.0);
                mids = mids.with_point(mid);
            }
            points = points.with_point(p);
            previous = Some(p);
        }
        Ok(Self {
            key: EdgeKey::of(edge),
            curve: std::sync::Arc::new(curve),
            ends,
            vertices: at_vertex,
            ends_tolerance,
            sampled: (points, mids, sag),
        })
    }
}

/// One boundary strand source: an edge as one face uses it, with the pcurve
/// that side of it.
pub(crate) struct BoundaryEdge {
    /// The edge occurrence's identity, for sharing paves across faces.
    pub(crate) node: EdgeKey,
    /// The curve in world space, shared by the faces the edge bounds.
    pub(crate) curve: std::sync::Arc<Curve>,
    /// The portion the edge covers, in the curve's parameter.
    pub(crate) crange: (f64, f64),
    /// The pcurve on this face's surface.
    pub(crate) pcurve: PlanarCurve,
    /// The portion of the pcurve, mapped proportionally from `crange`.
    pub(crate) prange: (f64, f64),
    /// The other side of a seam, where this edge is one.
    pub(crate) other_side: Option<(PlanarCurve, (f64, f64))>,
    /// The edge's two seam columns as its node states them, forward then
    /// reversed, kept where this face uses only one: on a surface several
    /// faces stand on, the edge between two of them is a seam of that
    /// surface, each face walking its own column.
    pub(crate) seam: Option<(PlanarCurve, PlanarCurve)>,
    /// The radius within which the edge honestly lies: a fitted rail can
    /// carry a few dozen microns of construction slop, and every filter that
    /// compares this edge's curve against exact geometry widens by it.
    pub(crate) tolerance: f64,
    /// How far the edge stands off this face: its curve from the face's
    /// own image of it, the pcurve on the surface, at most its tolerance.
    /// An edge's tolerance is one radius for both its faces, and it may be
    /// owed to one of them alone: a chord between a facet and the drum
    /// it was converted beside lies in the facet's plane and stands off
    /// the drum by its sag. On the facet the chord is exact, and what
    /// meets it there is read as exactly as the facet's other edges. A
    /// tolerance no face of the edge stands off by half of is owed to
    /// none, and is this too.
    pub(crate) standoff: f64,
    /// The radius the edge holds to on this face as measured: its stated
    /// tolerance, or the curve's distance from the face's own image of it
    /// where that is wider.
    pub(crate) residual: f64,
    /// Whether the face walks the edge along its curve (`true`) or back,
    /// keeping its material on the left about its outward normal; `None`
    /// where it walks it both ways, a seam.
    pub(crate) walked: Option<bool>,
    /// The edge as this face draws it, where it stands off the face by more
    /// than [`DRAWN_APART`]: the pcurve on the surface, over `prange`. A
    /// chord between a facet and the drum it was converted beside is drawn
    /// on the drum as the arc beneath it, and what meets the drum's trim
    /// meets that arc: a section running a hair outside the arc, read
    /// against the chord its sag away, runs along the edge it never
    /// reaches.
    pub(crate) image: Option<std::sync::Arc<Curve>>,
    /// The looser of the edge's two end vertices' own radii: the doubt an
    /// earlier junction recorded there, which the arrangement's node snap
    /// on this face must reach.
    pub(crate) ends_tolerance: f64,
    /// Each end vertex in space with its own recorded radius: a vertex that
    /// owns a span is a junction of this boolean too, and every strand end
    /// inside it names it.
    pub(crate) ends: [(Point, f64); 2],
    /// Each end's vertex where the vertex itself stands (see
    /// [`Occurrence::vertices`]).
    pub(crate) vertices: [Point; 2],
    /// A conservative box round the edge: its samples, grown by twice the
    /// sag measured between them and by its own radius.
    pub(crate) bound: ogeom_math::Aabb,
}

/// How far an edge must stand off a face before the face reads it as its
/// own image of it: a hundred confusions, past the rounding of any edge
/// whose curve and pcurve were made together.
pub(crate) const DRAWN_APART: f64 = 1e2;

/// How many of a loose edge's tolerances a crossing of it may stand along
/// it from the crossing the face across it finds: where the edge leans out
/// of the crossing face by less than a tenth of a radian, the two are not
/// told apart by the lean but by the touch, and held to this.
pub(crate) const LEAN_SPREAD: f64 = 10.0;

impl BoundaryEdge {
    /// The curve this face's boundary runs along, its range, and how far
    /// the face may stand off it: the image where there is one, exact on
    /// the face, or the edge's curve and its standoff.
    pub(crate) fn drawn(&self) -> (&Curve, (f64, f64), f64) {
        match &self.image {
            Some(image) => (image, ordered(self.prange), 0.0),
            None => (&self.curve, self.crange, self.standoff),
        }
    }

    /// A parameter of [`BoundaryEdge::drawn`] as the edge curve's: the
    /// image runs over the pcurve's range, mapped onto the curve's in
    /// proportion, as the face's strands draw it.
    pub(crate) fn onto_curve(&self, s: f64) -> f64 {
        if self.image.is_some() {
            rescale(s, self.prange, self.crange)
        } else {
            s
        }
    }

    /// A parameter of the edge curve as [`BoundaryEdge::drawn`]'s.
    pub(crate) fn onto_drawn(&self, t: f64) -> f64 {
        if self.image.is_some() {
            rescale(t, self.crange, self.prange)
        } else {
            t
        }
    }
}

/// A range low end first.
pub(crate) fn ordered(range: (f64, f64)) -> (f64, f64) {
    if range.0 <= range.1 {
        range
    } else {
        (range.1, range.0)
    }
}

/// A pole: an edge that bounds a face in parameter space and collapses to
/// a point in space.
///
/// A sphere's poles and a cone's apex have no curve to split, but they are
/// half the chart's boundary. Leave them out and the face's outline does
/// not close, and nothing can be arranged inside it.
pub(crate) struct PoleEdge {
    /// Where the pole sits in space.
    pub(crate) point: Point,
    /// The chart line the pole runs along.
    pub(crate) pcurve: PlanarCurve,
    pub(crate) prange: (f64, f64),
}

/// One face, gathered and vetted.
pub(crate) struct GFace {
    pub(crate) face: Shape,
    /// The surface in world space, placement applied.
    pub(crate) surface: SurfaceGeometry,
    /// A conservative world bound: the boundary edges' sampled extent plus
    /// the surface's own allowance: nothing for a plane, measured sampling
    /// slack for the ruled kinds whose rulings pin them to their boundary's
    /// hull, most of the diagonal for anything that can genuinely bulge, a
    /// dome past its equator. The poles join after. Gates the pair filter
    /// and refusals. The `bool-audit` diagnostics channel audits its conservatism.
    pub(crate) bound: ogeom_math::Aabb,
    /// What the pair tests read: [`GFace::bound`] cut down to the box the
    /// model keeps for the face, which stands off a dome by its tolerance
    /// where `bound` keeps most of its diagonal as allowance. The scales
    /// taken from a face's size still read `bound`.
    pub(crate) filter: ogeom_math::Aabb,
    /// The scale the marching chord derives from: deliberately *not* the
    /// filter box's diagonal, so the filter can tighten without silently
    /// tightening the marcher.
    pub(crate) chord_scale: f64,
    /// How far the face may stand off its surface, as it states: a face
    /// fitted to a mesh lies within its fit's deviation of it.
    pub(crate) tolerance: f64,
    pub(crate) edges: Vec<BoundaryEdge>,
    pub(crate) poles: Vec<PoleEdge>,
    /// The edges its rings walk out along and straight back ([`slits`]),
    /// left out of `edges`.
    pub(crate) slits: ogeom_core::FastSet<EdgeKey>,
    /// The face's boundary polylined in its chart within [`strand_sag`],
    /// both seam columns, ends welded: the trim [`chart_point_within`]
    /// tests against. Built the first time a probe asks and kept, since the
    /// same face is asked about once per piece that might lie on it. `None`
    /// where a pcurve would not polyline.
    outline: std::sync::OnceLock<Option<Vec<Vec<Point2>>>>,
    /// [`GFace::outline`] with its lines' boxes, for asking whether many
    /// points lie inside it.
    outline_trim: std::sync::OnceLock<Option<arrange::Trim>>,
    /// The holes left out of the arrangement (see [`aside`]): each lies
    /// inside one piece, which takes it back.
    pub(crate) holes_aside: std::sync::Arc<[aside::Hole]>,
    /// The surface the face's node stands on, where the face is rebuilt on
    /// it: a face holding edges it shares with faces set aside, which keep
    /// their pcurves on it.
    pub(crate) own_surface: Option<ogeom_topo::SurfaceId>,
    /// The stored surface and the placement it stands at: faces with the
    /// same carrier stand on one surface in space.
    pub(crate) carrier: (ogeom_topo::SurfaceId, Location),
    /// The face's trim sampled coarsely for folding a chart image inside
    /// it: [`face_trim_lines`], kept for the face's every sub-edge.
    trim_lines: std::sync::OnceLock<Vec<Vec<Point2>>>,
}

impl GFace {
    pub(crate) fn trim_lines(&self, tol: Tolerances) -> &[Vec<Point2>] {
        self.trim_lines.get_or_init(|| face_trim_lines(self, tol))
    }

    pub(crate) fn outline_trim(&self, tol: Tolerances) -> Option<&arrange::Trim> {
        self.outline_trim
            .get_or_init(|| {
                self.outline(tol)
                    .map(|lines| arrange::Trim::new(lines.to_vec()))
            })
            .as_ref()
    }

    pub(crate) fn outline(&self, tol: Tolerances) -> Option<&[Vec<Point2>]> {
        self.outline
            .get_or_init(|| {
                let mut lines: Vec<Vec<Point2>> = Vec::new();
                for e in &self.edges {
                    lines.push(
                        pcurve_polyline_within_sag(
                            &e.pcurve,
                            e.prange,
                            e.crange,
                            e.crange,
                            &self.surface,
                            tol,
                        )
                        .ok()?,
                    );
                    // A seam bounds the chart twice (once per column), and
                    // a trim test that sees only one side reads half the
                    // band as outside.
                    if let Some((other, orange)) = &e.other_side {
                        lines.push(
                            pcurve_polyline_within_sag(
                                other,
                                *orange,
                                e.crange,
                                e.crange,
                                &self.surface,
                                tol,
                            )
                            .ok()?,
                        );
                    }
                }
                lines.extend(self.pole_rows(tol).ok()?);
                weld_outline_ends(&mut lines, self, tol);
                bridge_outline_gaps(&mut lines, param_snap(&self.surface, tol));
                Some(lines)
            })
            .as_deref()
    }

    /// Each pole's row in the chart, end to end: the stretch of the
    /// outline a collapsed side stands for. A row across the chart's
    /// columns (a sphere's pole) is level with a containment ray, but one
    /// along them (a patch whose every row ends in one point) is crossed
    /// like any other side.
    pub(crate) fn pole_rows(&self, tol: Tolerances) -> OgeomResult<Vec<Vec<Point2>>> {
        self.poles
            .iter()
            .map(|pole| {
                Ok(vec![
                    pole.pcurve.point_at(pole.prange.0, tol)?,
                    pole.pcurve.point_at(pole.prange.1, tol)?,
                ])
            })
            .collect()
    }
}

/// Whether `shape` is a solid, or a compound of solids and nothing else:
/// disjoint lumps (a fuse of parts that do not meet, a fillet's wedges
/// taken together) bound their material as one solid's shells do.
pub(crate) fn is_solid_or_lumps(model: &Model, shape: &Shape) -> OgeomResult<bool> {
    Ok(match model.kind_of(shape)? {
        ShapeType::Solid => true,
        ShapeType::Compound => {
            let parts = model.children_of(shape)?;
            !parts.is_empty()
                && parts
                    .iter()
                    .map(|part| model.kind_of(part))
                    .collect::<OgeomResult<Vec<_>>>()?
                    .iter()
                    .all(|kind| *kind == ShapeType::Solid)
        }
        _ => false,
    })
}

/// Whether `shape` is a sheet: a face, a shell, or a compound of faces and
/// shells, with no solid behind any of it.
pub(crate) fn is_sheet(model: &Model, shape: &Shape) -> OgeomResult<bool> {
    Ok(match model.kind_of(shape)? {
        ShapeType::Face | ShapeType::Shell => true,
        ShapeType::Compound => {
            let parts = model.children_of(shape)?;
            !parts.is_empty()
                && parts
                    .iter()
                    .map(|part| model.kind_of(part))
                    .collect::<OgeomResult<Vec<_>>>()?
                    .iter()
                    .all(|kind| matches!(kind, ShapeType::Face | ShapeType::Shell))
        }
        _ => false,
    })
}

/// Every edge of a sheet's faces given its exact pcurve on the face where it
/// has none: a face built on a bare boundary is split in its chart like any
/// other.
pub(crate) fn trim_sheet_edges(
    model: &mut Model,
    sheet: &Shape,
    tol: Tolerances,
) -> OgeomResult<()> {
    for face in explore_unique(model, sheet, ShapeType::Face)? {
        ogeom_algo::build::trimmed_where_bare(model, &face, tol)?;
    }
    Ok(())
}

/// An argument: a solid, or a sheet where the operation takes one.
pub(crate) struct GSolid {
    pub(crate) solid: Shape,
    pub(crate) faces: Vec<GFace>,
    /// The faces set aside, passed through untouched, each with its place
    /// among the solid's faces.
    pub(crate) aside: Vec<(Shape, usize)>,
    /// The edges the gathered faces share with the faces set aside, each
    /// as it stands: rebuilt as itself, so the faces set aside still meet
    /// the rebuilt ones along it.
    pub(crate) shared: ogeom_core::FastMap<EdgeKey, Shape>,
    /// The vertices of the shared edges, each where it stands: the strands
    /// ending there end on it.
    pub(crate) shared_vertices: Vec<(Point, Shape)>,
    /// Each face's place in the solid, where faces are set aside.
    pub(crate) order: ogeom_core::FastMap<ogeom_topo::SameKey, usize>,
    /// What joins the faces set aside to the rest.
    pub(crate) beside: std::sync::Arc<aside::Beside>,
    /// Whether the solid is one shell already found closed.
    pub(crate) closed: bool,
}

/// What the general fuse makes of its two arguments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Operands {
    /// Two solids: each one's faces split by the other's and classified
    /// against it.
    Solids,
    /// A sheet and a solid: the sheet's faces split by the solid's and
    /// classified against it. The solid's faces only cut.
    SheetBySolid,
    /// Two sheets: the first one's faces split by the second's and left
    /// unclassified. The second one's faces only cut.
    SheetBySheet,
    /// A solid and itself: each face split by every face of the solid but
    /// itself and the faces it shares an edge or a corner with, and left
    /// unclassified. Only the first argument's pieces are made.
    Itself,
}

pub(crate) fn gather(
    model: &Model,
    solid: &Shape,
    sheet: bool,
    aside: Option<&aside::Aside>,
    tol: Tolerances,
) -> OgeomResult<GSolid> {
    let mut read_closed = false;
    if sheet {
        if !is_sheet(model, solid)? {
            ogeom_bail!(
                Construction,
                "a sheet argument is a face, a shell or a compound of those"
            );
        }
    } else {
        if !is_solid_or_lumps(model, solid)? {
            ogeom_bail!(Construction, "boolean arguments are solids");
        }
        let shells = explore_unique(model, solid, ShapeType::Shell)?;
        // A solid of one shell whose every edge its faces walk an even
        // number of times, as setting faces aside read it, is closed.
        read_closed =
            shells.len() == 1 && aside.is_some_and(|aside| !aside.is_empty() && aside.closed);
        for shell in shells {
            if !read_closed && !is_shell_closed(model, &shell)? {
                ogeom_bail!(Construction, "an open shell bounds no volume to operate on");
            }
        }
    }

    let mut faces = Vec::new();
    let mut shared = ogeom_core::FastMap::default();
    // Each edge occurrence read once for the faces on either side of it:
    // its world curve, its ends, and its sampled extent.
    let mut occurrences: ogeom_core::FastMap<(ogeom_topo::TShapeId, Location), Occurrence> =
        ogeom_core::FastMap::default();
    let aside = aside.filter(|aside| !aside.is_empty());
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        if aside.is_some_and(|aside| {
            !aside
                .gathered
                .contains_key(&ogeom_topo::SameKey(face.clone()))
        }) {
            continue;
        }
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(stored) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let placement = face.transform(model.datums())?;
        if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 {
            ogeom_bail!(
                NotDone,
                "a scaled placement changes a surface's parameterization out \
                 from under its pcurves; bake the scale before a boolean"
            );
        }
        let surface = stored.transformed(&placement, tol)?;
        let chart_snap = param_snap(&surface, tol);
        let surface_id = data.surface;
        let tolerance = data.tolerance.get();

        // A plane's holes clear of the other solid are left out: nothing
        // crosses them, and each goes back into the piece holding it.
        let holes_aside = aside
            .and_then(|aside| aside.holes.get(&ogeom_topo::SameKey(face.clone())))
            .cloned()
            .unwrap_or_else(|| std::sync::Arc::from(Vec::new()));
        let slits = slits(model, &face, &surface, surface_id, tol)?;
        // How the face walks each of its edges, senses composed: forward
        // along the curve, back, or both ways (a seam).
        let mut walks: ogeom_core::FastMap<EdgeKey, (bool, bool)> = ogeom_core::FastMap::default();
        for edge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            let forward = edge.orientation() == ogeom_topo::Orientation::Forward;
            let walk = walks.entry(EdgeKey::of(&edge)).or_default();
            if forward {
                walk.0 = true;
            } else {
                walk.1 = true;
            }
        }
        let mut face_edges = if holes_aside.is_empty() {
            explore_unique(model, &face, ShapeType::Edge)?
        } else {
            let wires = model.children_of(&face)?;
            let mut left = vec![false; wires.len()];
            for hole in holes_aside.iter() {
                left[hole.at] = true;
            }
            let mut seen = ogeom_core::FastSet::default();
            let mut edges = Vec::new();
            for (wire, _) in wires.iter().zip(&left).filter(|(_, left)| !**left) {
                for edge in explore_unique(model, wire, ShapeType::Edge)? {
                    if seen.insert(ogeom_topo::SameKey(edge.clone())) {
                        edges.push(edge);
                    }
                }
            }
            edges
        };
        face_edges.retain(|edge| !slits.contains(&EdgeKey::of(edge)));
        let mut shares = false;

        let mut edges = Vec::new();
        let mut samples: Vec<(ogeom_math::Aabb, ogeom_math::Aabb, f64)> = Vec::new();
        let mut poles = Vec::new();
        let mut pole_keys = Vec::new();
        for edge in face_edges {
            if let Some(aside) = aside
                && aside.edges.contains(&EdgeKey::of(&edge))
            {
                shared.insert(
                    EdgeKey::of(&edge),
                    edge.oriented(ogeom_topo::Orientation::Forward),
                );
                shares = true;
            }
            let Some(edge_node) = model.node(&edge) else {
                ogeom_bail!(Dangling, "edge is not in this model");
            };
            let NodeData::Edge(edge_data) = edge_node.data() else {
                ogeom_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = edge_data.curve3d() else {
                // A degenerate edge (a sphere's pole, a cone's apex) has
                // no extent to split, but it *does* bound the chart, and a
                // chart whose top is missing bounds nothing.
                if let Some(EdgeRepr::PCurve {
                    curve: pc, range, ..
                }) = edge_data.pcurve_for(surface_id, edge.location())
                {
                    let Some(planar) = model.geometry().pcurve(*pc) else {
                        ogeom_bail!(Dangling, "pcurve is not in this model");
                    };
                    let at = explore_unique(model, &edge, ShapeType::Vertex)?;
                    let Some(point) = at
                        .first()
                        .and_then(|v| model.node(v))
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                    else {
                        ogeom_bail!(Construction, "a pole with no vertex is nowhere");
                    };
                    let placed = edge.transform(model.datums())?.apply(point);
                    let (pcurve, prange) = straight_pole(planar, *range, tol)?;
                    poles.push(PoleEdge {
                        point: placed,
                        pcurve,
                        prange,
                    });
                    pole_keys.push(EdgeKey::of(&edge));
                }
                continue;
            };
            let occurrence = match occurrences.entry((edge.node(), edge.location().clone())) {
                ogeom_core::collections::hash_map::Entry::Occupied(known) => known.into_mut(),
                ogeom_core::collections::hash_map::Entry::Vacant(slot) => {
                    let Some(geometry) = model.geometry().curve(*curve) else {
                        ogeom_bail!(Dangling, "curve is not in this model");
                    };
                    let placed = geometry.transformed(&edge.transform(model.datums())?, tol)?;
                    slot.insert(Occurrence::of(model, &edge, placed, *range, tol)?)
                }
            };
            let world = occurrence.curve.clone();
            samples.push(occurrence.sampled);
            let (pcurve, prange, other_side, seam) =
                match edge_data.pcurve_for(surface_id, edge.location()) {
                    Some(EdgeRepr::PCurve {
                        curve: pc, range, ..
                    }) => {
                        let Some(planar) = model.geometry().pcurve(*pc) else {
                            ogeom_bail!(Dangling, "pcurve is not in this model");
                        };
                        (planar.clone(), *range, None, None)
                    }
                    Some(EdgeRepr::Seam {
                        forward,
                        reversed,
                        range,
                        ..
                    }) => {
                        let (Some(f), Some(r)) = (
                            model.geometry().pcurve(*forward),
                            model.geometry().pcurve(*reversed),
                        ) else {
                            ogeom_bail!(Dangling, "seam pcurve is not in this model");
                        };
                        (
                            f.clone(),
                            *range,
                            Some((r.clone(), *range)),
                            Some((f.clone(), r.clone())),
                        )
                    }
                    _ => ogeom_bail!(
                        Construction,
                        "an edge with no pcurve on its face cannot be split in \
                         that face's parameter space"
                    ),
                };
            let edge_tolerance = edge_data.tolerance.get();
            let mut standoff = edge_tolerance;
            let mut residual = edge_tolerance;
            if edge_tolerance > tol.confusion() * 10.0 {
                let gap = |pc: &PlanarCurve, pr: (f64, f64)| {
                    ogeom_algo::pcurve_fit::lifted_gap(
                        (&world, *range),
                        (pc, pr),
                        &surface,
                        false,
                        tol,
                    )
                    .unwrap_or(edge_tolerance)
                };
                let mut measured = gap(&pcurve, prange);
                if let Some((other, orange)) = &other_side {
                    measured = measured.max(gap(other, *orange));
                }
                standoff = measured.min(edge_tolerance);
                residual = measured.max(edge_tolerance);
            }

            edges.push(BoundaryEdge {
                node: occurrence.key,
                curve: world,
                crange: *range,
                pcurve,
                prange,
                other_side,
                seam,
                tolerance: edge_tolerance,
                standoff,
                residual,
                walked: match walks.get(&occurrence.key) {
                    Some((true, false)) => Some(true),
                    Some((false, true)) => Some(false),
                    _ => None,
                },
                image: None,
                ends_tolerance: occurrence.ends_tolerance,
                ends: occurrence.ends,
                vertices: occurrence.vertices,
                bound: ogeom_math::Aabb::EMPTY,
            });
        }
        if edges.is_empty() {
            ogeom_bail!(Construction, "a face with no boundary bounds nothing");
        }

        if aside.is_none() || !(shares || !holes_aside.is_empty()) {
            fold_seamed_rings(
                model,
                &face,
                &surface,
                &mut edges,
                (&poles, &pole_keys),
                &slits,
                tol,
            )?;
        }
        // A seam edge bounds the chart twice *only when this face wraps*:
        // when its other boundary reaches both columns, as a full drum's
        // rims do. A face that merely sits against the meridian (a sphere
        // octant whose boundary happens to be the seam) uses one column,
        // and feeding the far copy into the arrangement leaves a strand
        // nothing connects to. The decision is made here, once, by chart
        // connectivity, and everything downstream (the arrangement, the
        // trim tests, the rebuild) inherits it.
        let seam_indices: Vec<usize> = edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.other_side.is_some())
            .map(|(i, _)| i)
            .collect();
        for i in seam_indices {
            // Each other end in the chart, with how far its edge may stand
            // off where it meets the next: its own radius, or its ends'.
            let mut pool: Vec<(Point2, f64)> = Vec::new();
            for (j, e) in edges.iter().enumerate() {
                if j == i {
                    continue;
                }
                let reach = e.tolerance.max(e.ends_tolerance);
                for t in [e.prange.0, e.prange.1] {
                    pool.push((e.pcurve.point_at(t, tol)?, reach));
                }
                if let Some((other, orange)) = &e.other_side {
                    for t in [orange.0, orange.1] {
                        pool.push((other.point_at(t, tol)?, reach));
                    }
                }
            }
            for p in &poles {
                for t in [p.prange.0, p.prange.1] {
                    pool.push((p.pcurve.point_at(t, tol)?, 0.0));
                }
            }
            let own = edges[i].tolerance.max(edges[i].ends_tolerance);
            // Two ends meet within the chart snap, or as far apart
            // across it as both their radii reach in space, twice over: a
            // loop's pcurve, sewn shut in space within its ends' radius, can
            // stop a hair short of the column it closes on. Read across the
            // chart, not in space, where the seam's two columns are one.
            let meets = |q: Point2, at: Point2, reach: f64| -> bool {
                let apart = q.distance(at);
                if apart <= chart_snap {
                    return true;
                }
                // A hair, not a turn: at a pole a column's whole width
                // costs nothing in space.
                if apart > SEAM_HAIR {
                    return false;
                }
                let middle = Point2::new(f64::midpoint(q.x, at.x), f64::midpoint(q.y, at.y));
                let Ok((du, dv)) = surface.d1_at(middle.x, middle.y, tol) else {
                    return false;
                };
                let across = (du.magnitude() * (q.x - at.x)).hypot(dv.magnitude() * (q.y - at.y));
                across <= 2.0 * (reach + own) + tol.confusion()
            };
            let connected = |pc: &PlanarCurve, range: (f64, f64)| -> OgeomResult<bool> {
                for t in [range.0, range.1] {
                    let at = pc.point_at(t, tol)?;
                    if !pool.iter().any(|&(q, reach)| meets(q, at, reach)) {
                        return Ok(false);
                    }
                }
                Ok(true)
            };
            let e = &edges[i];
            let primary_connects = connected(&e.pcurve, e.prange)?;
            let (other_pc, orange) = e
                .other_side
                .clone()
                .unwrap_or_else(|| (e.pcurve.clone(), e.prange));
            let other_connects = connected(&other_pc, orange)?;
            match (primary_connects, other_connects) {
                // Both columns meet the rest of the boundary: the face
                // wraps, and both belong.
                (true, true) => {}
                // One column is this face's. The far copy would dangle.
                (true, false) => edges[i].other_side = None,
                (false, true) => {
                    let e = &mut edges[i];
                    e.pcurve = other_pc;
                    e.prange = orange;
                    e.other_side = None;
                }
                // Neither connects: leave both, and the arrangement's own
                // refusal names the failure.
                (false, false) => {}
            }
        }
        // Every ring is folded onto one branch of its chart, a face rebuilt
        // on its own surface beside faces set aside too. Such a face keeps
        // the pcurves of the edges it shares with them as stated. Along a
        // direction the surface repeats in, a stated pcurve whole periods
        // from where the fold reads it names the same curve, as it did in
        // the face given. Along one the chart only closes in, the two lie on
        // different lines of the closure, and the fuse gathers every face
        // instead.
        let stated: Vec<Option<Point2>> = edges
            .iter()
            .map(|e| {
                aside
                    .is_some_and(|aside| aside.edges.contains(&e.node))
                    .then(|| {
                        e.pcurve
                            .point_at(f64::midpoint(e.prange.0, e.prange.1), tol)
                    })
                    .transpose()
            })
            .collect::<OgeomResult<_>>()?;
        fold_rings(
            model,
            &face,
            &surface,
            &mut edges,
            (&poles, &pole_keys),
            &slits,
            tol,
        )?;
        let repeats = [surface.is_periodic_u(), surface.is_periodic_v()];
        for (e, at) in edges.iter().zip(&stated) {
            let Some(at) = at else { continue };
            let now = e
                .pcurve
                .point_at(f64::midpoint(e.prange.0, e.prange.1), tol)?;
            if ((now.x - at.x).abs() > chart_snap && !repeats[0])
                || ((now.y - at.y).abs() > chart_snap && !repeats[1])
            {
                ogeom_bail!(
                    NotDone,
                    "a face sharing edges with faces set aside has its ring stated \
                     across its chart's join on one of them"
                );
            }
        }
        let moved = bring_trim_home(&surface, &mut edges, &mut poles, tol)?;
        close_pole_joints(
            model,
            &face,
            &surface,
            (&edges, &slits),
            &mut poles,
            &pole_keys,
            tol,
        )?;
        // Rebuilt around edges it shares with faces set aside, the face
        // stays on its own surface, where those edges' pcurves are, in the
        // chart they are stated in.
        let own_surface = if aside.is_some() && (shares || !holes_aside.is_empty()) {
            if moved {
                ogeom_bail!(
                    NotDone,
                    "a face sharing edges with faces set aside has its trim stated \
                     a period off its chart"
                );
            }
            Some(surface_id)
        } else {
            None
        };
        let mut bound = ogeom_math::Aabb::EMPTY;
        // For a plane or a ruled surface the box is trusted to the
        // boundary's own hull, so the boundary's sampling slack must be
        // measured: how far the true edge sags from the 16-chord polyline,
        // read at each chord's midpoint and doubled for the sag's asymmetry.
        // A disc's rim bows out between its samples by more than a
        // cylinder's cap pokes past a wall it all but meets.
        let ruled = matches!(
            &surface,
            SurfaceGeometry::Plane(_) | SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_)
        );
        let mut slack = 0.0_f64;
        for (e, &(points, mids, sag)) in edges.iter_mut().zip(&samples) {
            if ruled {
                slack = slack.max(sag);
                bound = bound.union(&mids);
            }
            bound = bound.union(&points);
            e.bound = points
                .union(&mids)
                .expanded(sag + e.tolerance + tol.confusion() * 1e2);
        }
        // A plane never bulges past its boundary. A ruled surface (cylinder,
        // cone) cannot either: every surface point lies on a straight ruling
        // whose ends are on the boundary, so the face sits inside its
        // boundary's hull and only the boundary's own sampling slack is needed,
        // the plane's as much as the ruled surface's.
        // Anything else may genuinely bulge (a dome past its equator) and
        // keeps most of its own diagonal as allowance. The `bool-audit`
        // channel's audit holds every arm to conservatism.
        let bulge = match &surface {
            SurfaceGeometry::Plane(_) | SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_) => {
                slack
            }
            _ => bound.diagonal() * 0.75,
        };
        // The scale the marching chord is derived from, decoupled from the
        // filter box: the diagonal with a blanket three-quarter bulge for
        // every curved surface. The chord is a tolerance, not a bound, so
        // tightening the filter must not silently tighten the marcher,
        // which is what a single box feeding both does.
        let margin = tol.confusion() * 1e2;
        let chord_scale = match &surface {
            SurfaceGeometry::Plane(_) => bound.expanded(margin).diagonal(),
            _ => bound.expanded(bound.diagonal() * 0.75 + margin).diagonal(),
        };
        // A pole bounds the chart with no edge to sample: a cone drilled to
        // its apex reaches the apex, and a bound that omits it would let the
        // filter drop a pair the apex genuinely meets. It joins *after* the
        // bulge is taken from the boundary's own diagonal: the pole is an
        // exact point, and letting it stretch the diagonal would inflate a
        // dome's allowance by its own height over again.
        let mut bound = bound.expanded(bulge);
        for pole in &poles {
            bound = bound.with_point(pole.point);
        }
        let bound = bound.expanded(tol.confusion() * 1e2);
        let filter = kept_filter(model, &face, &surface, bound, tol);
        let carrier = (surface_id, face.location().clone());
        faces.push(GFace {
            poles,
            slits,
            face,
            surface,
            bound,
            filter,
            chord_scale,
            tolerance,
            edges,
            outline: std::sync::OnceLock::new(),
            outline_trim: std::sync::OnceLock::new(),
            trim_lines: std::sync::OnceLock::new(),
            holes_aside,
            own_surface,
            carrier,
        });
    }
    if faces.is_empty() {
        if sheet {
            ogeom_bail!(Construction, "a sheet with no faces has nothing to split");
        }
        ogeom_bail!(Construction, "a solid with no faces bounds nothing");
    }
    let mut shared_vertices = Vec::new();
    let mut seen = ogeom_core::FastSet::default();
    for edge in shared.values() {
        for vertex in model.children_of(edge)? {
            if !seen.insert(ogeom_topo::SameKey(vertex.clone())) {
                continue;
            }
            let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
                ogeom_bail!(Dangling, "vertex is not in this model");
            };
            let at = vertex.transform(model.datums())?.apply(data.point);
            shared_vertices.push((at, vertex.oriented(ogeom_topo::Orientation::Forward)));
        }
    }
    // An edge's tolerance is owed to the face whose image of it stands that
    // far off, and each other face reads the edge as it stands off that
    // face. One no face of it stands off by half as much (an edge stated
    // loose, a rail welded across a gap an earlier operation left) owes it
    // to none, and every face reads the edge to its tolerance.
    let mut owed: ogeom_core::FastMap<EdgeKey, f64> = ogeom_core::FastMap::default();
    for e in faces.iter().flat_map(|f: &GFace| f.edges.iter()) {
        let most = owed.entry(e.node).or_insert(0.0);
        *most = most.max(e.standoff);
    }
    for face in &mut faces {
        for e in &mut face.edges {
            if owed
                .get(&e.node)
                .is_some_and(|most| most * 2.0 < e.tolerance)
            {
                e.standoff = e.tolerance;
                continue;
            }
            // Standing off the face, the edge is drawn as the face draws it,
            // its pcurve read against its curve in proportion as the face's
            // strands draw it. Where the two are paced apart by more than
            // twice the edge's tolerance, the reach a loose edge's junction
            // takes, that reading would put the face's stops a pace off the
            // edge's splits, and the edge is read as it stands.
            // A periodic curve or pcurve (a circle, in space or in a chart)
            // is read through one window only, and parameters met along it
            // come back folded into its period: such an edge is left to its
            // curve.
            if e.standoff <= tol.confusion() * DRAWN_APART
                || e.other_side.is_some()
                || e.pcurve.is_periodic()
                || e.curve.is_periodic()
            {
                continue;
            }
            let paced = ogeom_algo::pcurve_fit::lifted_gap(
                (&e.curve, e.crange),
                (&e.pcurve, e.prange),
                &face.surface,
                true,
                tol,
            );
            if !paced.is_ok_and(|gap| gap <= e.tolerance * 2.0) {
                continue;
            }
            let (lo, hi) = ordered(e.prange);
            e.image = ogeom_geom::TrimmedCurve::new(
                Curve::OnSurface(Box::new(ogeom_geom::CurveOnSurface::new(
                    e.pcurve.clone(),
                    face.surface.clone(),
                ))),
                lo,
                hi,
                tol,
            )
            .ok()
            .map(|trimmed| std::sync::Arc::new(Curve::Trimmed(Box::new(trimmed))));
        }
    }
    Ok(GSolid {
        solid: solid.clone(),
        faces,
        aside: aside.map(|aside| aside.faces.clone()).unwrap_or_default(),
        shared,
        shared_vertices,
        order: aside
            .map(|aside| aside.gathered.clone())
            .unwrap_or_default(),
        beside: aside.map(|aside| aside.beside.clone()).unwrap_or_default(),
        closed: read_closed,
    })
}

/// `bound` cut down to the box the model keeps for `face`, widened by the
/// margin the box tests allow. The kept box is found exactly on the
/// elementary surfaces; on any other its interior extremes start from a
/// mesh, which may fall short of the surface by the mesh's chord, and that
/// chord is allowed for. A face whose box cannot be found keeps `bound`.
pub(crate) fn kept_filter(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    bound: ogeom_math::Aabb,
    tol: Tolerances,
) -> ogeom_math::Aabb {
    let Ok(kept) = ogeom_algo::face_bounds(model, face) else {
        return bound;
    };
    let exact = matches!(
        surface,
        SurfaceGeometry::Plane(_)
            | SurfaceGeometry::Cylinder(_)
            | SurfaceGeometry::Cone(_)
            | SurfaceGeometry::Sphere(_)
            | SurfaceGeometry::Torus(_)
    );
    let mut margin = tol.confusion() * 1e2;
    if !exact {
        margin += ogeom_mesh::Deflection::default().chord;
    }
    let filter = bound.intersection(&kept.expanded(margin));
    if filter.is_empty() { bound } else { filter }
}

/// The slits of a face: edges a ring walks out along and straight back, the
/// face on both sides, ending at a vertex no other edge of the ring reaches.
///
/// Such an edge bounds no material. The arrangement reads each strand as a
/// stretch of a cycle between nodes, and a slit's tip is a node of one
/// strand only. Left out of the face's boundary, the ring closes round the
/// region the face covers, and the pieces of it are rebuilt without it.
///
/// A ring is peeled at a joint where it takes the same edge, of two
/// distinct vertices, once each way, and where the joint is an ordinary
/// point of the chart. Peeling one slit may leave the edge it hangs from
/// walked out and back, and that edge is peeled in turn.
///
/// A seam to a pole is walked the same way, out along one column and back
/// along the other: at the pole the whole chart row is one point, and the
/// face wraps round it between the two columns. A joint is such a pole
/// where the surface's rate along a row through it, read at the row's ends
/// and middle, times the chart's whole width that way, stays within the
/// joint vertex's radius. A ring is never peeled there.
pub(crate) fn slits(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    surface_id: ogeom_topo::SurfaceId,
    tol: Tolerances,
) -> OgeomResult<ogeom_core::FastSet<EdgeKey>> {
    let mut slits = ogeom_core::FastSet::default();
    let ((u0, u1), (v0, v1)) = surface.domain();
    // Whether the joint where `edge` ends, as its ring walks it, is an
    // ordinary point of the chart, not on a row the surface collapses.
    let ordinary = |edge: &Shape| -> OgeomResult<bool> {
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            return Ok(false);
        };
        if !matches!(data.curve3d(), Some(EdgeRepr::Curve3d { .. })) {
            return Ok(false);
        }
        let Some((start, end)) = ogeom_algo::edge_vertices(model, edge)? else {
            return Ok(false);
        };
        if start.is_same(&end) {
            return Ok(false);
        }
        let Some(tip) = model.node(&end).and_then(|n| n.data().as_vertex()) else {
            return Ok(false);
        };
        let reach = tip.tolerance.get().max(tol.confusion());
        let (pcurve, range) = match data.pcurve_for(surface_id, edge.location()) {
            Some(EdgeRepr::PCurve { curve, range, .. }) => (curve, range),
            Some(EdgeRepr::Seam { forward, range, .. }) => (forward, range),
            _ => return Ok(false),
        };
        let Some(pcurve) = model.geometry().pcurve(*pcurve) else {
            return Ok(false);
        };
        let back = edge.orientation() == ogeom_topo::Orientation::Reversed;
        let at = pcurve.point_at(if back { range.0 } else { range.1 }, tol)?;
        for along_u in [true, false] {
            let (lo, hi) = if along_u { (u0, u1) } else { (v0, v1) };
            let width = hi - lo;
            let mut collapses = width.is_finite();
            for s in [0.0, 0.5, 1.0] {
                if !collapses {
                    break;
                }
                let q = if along_u {
                    Point2::new(lo + width * s, at.y)
                } else {
                    Point2::new(at.x, lo + width * s)
                };
                let Ok((du, dv)) = surface.d1_at(q.x, q.y, tol) else {
                    return Ok(false);
                };
                let rate = if along_u {
                    du.magnitude()
                } else {
                    dv.magnitude()
                };
                collapses = rate * width <= reach;
            }
            if collapses {
                return Ok(false);
            }
        }
        Ok(true)
    };
    let stored = face.oriented(ogeom_topo::Orientation::Forward);
    for wire in model.children_of(&stored)? {
        let mut ring = model.ordered_children_of(&wire)?;
        'peel: loop {
            let n = ring.len();
            for k in 0..n {
                let (out, back) = (&ring[k], &ring[(k + 1) % n]);
                if n < 2
                    || EdgeKey::of(out) != EdgeKey::of(back)
                    || out.orientation() == back.orientation()
                    || !ordinary(out)?
                {
                    continue;
                }
                slits.insert(EdgeKey::of(out));
                let (first, second) = (k.min((k + 1) % n), k.max((k + 1) % n));
                ring.remove(second);
                ring.remove(first);
                continue 'peel;
            }
            break;
        }
    }
    Ok(slits)
}

/// The rings of `face`, as it orients them, each without its slits.
pub(crate) fn rings_without(
    model: &Model,
    face: &Shape,
    slits: &ogeom_core::FastSet<EdgeKey>,
) -> OgeomResult<Vec<Vec<Shape>>> {
    model
        .children_of(face)?
        .iter()
        .map(|wire| {
            let mut ring = model.ordered_children_of(wire)?;
            ring.retain(|edge| !slits.contains(&EdgeKey::of(edge)));
            Ok(ring)
        })
        .collect()
}

/// Each of a face's rings walked in its chart, every edge's pcurve moved by
/// whole periods onto the branch that continues the ring, along each
/// direction the chart closes in, and each hole brought onto its outer
/// ring's branch.
///
/// Two faces sharing an edge share its pcurve, and on a periodic surface
/// it is stated on one branch: a quarter-round wall's straight edge on the
/// chart's closure may stand at `u = 0` while the arcs it meets end at
/// `u = 2 pi`. The ring is one closed curve on the surface, and in the
/// chart it jumps a period at one joint and back at another; the
/// arrangement reads two strands that end where nothing begins. The
/// mesher folds the same jumps back, edge by edge.
///
/// The walk starts on a ring's degenerate edge where it has one: a pole's
/// or an apex's chart line is the strand the arrangement draws there, and
/// is never moved, so the edges either side of it are carried onto its
/// branch, met in the chart, since every point of the line is the same
/// point of the surface. Failing that it starts on a seam the ring walks
/// once, a border with another face across the seam: read on the side the
/// ring walks it, it is never moved either, and the ring folded onto it
/// keeps that side alone. Elsewhere an edge is moved where its start stands
/// whole periods from the end walked to, unless the joint stands across a
/// degenerate row, and moved stands within both edges' radii of it on the
/// surface, measured through the middle of the chart segment between them.
/// A joint stands across a degenerate row where the surface's derivative
/// along the chart segment between its ends vanishes, not where the ends
/// merely lie close: a cylinder smaller than its edges' radii has ends a
/// period apart that close. A ring with no degenerate edge starts just
/// after a joint across an apex (two rulings meeting at one point stand
/// apart in the chart, and meet from any branch), and its first edge never
/// moves.
///
/// A ring is folded only where, folded, it closes where it began: one
/// winding round the surface keeps its jump. A ring walking a seam both
/// ways, or another edge twice, is left as stated, as is one walking a
/// seam once on a side it cannot tell. A hole whose chart extent lies
/// wholly beside its outer ring's is moved by the whole periods that put
/// it inside. On a chart that closes without repeating, a pcurve moves
/// only onto a stretch the chart covers.
pub(crate) fn fold_rings(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    edges: &mut [BoundaryEdge],
    (poles, pole_keys): (&[PoleEdge], &[EdgeKey]),
    slits: &ogeom_core::FastSet<EdgeKey>,
    tol: Tolerances,
) -> OgeomResult<()> {
    let chart_snap = param_snap(surface, tol);
    /// One edge as a ring walks it: a boundary edge, which may move; a
    /// degenerate one, which anchors the ring; or a seam the ring walks
    /// once, on the side it walks (an index into the ring's seam sides),
    /// which anchors it too; and whether the walk runs against its pcurve.
    #[derive(Clone, Copy)]
    enum Step {
        Edge(usize, bool),
        Pole(usize, bool),
        Seam(usize, bool),
    }
    let ((ua, ub), (va, vb)) = surface.domain();
    let periods = [
        ((surface.is_periodic_u() || surface.is_closed_u(tol)) && ub > ua).then_some(ub - ua),
        ((surface.is_periodic_v() || surface.is_closed_v(tol)) && vb > va).then_some(vb - va),
    ];
    if periods.iter().all(Option::is_none) {
        return Ok(());
    }
    let periodic = [surface.is_periodic_u(), surface.is_periodic_v()];
    // On a chart that closes without repeating, a pcurve is moved only onto
    // a stretch the chart covers: an edge on the closure, from one of its
    // two lines to the other.
    let inside = |pc: &PlanarCurve, range: (f64, f64)| -> OgeomResult<bool> {
        for k in 0..=2 {
            let p = pc.point_at(range.0 + (range.1 - range.0) * f64::from(k) / 2.0, tol)?;
            if (!periodic[0] && (p.x < ua - chart_snap || p.x > ub + chart_snap))
                || (!periodic[1] && (p.y < va - chart_snap || p.y > vb + chart_snap))
            {
                return Ok(false);
            }
        }
        Ok(true)
    };
    // The whole periods that bring `from` nearest `to`, along each closing
    // direction.
    let whole = |to: Point2, from: Point2| -> ogeom_math::Vector2 {
        let mut shift = ogeom_math::Vector2::new(0.0, 0.0);
        if let Some(period) = periods[0] {
            shift.x = ((to.x - from.x) / period).round() * period;
        }
        if let Some(period) = periods[1] {
            shift.y = ((to.y - from.y) / period).round() * period;
        }
        shift
    };
    let mut index: ogeom_core::FastMap<EdgeKey, (usize, usize)> = ogeom_core::FastMap::default();
    for (i, e) in edges.iter().enumerate() {
        index.insert(e.node, (i, 0));
    }
    let wires = rings_without(model, face, slits)?;
    for edge in wires.iter().flatten() {
        if let Some(slot) = index.get_mut(&EdgeKey::of(edge)) {
            slot.1 += 1;
        }
    }
    // How far apart two chart points stand on the surface, read through
    // the middle of the segment between them.
    let apart = |a: Point2, b: Point2| -> Option<f64> {
        let mid = Point2::new(f64::midpoint(a.x, b.x), f64::midpoint(a.y, b.y));
        let (p, m, q) = (
            surface.point_at(a.x, a.y, tol).ok()?,
            surface.point_at(mid.x, mid.y, tol).ok()?,
            surface.point_at(b.x, b.y, tol).ok()?,
        );
        Some(p.distance(m).max(m.distance(q)))
    };
    // Each ring's boundary edges where it was walked: the rings a hole may
    // be brought home among.
    let mut walked: Vec<Option<Vec<usize>>> = Vec::with_capacity(wires.len());
    for ring in &wires {
        let mut walk = Vec::with_capacity(ring.len());
        let mut kept = false;
        // Each seam the ring walks once, a border with another face rather
        // than the face's own closure: the edge and the side walked.
        let mut sides: Vec<(usize, PlanarCurve, (f64, f64))> = Vec::new();
        for edge in ring {
            let back = edge.orientation() == ogeom_topo::Orientation::Reversed;
            let key = EdgeKey::of(edge);
            if let Some(&(i, uses)) = index.get(&key) {
                let e = &edges[i];
                if uses == 1 && e.seam.is_some() {
                    match walked_seam_side(e, back, tol)? {
                        Some((pcurve, range)) => {
                            walk.push(Step::Seam(sides.len(), back));
                            sides.push((i, pcurve, range));
                        }
                        None => kept = true,
                    }
                    continue;
                }
                kept |= uses > 1 || e.seam.is_some() || e.other_side.is_some();
                walk.push(Step::Edge(i, back));
            } else if let Some(j) = pole_keys.iter().position(|k| *k == key) {
                walk.push(Step::Pole(j, back));
            }
        }
        let moving: Vec<usize> = walk
            .iter()
            .filter_map(|s| match s {
                Step::Edge(i, _) => Some(*i),
                Step::Pole(..) | Step::Seam(..) => None,
            })
            .collect();
        if kept || moving.is_empty() {
            walked.push(None);
            continue;
        }
        // A ring of one edge has no joint to fold, and as a hole may still
        // stand on another branch from its outer ring.
        if walk.len() < 2 {
            walked.push(Some(moving));
            continue;
        }
        let ends_of = |step: Step| -> OgeomResult<(Point2, Point2)> {
            let (pc, range, back) = match step {
                Step::Edge(i, back) => (&edges[i].pcurve, edges[i].prange, back),
                Step::Pole(j, back) => (&poles[j].pcurve, poles[j].prange, back),
                Step::Seam(k, back) => (&sides[k].1, sides[k].2, back),
            };
            let (a, b) = (pc.point_at(range.0, tol)?, pc.point_at(range.1, tol)?);
            Ok(if back { (b, a) } else { (a, b) })
        };
        let radius = |step: Step| match step {
            Step::Edge(i, _) => edges[i].tolerance.max(edges[i].ends_tolerance),
            Step::Seam(k, _) => edges[sides[k].0]
                .tolerance
                .max(edges[sides[k].0].ends_tolerance),
            Step::Pole(..) => 0.0,
        };
        let is_pole = |step: Step| matches!(step, Step::Pole(..));
        // Whether a joint meets: on the surface, within both edges' radii,
        // or beside a degenerate edge, whose line stands for one point, in
        // the chart within the same reach. There the chart shrinks to
        // nothing along the line, and the arrangement welds a strand's end
        // there at its length in space, unscaled, no finer than its own
        // node snap: the fold carries an edge onto the pole's line only
        // where that weld then meets it.
        let meets = |a: Point2, b: Point2, step: Step, prev: Step, within: f64| -> bool {
            if is_pole(step) || is_pole(prev) {
                a.distance(b) <= within.max(chart_snap)
            } else {
                apart(a, b).is_some_and(|d| d <= within)
            }
        };
        // Whether a joint's two ends, apart in the chart, stand across a
        // degenerate row: the surface does not move along the chart segment
        // between them, its derivative along it vanishing at both ends and
        // in the middle, so that the whole gap moves it by under the
        // confusion distance. A small surface's ends a period apart lie
        // within the edges' radii without its derivative vanishing, and are
        // not such a joint.
        let across_a_degenerate_row = |a: Point2, b: Point2| -> bool {
            let gap = b - a;
            let length = gap.magnitude();
            if length <= chart_snap {
                return false;
            }
            let middle = Point2::new(f64::midpoint(a.x, b.x), f64::midpoint(a.y, b.y));
            [a, middle, b].iter().all(|p| {
                surface
                    .d1_at(p.x, p.y, tol)
                    .is_ok_and(|(du, dv)| (du * gap.x + dv * gap.y).magnitude() <= tol.confusion())
            })
        };
        let anchor = |step: &Step| matches!(step, Step::Pole(..) | Step::Seam(..));
        let start = if let Some(k) = walk.iter().position(|s| is_pole(*s)) {
            k
        } else if let Some(k) = walk.iter().position(anchor) {
            k
        } else {
            let mut start = 0;
            for k in 0..walk.len() {
                let (step, prev) = (walk[k], walk[(k + walk.len() - 1) % walk.len()]);
                let (end, from) = (ends_of(prev)?.1, ends_of(step)?.0);
                if across_a_degenerate_row(end, from) {
                    start = k;
                    break;
                }
            }
            start
        };
        walk.rotate_left(start);
        let (begin, mut last) = ends_of(walk[0])?;
        let mut moves: Vec<(usize, PlanarCurve)> = Vec::new();
        let mut broken = false;
        for k in 1..walk.len() {
            let (step, prev) = (walk[k], walk[k - 1]);
            let (from, to) = ends_of(step)?;
            let within = 2.0 * (radius(step) + radius(prev)) + tol.confusion();
            let Step::Edge(i, _) = step else {
                // A degenerate edge stays where it is stated, and a seam
                // walked once on the side walked: the ring must reach it
                // there.
                broken |= !meets(last, from, step, prev, within);
                last = to;
                continue;
            };
            let shift = whole(last, from);
            let moved = Point2::new(from.x + shift.x, from.y + shift.y);
            // A whole period between two ends is a jump of branch, unless
            // the joint stands across a degenerate row; beside a degenerate
            // edge it is read in the chart.
            let jumps = if is_pole(step) || is_pole(prev) {
                !meets(last, from, step, prev, within)
            } else {
                !across_a_degenerate_row(last, from)
            };
            let fold = (shift.x != 0.0 || shift.y != 0.0)
                && jumps
                && meets(last, moved, step, prev, within);
            last = to;
            if fold {
                let by = ogeom_math::Transform2::translation(shift);
                let pcurve = edges[i].pcurve.transformed(&by, tol)?;
                if inside(&pcurve, edges[i].prange)? {
                    moves.push((i, pcurve));
                    last = Point2::new(to.x + shift.x, to.y + shift.y);
                }
            }
        }
        // Folded, the ring closes where it began, or it winds round the
        // surface and stays as stated.
        let (first, closing) = (walk[0], walk[walk.len() - 1]);
        let within = 2.0 * (radius(first) + radius(closing)) + tol.confusion();
        if !broken && !moves.is_empty() && meets(last, begin, first, closing, within) {
            for (i, pcurve) in moves {
                edges[i].pcurve = pcurve;
            }
            // Read on the side walked, the seam bounds the folded ring on
            // that column alone.
            for (i, pcurve, range) in sides {
                edges[i].pcurve = pcurve;
                edges[i].prange = range;
                edges[i].other_side = None;
            }
        }
        walked.push((!walk.iter().any(anchor)).then_some(moving));
    }
    // A hole on another branch from its outer ring is brought onto it.
    let extent = |edges: &[BoundaryEdge], ring: &[usize]| -> OgeomResult<(Point2, Point2)> {
        let mut low = Point2::new(f64::INFINITY, f64::INFINITY);
        let mut high = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
        for &i in ring {
            let e = &edges[i];
            for k in 0..=4 {
                let t = e.prange.0 + (e.prange.1 - e.prange.0) * f64::from(k) / 4.0;
                let p = e.pcurve.point_at(t, tol)?;
                low = Point2::new(low.x.min(p.x), low.y.min(p.y));
                high = Point2::new(high.x.max(p.x), high.y.max(p.y));
            }
        }
        Ok((low, high))
    };
    let Some(outer) = wires.first() else {
        return Ok(());
    };
    let outer: Vec<usize> = outer
        .iter()
        .filter_map(|e| index.get(&EdgeKey::of(e)).map(|&(i, _)| i))
        .collect();
    if outer.is_empty() {
        return Ok(());
    }
    let (olow, ohigh) = extent(edges, &outer)?;
    for hole in walked.iter().skip(1).flatten() {
        let (low, high) = extent(edges, hole)?;
        let centre = Point2::new(f64::midpoint(low.x, high.x), f64::midpoint(low.y, high.y));
        let ocentre = Point2::new(
            f64::midpoint(olow.x, ohigh.x),
            f64::midpoint(olow.y, ohigh.y),
        );
        let shift = whole(ocentre, centre);
        let beside =
            |lo: f64, hi: f64, olo: f64, ohi: f64| hi < olo - chart_snap || lo > ohi + chart_snap;
        let within =
            |lo: f64, hi: f64, olo: f64, ohi: f64| lo >= olo - chart_snap && hi <= ohi + chart_snap;
        let along = |s: f64, lo: f64, hi: f64, olo: f64, ohi: f64| {
            s == 0.0 || (beside(lo, hi, olo, ohi) && within(lo + s, hi + s, olo, ohi))
        };
        if (shift.x == 0.0 && shift.y == 0.0)
            || !along(shift.x, low.x, high.x, olow.x, ohigh.x)
            || !along(shift.y, low.y, high.y, olow.y, ohigh.y)
        {
            continue;
        }
        let by = ogeom_math::Transform2::translation(shift);
        let mut moved = Vec::with_capacity(hole.len());
        for &i in hole {
            let pcurve = edges[i].pcurve.transformed(&by, tol)?;
            if !inside(&pcurve, edges[i].prange)? {
                moved.clear();
                break;
            }
            moved.push((i, pcurve));
        }
        for (i, pcurve) in moved {
            edges[i].pcurve = pcurve;
        }
    }
    Ok(())
}

/// Which side of a seam an occurrence walks: the one that keeps the face
/// on the ring's left, as every ring of a face stored forward keeps it. A
/// seam's two sides bound the window the face wraps, so the face lies
/// towards the other side, and the occurrence runs along the side from
/// which that lies on its left. `None` where the sides do not stand apart
/// across the occurrence's run.
pub(crate) fn walked_seam_side(
    e: &BoundaryEdge,
    back: bool,
    tol: Tolerances,
) -> OgeomResult<Option<(PlanarCurve, (f64, f64))>> {
    let Some((other, orange)) = &e.other_side else {
        return Ok(None);
    };
    let (a, b) = (
        e.pcurve.point_at(e.prange.0, tol)?,
        e.pcurve.point_at(e.prange.1, tol)?,
    );
    let run = if back { a - b } else { b - a };
    let mid = |pc: &PlanarCurve, r: (f64, f64)| pc.point_at(f64::midpoint(r.0, r.1), tol);
    let towards = mid(other, *orange)? - mid(&e.pcurve, e.prange)?;
    let left = run.x * towards.y - run.y * towards.x;
    Ok(if left > 0.0 {
        Some((e.pcurve.clone(), e.prange))
    } else if left < 0.0 {
        Some((other.clone(), *orange))
    } else {
        None
    })
}

/// A ring walking a seam both ways brought onto the seam's branch of a
/// closed chart: each seam occurrence on the side the ring walks (see
/// [`walked_seam_side`]), and every other edge moved by the whole periods
/// that make it start where the ring has reached. A file may state a rim
/// a period or more off the window its seam bounds; the arrangement
/// draws both seam columns and needs the rims between them.
///
/// The ring is moved only where, so walked, every joint meets and it
/// closes where it began, and every moved pcurve stays on a stretch the
/// chart covers. Consecutive occurrences of one seam meet across the
/// closure, a whole period apart, where the ring crosses the join. A ring
/// with any other edge walked twice is left as stated.
pub(crate) fn fold_seamed_rings(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    edges: &mut [BoundaryEdge],
    (poles, pole_keys): (&[PoleEdge], &[EdgeKey]),
    slits: &ogeom_core::FastSet<EdgeKey>,
    tol: Tolerances,
) -> OgeomResult<()> {
    let chart_snap = param_snap(surface, tol);
    #[derive(Clone)]
    enum Step {
        Seam(usize, Point2, Point2),
        Edge(usize, Point2, Point2),
        Pole(Point2, Point2),
    }
    let ((ua, ub), (va, vb)) = surface.domain();
    let periods = [
        ((surface.is_periodic_u() || surface.is_closed_u(tol)) && ub > ua).then_some(ub - ua),
        ((surface.is_periodic_v() || surface.is_closed_v(tol)) && vb > va).then_some(vb - va),
    ];
    if periods.iter().all(Option::is_none) {
        return Ok(());
    }
    let periodic = [surface.is_periodic_u(), surface.is_periodic_v()];
    let whole = |to: Point2, from: Point2| -> ogeom_math::Vector2 {
        let mut shift = ogeom_math::Vector2::new(0.0, 0.0);
        if let Some(period) = periods[0] {
            shift.x = ((to.x - from.x) / period).round() * period;
        }
        if let Some(period) = periods[1] {
            shift.y = ((to.y - from.y) / period).round() * period;
        }
        shift
    };
    let inside = |pc: &PlanarCurve, range: (f64, f64)| -> OgeomResult<bool> {
        for k in 0..=2 {
            let p = pc.point_at(range.0 + (range.1 - range.0) * f64::from(k) / 2.0, tol)?;
            if (!periodic[0] && (p.x < ua - chart_snap || p.x > ub + chart_snap))
                || (!periodic[1] && (p.y < va - chart_snap || p.y > vb + chart_snap))
            {
                return Ok(false);
            }
        }
        Ok(true)
    };
    let apart = |a: Point2, b: Point2| -> Option<f64> {
        let mid = Point2::new(f64::midpoint(a.x, b.x), f64::midpoint(a.y, b.y));
        let (p, m, q) = (
            surface.point_at(a.x, a.y, tol).ok()?,
            surface.point_at(mid.x, mid.y, tol).ok()?,
            surface.point_at(b.x, b.y, tol).ok()?,
        );
        Some(p.distance(m).max(m.distance(q)))
    };
    // A seam's two sides stand one period apart across the direction its
    // chart closes in; a file may state the far side whole periods
    // further. On a chart that closes without repeating, only one of the
    // two copies a period from the near side lies in the window: the far
    // side is moved there, exactly, by whole periods.
    for e in edges.iter_mut() {
        let Some((other, orange)) = &e.other_side else {
            continue;
        };
        let mid = |pc: &PlanarCurve, r: (f64, f64)| pc.point_at(f64::midpoint(r.0, r.1), tol);
        let (near, far) = (mid(&e.pcurve, e.prange)?, mid(other, *orange)?);
        let mut shift = ogeom_math::Vector2::new(0.0, 0.0);
        for (axis, period) in periods.iter().enumerate() {
            let Some(period) = *period else { continue };
            if periodic[axis] {
                continue;
            }
            let (from, to, lo, hi) = if axis == 0 {
                (near.x, far.x, ua, ub)
            } else {
                (near.y, far.y, va, vb)
            };
            let turns = ((to - from) / period).round();
            if turns.abs() < 2.0 || ((to - from) - turns * period).abs() > chart_snap {
                continue;
            }
            let inside = |x: f64| x >= lo - chart_snap && x <= hi + chart_snap;
            let copies: Vec<f64> = [-1.0, 1.0]
                .into_iter()
                .filter(|s| inside(from + s * period))
                .collect();
            if let [s] = copies.as_slice() {
                let by = (from + s * period) - (from + turns * period);
                if axis == 0 {
                    shift.x = by;
                } else {
                    shift.y = by;
                }
            }
        }
        if shift.x != 0.0 || shift.y != 0.0 {
            let by = ogeom_math::Transform2::translation(shift);
            let moved = other.transformed(&by, tol)?;
            let range = *orange;
            e.other_side = Some((moved, range));
        }
    }
    let mut uses: ogeom_core::FastMap<EdgeKey, (usize, usize)> = ogeom_core::FastMap::default();
    for (i, e) in edges.iter().enumerate() {
        uses.insert(e.node, (i, 0));
    }
    let wires = rings_without(
        model,
        &face.oriented(ogeom_topo::Orientation::Forward),
        slits,
    )?;
    for edge in wires.iter().flatten() {
        if let Some(slot) = uses.get_mut(&EdgeKey::of(edge)) {
            slot.1 += 1;
        }
    }
    'rings: for ring in &wires {
        let mut walk: Vec<Step> = Vec::with_capacity(ring.len());
        for edge in ring {
            let back = edge.orientation() == ogeom_topo::Orientation::Reversed;
            let key = EdgeKey::of(edge);
            let ends = |pc: &PlanarCurve, r: (f64, f64)| -> OgeomResult<(Point2, Point2)> {
                let (a, b) = (pc.point_at(r.0, tol)?, pc.point_at(r.1, tol)?);
                Ok(if back { (b, a) } else { (a, b) })
            };
            if let Some(&(i, n)) = uses.get(&key) {
                let e = &edges[i];
                if e.other_side.is_some() && n == 2 {
                    let Some((side, range)) = walked_seam_side(e, back, tol)? else {
                        continue 'rings;
                    };
                    let (a, b) = ends(&side, range)?;
                    walk.push(Step::Seam(i, a, b));
                } else if n == 1 && e.other_side.is_none() && e.seam.is_none() {
                    let (a, b) = ends(&e.pcurve, e.prange)?;
                    walk.push(Step::Edge(i, a, b));
                } else {
                    continue 'rings;
                }
            } else if let Some(j) = pole_keys.iter().position(|k| *k == key) {
                let (a, b) = ends(&poles[j].pcurve, poles[j].prange)?;
                walk.push(Step::Pole(a, b));
            } else {
                continue 'rings;
            }
        }
        let Some(start) = walk.iter().position(|s| matches!(s, Step::Seam(..))) else {
            continue;
        };
        walk.rotate_left(start);
        let radius = |s: &Step| match s {
            Step::Seam(i, ..) | Step::Edge(i, ..) => {
                edges[*i].tolerance.max(edges[*i].ends_tolerance)
            }
            Step::Pole(..) => 0.0,
        };
        let meets = |a: Point2, b: Point2, step: &Step, prev: &Step| -> bool {
            if matches!(step, Step::Pole(..)) || matches!(prev, Step::Pole(..)) {
                a.distance(b) <= SEAM_HAIR
            } else {
                let within = 2.0 * (radius(step) + radius(prev)) + tol.confusion();
                a.distance(b) <= chart_snap || apart(a, b).is_some_and(|d| d <= within)
            }
        };
        let Step::Seam(_, begin, mut last) = walk[0] else {
            continue;
        };
        let mut moves: Vec<(usize, PlanarCurve)> = Vec::new();
        for k in 1..walk.len() {
            let (step, prev) = (&walk[k], &walk[k - 1]);
            match step {
                Step::Edge(i, from, to) => {
                    let shift = whole(last, *from);
                    let moved = Point2::new(from.x + shift.x, from.y + shift.y);
                    if !meets(last, moved, step, prev) {
                        continue 'rings;
                    }
                    if shift.x != 0.0 || shift.y != 0.0 {
                        let by = ogeom_math::Transform2::translation(shift);
                        let pcurve = edges[*i].pcurve.transformed(&by, tol)?;
                        if !inside(&pcurve, edges[*i].prange)? {
                            continue 'rings;
                        }
                        moves.push((*i, pcurve));
                    }
                    last = Point2::new(to.x + shift.x, to.y + shift.y);
                }
                Step::Seam(i, from, to) => {
                    let across = matches!(prev, Step::Seam(j, ..) if j == i) && {
                        let shift = whole(last, *from);
                        (shift.x != 0.0 || shift.y != 0.0)
                            && meets(
                                last,
                                Point2::new(from.x + shift.x, from.y + shift.y),
                                step,
                                prev,
                            )
                    };
                    if !across && !meets(last, *from, step, prev) {
                        continue 'rings;
                    }
                    last = *to;
                }
                Step::Pole(from, to) => {
                    if !meets(last, *from, step, prev) {
                        continue 'rings;
                    }
                    last = *to;
                }
            }
        }
        let closing = &walk[walk.len() - 1];
        let across = matches!((&walk[0], closing), (Step::Seam(i, ..), Step::Seam(j, ..)) if i == j)
            && {
                let shift = whole(last, begin);
                (shift.x != 0.0 || shift.y != 0.0)
                    && meets(
                        last,
                        Point2::new(begin.x + shift.x, begin.y + shift.y),
                        &walk[0],
                        closing,
                    )
            };
        if !across && !meets(last, begin, &walk[0], closing) {
            continue;
        }
        for (i, pcurve) in moves {
            edges[i].pcurve = pcurve;
        }
    }
    Ok(())
}

/// The stretches of a pole row a face's ring crosses with no edge there:
/// two edges meeting at a cone's apex or a sphere's pole stand apart in the
/// chart, at the columns they arrive along, and the row between them is
/// one point in space. Each such joint becomes a pole of the face, the
/// chart line from one edge's end to the next one's start, as a degenerate
/// edge there would state it.
///
/// A joint closes this way only where that line is the pole row the face
/// itself bounds:
/// - both ends lift to one point, within the joint vertex's radius;
/// - the row is that point: the surface's rate along it, read at its ends
///   and middle, times the chart's whole width that way, and its rate
///   across times the line's offset, stay within that radius;
/// - the line is shorter than the chart's period that way, and the face
///   lies on its left, as every ring of a face stored forward keeps it:
///   towards the middle of the face's own trim.
///
/// Any other joint is left for the arrangement to snap or refuse by name.
/// A ring that walks a seam both ways reads each occurrence on the side it
/// walks ([`walked_seam_side`]).
pub(crate) fn close_pole_joints(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    (edges, slits): (&[BoundaryEdge], &ogeom_core::FastSet<EdgeKey>),
    poles: &mut Vec<PoleEdge>,
    pole_keys: &[EdgeKey],
    tol: Tolerances,
) -> OgeomResult<()> {
    let chart_snap = param_snap(surface, tol);
    let index: ogeom_core::FastMap<EdgeKey, usize> =
        edges.iter().enumerate().map(|(i, e)| (e.node, i)).collect();
    // The middle of the face's trim in the chart.
    let (mut low, mut high) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for e in edges {
        for k in 0..=4 {
            let p = e.pcurve.point_at(
                e.prange.0 + (e.prange.1 - e.prange.0) * f64::from(k) / 4.0,
                tol,
            )?;
            low = Point2::new(low.x.min(p.x), low.y.min(p.y));
            high = Point2::new(high.x.max(p.x), high.y.max(p.y));
        }
    }
    let centre = Point2::new(f64::midpoint(low.x, high.x), f64::midpoint(low.y, high.y));
    let ((u0, u1), (v0, v1)) = surface.domain();
    let closes = [
        surface.is_periodic_u() || surface.is_closed_u(tol),
        surface.is_periodic_v() || surface.is_closed_v(tol),
    ];
    for ring in rings_without(
        model,
        &face.oriented(ogeom_topo::Orientation::Forward),
        slits,
    )? {
        // Each occurrence's edge, whether the ring walks it backwards and,
        // for a seam the ring walks both ways, the side it walks; `None`
        // for a degenerate edge, which stands between its neighbours
        // already.
        type Walked = (usize, bool, Option<(PlanarCurve, (f64, f64))>);
        let mut walk: Vec<Option<Walked>> = Vec::with_capacity(ring.len());
        for edge in &ring {
            let key = EdgeKey::of(edge);
            if pole_keys.contains(&key) {
                walk.push(None);
                continue;
            }
            let Some(&i) = index.get(&key) else {
                walk.clear();
                break;
            };
            let back = edge.orientation() == ogeom_topo::Orientation::Reversed;
            let side = if edges[i].other_side.is_some() {
                let Some(side) = walked_seam_side(&edges[i], back, tol)? else {
                    walk.clear();
                    break;
                };
                Some(side)
            } else {
                None
            };
            walk.push(Some((i, back, side)));
        }
        let n = walk.len();
        for k in 0..n {
            let (Some((a, a_back, a_side)), Some((b, b_back, b_side))) =
                (&walk[k], &walk[(k + 1) % n])
            else {
                continue;
            };
            // The joint: where the ring leaves one edge and takes up the
            // next, in the chart and in space, with the joint vertex's
            // radius as each edge states it.
            let end = |i: usize,
                       back: bool,
                       side: &Option<(PlanarCurve, (f64, f64))>,
                       last: bool|
             -> OgeomResult<(Point2, (Point, f64))> {
                let e = &edges[i];
                let (pcurve, range) = side
                    .as_ref()
                    .map_or((&e.pcurve, e.prange), |(c, r)| (c, *r));
                let at_start = back == last;
                let t = if at_start { range.0 } else { range.1 };
                Ok((pcurve.point_at(t, tol)?, e.ends[usize::from(!at_start)]))
            };
            let ((from, (pa, ra)), (to, (pb, rb))) = (
                end(*a, *a_back, a_side, true)?,
                end(*b, *b_back, b_side, false)?,
            );
            if from.distance(to) <= chart_snap {
                continue;
            }
            let reach = ra.max(rb).max(tol.confusion());
            if pa.distance(pb) > reach {
                continue;
            }
            let along_u = (to.x - from.x).abs() >= (to.y - from.y).abs();
            let (length, offset, width, closing) = if along_u {
                (
                    (to.x - from.x).abs(),
                    (to.y - from.y).abs(),
                    (u1 - u0).abs(),
                    closes[0],
                )
            } else {
                (
                    (to.y - from.y).abs(),
                    (to.x - from.x).abs(),
                    (v1 - v0).abs(),
                    closes[1],
                )
            };
            if closing && length >= width - chart_snap {
                continue;
            }
            let mut collapses = true;
            for s in [0.0, 0.5, 1.0] {
                let q = Point2::new(from.x + (to.x - from.x) * s, from.y + (to.y - from.y) * s);
                let Ok((du, dv)) = surface.d1_at(q.x, q.y, tol) else {
                    collapses = false;
                    break;
                };
                let (rate, across) = if along_u {
                    (du.magnitude(), dv.magnitude())
                } else {
                    (dv.magnitude(), du.magnitude())
                };
                if !(rate * width <= reach && across * offset <= reach) {
                    collapses = false;
                    break;
                }
            }
            if !collapses {
                continue;
            }
            // The face on the line's left.
            let run = to - from;
            let middle = Point2::midpoint(from, to);
            let towards = centre - middle;
            if run.x * towards.y - run.y * towards.x <= 0.0 {
                continue;
            }
            let Ok(at) = surface.point_at(from.x, from.y, tol) else {
                continue;
            };
            let line = ogeom_geom::Line2d::segment(from, to, tol)?;
            let prange = ogeom_geom::Curve2d::domain(&line);
            poles.push(PoleEdge {
                point: at,
                pcurve: line.into(),
                prange,
            });
        }
    }
    Ok(())
}

/// A face's pcurves moved by whole periods, all alike, so the middle of
/// their extent lies in the chart's own window along each periodic
/// direction.
///
/// A periodic face's trim may sit any number of turns along its chart: a
/// partial drum's boundary can arrive two turns below zero. The paving and
/// the arrangement fold a section's image into the chart and try it at
/// most one period either side, and against a trim two turns off no try
/// lands inside it. Moved home, the trim straddles the window's edge at
/// most, which one period's try reaches. The face is the same face: every
/// pcurve moves by the same whole period. Whether they moved.
pub(crate) fn bring_trim_home(
    surface: &SurfaceGeometry,
    edges: &mut [BoundaryEdge],
    poles: &mut [PoleEdge],
    tol: Tolerances,
) -> OgeomResult<bool> {
    if !(surface.is_periodic_u() || surface.is_periodic_v()) {
        return Ok(false);
    }
    let mut low = Point2::new(f64::INFINITY, f64::INFINITY);
    let mut high = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut take = |pc: &PlanarCurve, range: (f64, f64)| -> OgeomResult<()> {
        for k in 0..=4 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 4.0;
            let p = pc.point_at(t, tol)?;
            low = Point2::new(low.x.min(p.x), low.y.min(p.y));
            high = Point2::new(high.x.max(p.x), high.y.max(p.y));
        }
        Ok(())
    };
    for e in edges.iter() {
        take(&e.pcurve, e.prange)?;
        if let Some((other, range)) = &e.other_side {
            take(other, *range)?;
        }
    }
    for p in poles.iter() {
        take(&p.pcurve, p.prange)?;
    }
    let ((ua, ub), (va, vb)) = surface.domain();
    // A trim within a period of the window is reached already, and stays
    // where it was stated.
    let turns = |periodic: bool, lo: f64, hi: f64, start: f64, end: f64| -> f64 {
        let span = end - start;
        if !periodic
            || span <= 0.0
            || !lo.is_finite()
            || !hi.is_finite()
            || (lo >= start - span && hi <= end + span)
        {
            return 0.0;
        }
        -((f64::midpoint(lo, hi) - start) / span).floor() * span
    };
    let shift = ogeom_math::Vector2::new(
        turns(surface.is_periodic_u(), low.x, high.x, ua, ub),
        turns(surface.is_periodic_v(), low.y, high.y, va, vb),
    );
    if shift.x == 0.0 && shift.y == 0.0 {
        return Ok(false);
    }
    let moved = ogeom_math::Transform2::translation(shift);
    for e in edges.iter_mut() {
        e.pcurve = e.pcurve.transformed(&moved, tol)?;
        if let Some((other, _)) = &mut e.other_side {
            *other = other.transformed(&moved, tol)?;
        }
        if let Some((forward, reversed)) = &mut e.seam {
            *forward = forward.transformed(&moved, tol)?;
            *reversed = reversed.transformed(&moved, tol)?;
        }
    }
    for p in poles.iter_mut() {
        p.pcurve = p.pcurve.transformed(&moved, tol)?;
    }
    Ok(true)
}

/// A parameter brought onto the turn its edge actually covers.
///
/// A periodic curve's own domain and an edge's range over it need not agree
/// on which turn to count from: a sphere's seam runs its half meridian over
/// `[-π/2, π/2]`, while every crossing found on it comes back in `[0, 2π)`.
/// Left alone, a crossing on the negative half arrives a turn on, sits
/// outside the range, and is discarded, so the seam never splits where a
/// section genuinely meets it, and the arrangement finds the chain hanging.
pub(crate) fn onto_range(t: f64, curve: &Curve, range: (f64, f64), tol: Tolerances) -> f64 {
    if !curve.is_periodic() {
        return t;
    }
    let (a, b) = curve.domain();
    let period = b - a;
    if period <= 0.0 {
        return t;
    }
    for k in [0.0, -1.0, 1.0, -2.0, 2.0] {
        let shifted = period.mul_add(k, t);
        if shifted >= range.0 - tol.parametric() && shifted <= range.1 + tol.parametric() {
            return shifted;
        }
    }
    t
}

/// A stretch of an overlap on the first curve, and the second curve's
/// parameter at either end where the second curve's range ends it.
pub(crate) type OverlapPart = ((f64, f64), [Option<f64>; 2]);

/// The parts of an overlap that fall within the second curve's own bounded
/// range, stated in the *first* curve's parameter, in increasing order,
/// each with the second curve's parameter at either end where the range's
/// own end, not the overlap's, bounds it.
///
/// An overlap is between two curves. An edge covers only part of its curve.
/// The correspondence the overlap states is affine, so the edge's range
/// carries across as an interval, and on a periodic curve it carries across
/// up to whole turns. A full turn of the overlap meets an edge whose range
/// starts elsewhere on the circle in two stretches, one either side of the
/// edge's own start, and both are returned. Empty where the edge's own
/// stretch and the overlap do not meet.
pub(crate) fn overlap_within(
    overlap: &ogeom_intersect::Overlap,
    range: (f64, f64),
    curve: &Curve,
    tol: Tolerances,
) -> Vec<OverlapPart> {
    let ordered = |r: (f64, f64)| if r.0 <= r.1 { r } else { (r.1, r.0) };
    let span_a = overlap.on_a.1 - overlap.on_a.0;
    let span_b = overlap.on_b.1 - overlap.on_b.0;
    if span_b.abs() <= f64::MIN_POSITIVE {
        return Vec::new();
    }
    let to_a = |t: f64| overlap.on_a.0 + span_a * (t - overlap.on_b.0) / span_b;
    let (from_start, from_end) = (to_a(range.0), to_a(range.1));
    let ((wlo, elo), (whi, ehi)) = if from_start <= from_end {
        ((from_start, range.0), (from_end, range.1))
    } else {
        ((from_end, range.1), (from_start, range.0))
    };
    let (lo, hi) = ordered(overlap.on_a);
    let period = if curve.is_periodic() {
        let (a, b) = curve.domain();
        // The turn measured in the *first* curve's parameter, which is where
        // the intersection is being taken.
        (b - a) * (span_a / span_b).abs()
    } else {
        0.0
    };
    let turns: &[f64] = if period == 0.0 {
        &[0.0]
    } else {
        &[-2.0, -1.0, 0.0, 1.0, 2.0]
    };
    turns
        .iter()
        .map(|&k| {
            let (start, end) = (period.mul_add(k, wlo), period.mul_add(k, whi));
            (
                (lo.max(start), hi.min(end)),
                [(start > lo).then_some(elo), (end < hi).then_some(ehi)],
            )
        })
        .filter(|((a, b), _)| b - a > tol.parametric())
        .collect()
}

/// A parameter carried proportionally from one range to another.
pub(crate) fn rescale(t: f64, from: (f64, f64), to: (f64, f64)) -> f64 {
    let span = from.1 - from.0;
    if span.abs() <= f64::MIN_POSITIVE {
        return to.0;
    }
    to.0 + (to.1 - to.0) * (t - from.0) / span
}

/// Whether a curve is exact: a line or a conic, or a trim of one (an arc
/// stated as a stretch of its circle is as exact as its basis).
pub(crate) fn exact(c: &Curve) -> bool {
    match c {
        Curve::Line(_)
        | Curve::Circle(_)
        | Curve::Ellipse(_)
        | Curve::Hyperbola(_)
        | Curve::Parabola(_) => true,
        Curve::Trimmed(t) => exact(t.basis()),
        _ => false,
    }
}

/// How far `curve`, held to `tolerance`, may stand off the edge `e` and
/// still be read as running along it.
///
/// A fitted curve wobbles about the edge it runs along by far more than its
/// stated budget, and is measured loosely. Two exact curves are held to
/// their own tolerances, as an exact section is: an edge that parts from
/// another by a few microns over its length (a wall leaning off a pad's
/// side) runs beside it, not along it, and read as along it the strands
/// that meet it are left a few microns short.
pub(crate) fn measuring_width(
    curve: &Curve,
    tolerance: f64,
    e: &BoundaryEdge,
    tol: Tolerances,
) -> f64 {
    let floor = if exact(curve) && exact(&e.curve) {
        tol.confusion() * 10.0
    } else {
        tol.confusion() * 1e3
    };
    (tolerance.max(e.tolerance) * 2.0).max(floor)
}

/// Whether the distance `off` grows in proportion from where two curves
/// meet, at `meet` (a parameter and the distance there), out to `t`: the
/// two crossing at an angle there, not parting tangentially.
///
/// Growing as the k-th power from the meeting, the distance halfway out
/// stands at 2^-k of the distance at the end: past 2^-3/2, half way between
/// the crossing's 1/2 and the tangent parting's 1/4, the two cross.
pub(crate) fn grows_in_proportion(
    off: impl Fn(f64) -> OgeomResult<f64>,
    meet: (f64, f64),
    t: f64,
) -> OgeomResult<bool> {
    let rise = off(t)? - meet.1;
    let half = off(f64::midpoint(meet.0, t))? - meet.1;
    Ok(rise > 0.0 && half >= rise * core::f64::consts::FRAC_1_SQRT_2 / 2.0)
}

/// Whether `curve`, over `crange` and held to `tolerance`, crosses the edge
/// `e` at a shallow turn where the two meet at `on_a` on it, rather than
/// touching it.
///
/// Crossing at a turn whose sine is `s`, the two stand within the measuring
/// width `w` of each other for `w / s` either side of the crossing, and
/// past it their distance grows in proportion on both sides. Touching, it
/// grows as the square of the way from the touch, and a crossing found
/// along the touch is noise in the valley between them: from there the
/// distance grows on at most one side in proportion. Bending apart by `k`,
/// the difference of their curvature vectors, the two part by `k x² / 2`
/// over a way `x` besides the turn's `s x`, which over `w / s` stays below
/// half the width only where `s² > k w`: a shallower turn is not told from
/// a touch at the width, and reads as one.
pub(crate) fn crosses_shallowly(
    curve: &Curve,
    crange: (f64, f64),
    tolerance: f64,
    e: &BoundaryEdge,
    (on_a, on_b): (f64, f64),
    tol: Tolerances,
) -> OgeomResult<bool> {
    let bending =
        |c: &Curve, t: f64| -> OgeomResult<Option<(ogeom_math::Vector, ogeom_math::Vector)>> {
            let d = c.derivatives_at(t, 2, tol)?;
            let speed = d[1].magnitude();
            if speed <= tol.confusion() {
                return Ok(None);
            }
            let unit = d[1] / speed;
            let normal = d[2] - unit * unit.dot(d[2]);
            Ok(Some((d[1], normal / (speed * speed))))
        };
    let (Some((da, ka)), Some((db, kb))) = (bending(curve, on_a)?, bending(&e.curve, on_b)?) else {
        return Ok(false);
    };
    let (ma, mb) = (da.magnitude(), db.magnitude());
    let turn = da.cross(db).magnitude() / (ma * mb);
    let width = measuring_width(curve, tolerance, e, tol);
    if turn * turn <= (ka - kb).magnitude() * width {
        return Ok(false);
    }
    let step = width / turn / ma;
    let off = |t: f64| -> OgeomResult<f64> {
        distance_to_edge_curve(&e.curve, e.crange, curve.point_at(t, tol)?, tol)
    };
    let meet = (on_a, off(on_a)?);
    for t in [(on_a - step).max(crange.0), (on_a + step).min(crange.1)] {
        if (t - on_a).abs() <= tol.parametric() || !grows_in_proportion(off, meet, t)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The stretches of `curve` over `crange` that lie along the edge `e`
/// within the pair's honesty, as overlaps with the edge's own parameters at
/// their ends: the measured twin of the intersector's closed-form
/// coincidence, for a fitted curve on an exact one.
///
/// With `parted`, the sum of the radii both curves hold to, a fitted
/// stretch ends where the two part by it (see below); without, where they
/// part by the measuring width.
pub(crate) fn measured_overlaps(
    curve: &Curve,
    crange: (f64, f64),
    tolerance: f64,
    e: &BoundaryEdge,
    parted: Option<f64>,
    tol: Tolerances,
) -> OgeomResult<Vec<ogeom_intersect::Overlap>> {
    const SAMPLES: usize = 48;
    let width = measuring_width(curve, tolerance, e, tol);
    let mut near: Vec<bool> = Vec::with_capacity(SAMPLES + 1);
    for i in 0..=SAMPLES {
        #[allow(clippy::cast_precision_loss)]
        let t = crange.0 + (crange.1 - crange.0) * (i as f64) / (SAMPLES as f64);
        let p = curve.point_at(t, tol)?;
        near.push(
            e.bound.expanded(width).contains(p)
                && distance_to_edge_curve(&e.curve, e.crange, p, tol)? <= width,
        );
    }
    let mut out = Vec::new();
    let mut i = 0;
    while i <= SAMPLES {
        if !near[i] {
            i += 1;
            continue;
        }
        let start = i;
        while i < SAMPLES && near[i + 1] {
            i += 1;
        }
        let end = i;
        i += 1;
        // Three samples along is a stretch. Fewer is a crossing's blur.
        if end - start < 2 {
            continue;
        }
        #[allow(clippy::cast_precision_loss)]
        let at = |k: usize| crange.0 + (crange.1 - crange.0) * (k as f64) / (SAMPLES as f64);
        // Each end refined by bisection into the neighbouring sample gap.
        let refine = |inside: f64, outside: f64| -> OgeomResult<f64> {
            let (mut a, mut b) = (inside, outside);
            for _ in 0..24 {
                let m = f64::midpoint(a, b);
                let p = curve.point_at(m, tol)?;
                if e.bound.expanded(width).contains(p)
                    && distance_to_edge_curve(&e.curve, e.crange, p, tol)? <= width
                {
                    a = m;
                } else {
                    b = m;
                }
            }
            Ok(a)
        };
        let lo = if start == 0 {
            at(0)
        } else {
            refine(at(start), at(start - 1))?
        };
        let hi = if end == SAMPLES {
            at(SAMPLES)
        } else {
            refine(at(end), at(end + 1))?
        };
        // A stretch ending at one of the edge's own ends ends at that end's
        // vertex, not a measuring width past it: within the width of the
        // vertex the curve is still near the edge's end, and read as along
        // the edge there the contact splits a hair beyond the vertex the
        // face across it keeps.
        let at_vertex = |t: f64| -> OgeomResult<f64> {
            let p = curve.point_at(t, tol)?;
            for end in [e.crange.0, e.crange.1] {
                let vertex = e.curve.point_at(end, tol)?;
                if vertex.distance(p) > width * 2.0 {
                    continue;
                }
                let foot = ogeom_algo::project_on_curve(curve, vertex, 64, tol)?;
                let (clo, chi) = (crange.0.min(crange.1), crange.0.max(crange.1));
                if foot.distance <= width && foot.parameter >= clo && foot.parameter <= chi {
                    return Ok(foot.parameter);
                }
            }
            Ok(t)
        };
        let lo = if start == 0 { lo } else { at_vertex(lo)? };
        let hi = if end == SAMPLES { hi } else { at_vertex(hi)? };
        if hi - lo <= tol.parametric() {
            continue;
        }
        // Two exact curves that lie on one another do so over the whole of
        // what they share, so a stretch of them ends where one of them
        // does. One that ends inside both is two curves converging and
        // parting (lines meeting at a vertex at a tiny angle), near each
        // other there only by the measuring width, and is no overlap.
        let ends_there = |t: f64, own_end: bool| -> OgeomResult<bool> {
            if own_end {
                return Ok(true);
            }
            let p = curve.point_at(t, tol)?;
            for end in [e.crange.0, e.crange.1] {
                if e.curve.point_at(end, tol)?.distance(p) <= width * 2.0 {
                    return Ok(true);
                }
            }
            Ok(false)
        };
        let (ends_lo, ends_hi) = (ends_there(lo, start == 0)?, ends_there(hi, end == SAMPLES)?);
        // A fitted curve that ends on the edge and leaves it tangentially
        // from there, nearest it at its own end and apart by more than both
        // state halfway out, touches the edge at that end and runs beside
        // it: two surfaces tangent along a line, met by a third where the
        // line crosses it, cut the third in two curves tangent there (a
        // round's spandrel cut back by a drum, its edges on the drum along
        // the wall and along the round). Read as along the edge out to the
        // measuring width, it split the edge there and nothing split the
        // edge's twin on the curve's own solid.
        if !(exact(curve) && exact(&e.curve))
            && touches_from_own_end(
                curve,
                (lo, hi),
                (start == 0, end == SAMPLES),
                (ends_lo, ends_hi),
                tolerance,
                e,
                tol,
            )?
        {
            continue;
        }
        let (lo, hi) = if exact(curve) && exact(&e.curve) {
            if !(ends_lo && ends_hi) {
                continue;
            }
            (lo, hi)
        } else if ends_lo && ends_hi {
            (lo, hi)
        } else {
            // A fitted curve may come in along an edge tangentially and
            // leave the loose width inside both, its distance growing as the
            // square of the way from where they meet or slower. One that
            // leaves it growing in proportion, the two crossing at an angle
            // (a band's section across a drum, crossing the drum's own edge
            // at a shallow one), is near it only by the width: the stretch
            // is a crossing and no overlap.
            let off = |t: f64| -> OgeomResult<f64> {
                distance_to_edge_curve(&e.curve, e.crange, curve.point_at(t, tol)?, tol)
            };
            let mut meet = (lo, off(lo)?);
            for k in 1..=SAMPLES {
                #[allow(clippy::cast_precision_loss)]
                let t = lo + (hi - lo) * (k as f64) / (SAMPLES as f64);
                let d = off(t)?;
                if d < meet.1 {
                    meet = (t, d);
                }
            }
            let mut crossing = false;
            for (t, ends) in [(lo, ends_lo), (hi, ends_hi)] {
                if !ends {
                    crossing |= grows_in_proportion(off, meet, t)?;
                }
            }
            if crossing {
                continue;
            }
            (lo, hi)
        };
        if hi - lo <= tol.parametric() {
            continue;
        }
        // A fitted curve near an edge within the loose width may still
        // stand apart from it by more than the radii both hold to, and what
        // lies between them there is geometry, a lune the outline must
        // keep. Given those radii, the stretch is the runs along which the
        // two stand within them, each ended where they part by them: the
        // end on this curve lies within them of the other. A stretch that
        // never comes that near keeps the width's ends.
        let off = |t: f64| -> OgeomResult<f64> {
            distance_to_edge_curve(&e.curve, e.crange, curve.point_at(t, tol)?, tol)
        };
        let runs = match parted {
            Some(within) if !(exact(curve) && exact(&e.curve)) => {
                runs_within(&off, (lo, hi), within, SAMPLES)?.unwrap_or_else(|| vec![(lo, hi)])
            }
            _ => vec![(lo, hi)],
        };
        for (lo, hi) in runs {
            if hi - lo <= tol.parametric() {
                continue;
            }
            // An end carried to the other side is carried within the edge's
            // own window: a closed curve's two ends stand at one point, and the
            // edge covering one end of it does not cover the other.
            let on_edge = |t: f64| -> OgeomResult<f64> {
                let p = curve.point_at(t, tol)?;
                if parted.is_some() {
                    return Ok(nearest_on_edge_curve(&e.curve, e.crange, p, tol)?.1);
                }
                let foot = ogeom_algo::project_on_curve(&e.curve, p, 64, tol)?;
                Ok(onto_range(foot.parameter, &e.curve, e.crange, tol))
            };
            ogeom_core::diag!(
                Channel::Wire,
                "CONTACT along edge {} measured over ({lo:.6}, {hi:.6}) within {width:.2e}",
                e.node.index()
            );
            out.push(ogeom_intersect::Overlap {
                on_a: (lo, hi),
                on_b: (on_edge(lo)?, on_edge(hi)?),
            });
        }
    }
    // The samples are spaced along the whole curve, and an edge lying
    // along it over less than three spacings gives no stretch. Measured
    // from the edge's side instead: where every sample of the edge lies
    // within the width of the curve, the edge runs along it end to end,
    // and the stretch is the curve's between the feet of its two ends.
    if out.is_empty()
        && let Some(overlap) = edge_along_curve(curve, crange, e, width, tol)?
    {
        out.push(overlap);
    }
    Ok(out)
}

/// Whether the stretch `(lo, hi)` of the fitted `curve` near `e` starts
/// at one of the curve's own ends lying on `e` (inside it, or at one of
/// its ends, a vertex both leave tangentially) and ends inside both, the
/// two nearest at that own end and standing apart halfway along by more
/// than both curves state: a touch at the curve's end, not a run along
/// the edge (see [`measured_overlaps`]). A fitted copy of the edge that
/// ends with it stays within those radii, and is a run.
pub(crate) fn touches_from_own_end(
    curve: &Curve,
    (lo, hi): (f64, f64),
    (own_lo, own_hi): (bool, bool),
    (ends_lo, ends_hi): (bool, bool),
    tolerance: f64,
    e: &BoundaryEdge,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let (from, to) = if own_lo && !ends_hi {
        (lo, hi)
    } else if own_hi && !ends_lo {
        (hi, lo)
    } else {
        return Ok(false);
    };
    let off = |t: f64| -> OgeomResult<f64> {
        distance_to_edge_curve(&e.curve, e.crange, curve.point_at(t, tol)?, tol)
    };
    const STEPS: usize = 16;
    let at_own = off(from)?;
    for k in 1..=STEPS {
        #[allow(clippy::cast_precision_loss)]
        let t = from + (to - from) * (k as f64) / (STEPS as f64);
        if off(t)? < at_own {
            return Ok(false);
        }
    }
    Ok(off(f64::midpoint(from, to))? > tolerance + e.tolerance)
}

/// The radius an edge holds to where stretches end at it: its widest
/// measured on any face of it, or its own tolerance.
pub(crate) fn held_radius(
    held: &ogeom_core::FastMap<EdgeKey, (f64, &BoundaryEdge)>,
    node: EdgeKey,
    tolerance: f64,
) -> f64 {
    held.get(&node)
        .map_or(tolerance, |&(r, _)| r)
        .max(tolerance)
}

/// A shallow crossing that stands at one curve's end vertex.
pub(crate) enum End {
    /// At an end of the edge: the contact splits at the vertex's foot on
    /// it, the vertex standing `gap` off it, and the edge does not split.
    OfEdge { on_contact: f64, gap: f64 },
    /// At an end of the contact: the edge splits at its foot.
    OfContact { on_edge: f64, gap: f64 },
}

/// Where a shallow crossing of the contact and the edge `e`, at `on` on
/// each, stands at an end vertex of one of them.
///
/// Two curves crossing at a turn whose sine is `s` stand within `width`,
/// the sum of the radii they hold to, of each other for `width / s` either
/// side of the crossing, and the crossing is known along them only that
/// far. A vertex of either curve within that way of it, which the other
/// curve passes within `width` of, is where they cross: found once from
/// each pair of edges the vertex parts, the crossing is two points either
/// side of the vertex, and one curve would split twice where the other
/// does once.
pub(crate) fn shallow_at_end(
    contact: &ContactRec,
    e: &BoundaryEdge,
    on: (f64, f64),
    width: f64,
    tol: Tolerances,
) -> OgeomResult<Option<End>> {
    let (da, db) = (contact.curve.d1_at(on.0, tol)?, e.curve.d1_at(on.1, tol)?);
    let (ma, mb) = (da.magnitude(), db.magnitude());
    if ma <= tol.confusion() || mb <= tol.confusion() {
        return Ok(None);
    }
    let turn = da.cross(db).magnitude() / (ma * mb);
    if turn <= f64::EPSILON {
        return Ok(None);
    }
    let along = width / turn;
    let at = contact.curve.point_at(on.0, tol)?;
    let mut best: Option<(f64, End)> = None;
    for end in [e.crange.0, e.crange.1] {
        let vertex = e.curve.point_at(end, tol)?;
        let away = vertex.distance(at);
        if away > along {
            continue;
        }
        let (gap, t) = nearest_on_edge_curve(&contact.curve, contact.crange, vertex, tol)?;
        if gap <= width && best.as_ref().is_none_or(|b| away < b.0) {
            best = Some((away, End::OfEdge { on_contact: t, gap }));
        }
    }
    for end in [contact.crange.0, contact.crange.1] {
        let vertex = contact.curve.point_at(end, tol)?;
        let away = vertex.distance(at);
        if away > along {
            continue;
        }
        let (gap, t) = nearest_on_edge_curve(&e.curve, e.crange, vertex, tol)?;
        if gap <= width && best.as_ref().is_none_or(|b| away < b.0) {
            best = Some((away, End::OfContact { on_edge: t, gap }));
        }
    }
    Ok(best.map(|(_, end)| end))
}

/// The runs of `span` along which `off` stays within `within`, each end
/// bisected to where it rises past it, from `samples` even samples; `None`
/// where no sample comes that near.
pub(crate) fn runs_within(
    off: &impl Fn(f64) -> OgeomResult<f64>,
    span: (f64, f64),
    within: f64,
    samples: usize,
) -> OgeomResult<Option<Vec<(f64, f64)>>> {
    #[allow(clippy::cast_precision_loss)]
    let at = |k: usize| span.0 + (span.1 - span.0) * (k as f64) / (samples as f64);
    let near: Vec<bool> = (0..=samples)
        .map(|k| off(at(k)).map(|d| d <= within))
        .collect::<OgeomResult<_>>()?;
    if !near.contains(&true) {
        return Ok(None);
    }
    let parting = |inside: f64, outside: f64| -> OgeomResult<f64> {
        let (mut a, mut b) = (inside, outside);
        for _ in 0..24 {
            let m = f64::midpoint(a, b);
            if off(m)? <= within {
                a = m;
            } else {
                b = m;
            }
        }
        Ok(a)
    };
    let mut runs = Vec::new();
    let mut k = 0;
    while k <= samples {
        if !near[k] {
            k += 1;
            continue;
        }
        let first = k;
        while k < samples && near[k + 1] {
            k += 1;
        }
        let lo = if first == 0 {
            span.0
        } else {
            parting(at(first), at(first - 1))?
        };
        let hi = if k == samples {
            span.1
        } else {
            parting(at(k), at(k + 1))?
        };
        runs.push((lo, hi));
        k += 1;
    }
    Ok(Some(runs))
}

/// The stretches the contact shares with the target edge `e`, each ended
/// where the two curves part by the sum of the radii they hold to (within
/// it, each could be the other), measured once for
/// the pair and carried to whichever of the two is met as the contact
/// later.
///
/// Measured from each curve's side, a tangential parting ends on that
/// curve, and the two ends stand apart along the curves by as much as the
/// curves stay within that radius over, a way no vertex weld reaches. So
/// the end is measured once, on the curve held to the tighter radius
/// (its points lie nearest the truth), and carried onto the other by
/// projection: both sides split there.
pub(crate) fn shared_overlaps(
    contact: &ContactRec,
    e: &BoundaryEdge,
    held: &ogeom_core::FastMap<EdgeKey, (f64, &BoundaryEdge)>,
    shared: &mut ogeom_core::FastMap<(EdgeKey, EdgeKey), Vec<ogeom_intersect::Overlap>>,
    tol: Tolerances,
) -> OgeomResult<Vec<ogeom_intersect::Overlap>> {
    if let Some(found) = shared.get(&(contact.node, e.node)) {
        return Ok(found.clone());
    }
    let swap = |overlaps: &[ogeom_intersect::Overlap]| -> Vec<ogeom_intersect::Overlap> {
        overlaps
            .iter()
            .map(|o| ogeom_intersect::Overlap {
                on_a: o.on_b,
                on_b: o.on_a,
            })
            .collect()
    };
    let own = held.get(&contact.node).copied();
    let held_by_contact = held_radius(held, contact.node, contact.tolerance);
    let held_by_edge = held_radius(held, e.node, e.tolerance);
    let within = (held_by_contact + held_by_edge).max(tol.confusion());
    // The contact's own record, where it covers just what the contact does,
    // to be measured against from the edge's side.
    let owner = own
        .map(|(_, owner)| owner)
        .filter(|owner| owner.crange == contact.crange && *owner.curve == contact.curve);
    let from_contact = measured_overlaps(
        &contact.curve,
        contact.crange,
        contact.tolerance,
        e,
        Some(within),
        tol,
    )?;
    let from_edge = match owner {
        Some(owner) => swap(&measured_overlaps(
            &e.curve,
            e.crange,
            e.tolerance,
            owner,
            Some(within),
            tol,
        )?),
        None => Vec::new(),
    };
    // Each side is sampled along its own curve, and a stretch shorter than
    // three of one side's samples goes unseen from that side while the
    // other's may see it. A stretch seen from both is kept as the tighter
    // side measures it; one seen from one side only, as that side does.
    let (mut found, other) = if held_by_edge < held_by_contact {
        (from_edge, from_contact)
    } else {
        (from_contact, from_edge)
    };
    let ordered = |r: (f64, f64)| if r.0 <= r.1 { r } else { (r.1, r.0) };
    let unseen: Vec<ogeom_intersect::Overlap> = other
        .into_iter()
        .filter(|o| {
            let (lo, hi) = ordered(o.on_a);
            found.iter().all(|f| {
                let (flo, fhi) = ordered(f.on_a);
                hi < flo || lo > fhi
            })
        })
        .collect();
    found.extend(unseen);
    shared.insert((e.node, contact.node), swap(&found));
    shared.insert((contact.node, e.node), found.clone());
    Ok(found)
}

/// The stretch of `curve` over `crange` that the edge `e` runs along end to
/// end within `width`, as an overlap with the edge's own parameters at its
/// ends; `None` where some point of the edge stands farther off.
pub(crate) fn edge_along_curve(
    curve: &Curve,
    crange: (f64, f64),
    e: &BoundaryEdge,
    width: f64,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_intersect::Overlap>> {
    const SAMPLES: u32 = 48;
    let at =
        |i: u32| (e.crange.1 - e.crange.0).mul_add(f64::from(i) / f64::from(SAMPLES), e.crange.0);
    // Its ends and middle first: an edge that only crosses the curve
    // leaves it at one of them.
    let order = [0, SAMPLES, SAMPLES / 2]
        .into_iter()
        .chain((1..SAMPLES).filter(|&i| i != SAMPLES / 2));
    let mut feet = [0.0; 2];
    for i in order {
        let p = e.curve.point_at(at(i), tol)?;
        let (off, t) = nearest_on_edge_curve(curve, crange, p, tol)?;
        if off > width {
            return Ok(None);
        }
        if i == 0 {
            feet[0] = t;
        } else if i == SAMPLES {
            feet[1] = t;
        }
    }
    let (on_a, on_b) = if feet[0] <= feet[1] {
        ((feet[0], feet[1]), e.crange)
    } else {
        ((feet[1], feet[0]), (e.crange.1, e.crange.0))
    };
    Ok((on_a.1 - on_a.0 > tol.parametric()).then_some(ogeom_intersect::Overlap { on_a, on_b }))
}

/// The pcurve's course over a sub-range of the *curve's* parameters.
pub(crate) fn pcurve_polyline(
    pcurve: &PlanarCurve,
    prange: (f64, f64),
    crange: (f64, f64),
    sub: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Vec<Point2>> {
    let sagged = matches!(surface, SurfaceGeometry::Plane(_));
    chart_polyline(pcurve, prange, crange, sub, surface, sagged, tol)
}

/// [`pcurve_polyline`] with its chords held within [`strand_sag`] of the
/// pcurve on any surface, not only on a plane.
pub(crate) fn pcurve_polyline_within_sag(
    pcurve: &PlanarCurve,
    prange: (f64, f64),
    crange: (f64, f64),
    sub: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Vec<Point2>> {
    chart_polyline(pcurve, prange, crange, sub, surface, true, tol)
}

/// The pcurve's course over a sub-range of the curve's parameters, its
/// chords within [`strand_sag`] where `sagged`.
pub(crate) fn chart_polyline(
    pcurve: &PlanarCurve,
    prange: (f64, f64),
    crange: (f64, f64),
    sub: (f64, f64),
    surface: &SurfaceGeometry,
    sagged: bool,
    tol: Tolerances,
) -> OgeomResult<Vec<Point2>> {
    let lo = rescale(sub.0, crange, prange);
    let hi = rescale(sub.1, crange, prange);
    // Enough samples that the first step approximates the tangent and the
    // scanline interior test has a faithful outline. Straight pcurves get
    // two points. Everything else a fixed fine sampling.
    let mut count = match pcurve {
        PlanarCurve::Line(_) => 1,
        _ => {
            let span = (hi - lo).abs().max(SCAFFOLD_CHORD);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let n = (span / SCAFFOLD_CHORD).sqrt().ceil() as usize;
            n.clamp(8, 256)
        }
    };
    // And never a step [`unwrap_polyline`] could read as a period jump. A
    // straight pcurve is faithful at two points (a coaxial rim's image on
    // a cylinder is exactly a line across the chart), but one that crosses
    // more than half the turn looks to that reader like a wrap, and the
    // strand comes back running the complementary way: a bore's mouth
    // blended where the mouth circle starts anywhere but the wall's own
    // seam. The chart's span decides the count, not the parameter's.
    let ((ua, ub), (va, vb)) = surface.domain();
    let (first, last) = (pcurve.point_at(lo, tol)?, pcurve.point_at(hi, tol)?);
    for (periodic, period, reach) in [
        (surface.is_periodic_u(), ub - ua, (last.x - first.x).abs()),
        (surface.is_periodic_v(), vb - va, (last.y - first.y).abs()),
    ] {
        if !periodic || period <= 0.0 || !reach.is_finite() {
            continue;
        }
        let steps = (reach / (period / 3.0)).ceil();
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let steps = steps.clamp(1.0, 4096.0) as usize;
        count = count.max(steps);
    }
    if sagged && !matches!(pcurve, PlanarCurve::Line(_)) {
        let within = strand_sag(surface, tol);
        count = match circle_sag_count(pcurve, lo, hi, count, within) {
            Some(count) => count,
            None => sag_count(&|t| pcurve.point_at(t, tol), lo, hi, count, within)?.0,
        };
    }
    let mut out = Vec::with_capacity(count + 1);
    for i in 0..=count {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * i as f64 / count as f64;
        out.push(pcurve.point_at(t, tol)?);
    }
    Ok(out)
}
