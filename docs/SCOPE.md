# Scope

What a complete CAD kernel has to do, where each capability lands, and in what
order.

**Scope is set by the capability, not by any host application's demand for it.**
An earlier draft of this document tiered work by how heavily one large
application exercised each area. That was the wrong criterion and it produced
wrong answers: it marked hidden-line removal "deferrable" because only two of
that application's modules called it, and product structure "probably never"
because that application has its own document model. Both are things a CAD
kernel is expected to have. openGeometry is meant to be a complete, independent
kernel that any application can build on, so the target is functional
completeness.

Empirical usage data still has a use — it tells us which parts get hammered in
practice, and therefore what to make fast and get right first. That is a
*sequencing* input, not a *scope* input. See
[Appendix: usage evidence](#appendix-usage-evidence).

**Independence.** Nothing in this repository depends on, vendors, links against,
or commits any existing CAD kernel. The *Elsewhere* notes below say what a
capability is conventionally called, because that vocabulary is how the field
talks about itself and pretending otherwise helps nobody. They are a glossary,
not a dependency. See `CONTRIBUTING.md`.

**Interchangeability is a goal. A shared codebase is not, and cannot be.** The
target is that an application built on the conventional kernel can move to this
one with mechanical work rather than a redesign: the same concepts, the same
operation vocabulary, the same semantics, and an API of comparable shape. That
is what the *Elsewhere* notes are for, and it is why this document is organised
around the same capability boundaries.

What it licenses: matching documented and observable behaviour; keeping
operation names and argument meanings recognisable where doing so does not make
the Rust worse; implementing the conventional kernel's *published file formats*,
which is interoperation exactly as STEP is.

What it does not licence, and this is a constraint rather than a preference:
copying, translating or transliterating another kernel's source. The licences
run one way and it is not this way. Everything here is implemented from
published specifications, from the literature, and from tests — see
`CONTRIBUTING.md`, which says so normatively.

What it costs: nothing, in the two places `DATA_MODEL.md` deliberately diverges.
Stable provenance (§8) is *additive* — the conventional history API is
implemented and required alongside it, so an application that wants only history
gets it. The predicate trait (§9) is internal and has no API surface at all.

Where it does bite is **sequencing**: a capability the conventional kernel has,
and that applications therefore call, is worth more early than one it lacks.
That is a priority input, not a scope input — everything listed here is still in.

---

## Status key

| | |
|---|---|
| **P1** | foundation, geometry, topology |
| **P2** | construction, analysis, tessellation |
| **P3** | intersection, booleans, healing |
| **P4** | blending, offsetting, features |
| **P5** | data exchange, document model, drawings |
| **P6** | sketching, selection, recognition |

Phases order by dependency, not by importance. Everything listed is in scope.

---

## Milestones

Scale has defeated funded, full-time teams at this (see the prior-art note in
`README.md`). The answer is not optimism about the schedule — it is that **every
milestone below is independently useful**, so stopping at any one of them leaves
something worth having rather than a half-built kernel.

A milestone is **closed** when every criterion under it holds and
`./tools/check.sh` is green. The criteria are deliberately things that can be
*checked* rather than judged: a count, a comparison against a closed form, a
round trip that is or is not the identity. "It looks right" is not a criterion,
because in this domain it is worth nothing.

A milestone is not the whole of its phase. A phase is a body of capability; a
milestone is the point at which enough of it works to be worth shipping. What is
left over is listed under the milestone, and tracked in
[Deferred implementation details](#deferred-implementation-details).

### M1 — a B-rep data library · P1 · §1–4 · **closed**

- Every invariant in `DATA_MODEL.md` implemented, with property tests over the
  *laws* rather than over examples: orientation composition, location chain
  composition and inversion, the identity trichotomy, tolerance containment.
- Curves and surfaces reachable through the adaptor traits, so no algorithm
  names a concrete geometry type.
- A document round-trips through the native format: written, read, and written
  again gives the same bytes, with arena handles and entity identities
  preserved rather than renumbered.

### M2 — a modeling and tessellation library · P2 · §5, §6, §12 · **closed**

- Every primitive — box, cylinder, cone, sphere, torus, wedge — and both sweeps
  build a solid whose shell is closed and whose *mesh agrees with its topology*.
  That last is a separate check from validity and catches what validity cannot.
- Mass properties converge on the closed form, from the inscribed side, within
  the deflection each result reports.
- The same solid built two ways agrees: a revolved rectangle against
  `make_cylinder`, a revolved disc against `make_torus`, face for face.
- Point classification against a face and against a solid, with the band it
  cannot decide reported as `On` rather than guessed.
- Build, measure, check, tessellate and export from the command line.

Still owed inside P2, none of it blocking the above: whole-*shape* NURBS
conversion and whole-shape affine transforms (the geometry-level conversions
they rest on are done), and mesh → B-rep for *curved* regions, which is
canonical recognition and lands in §9.

### M3 — a CAD kernel · P3 · §7–9 · **closed, residuals listed**

- **The gate comes first.** ✓ Closed: the intersector is measured against
  analytic ground truth (worst 9.2e-15 over twelve pairs), a completeness
  instrument with a working negative control, literature-derived hard cases,
  and the corpus of real exchange files the reader now imports.
- General fuse with fuse, common, cut and section as selection predicates
  over one result. ✓ Closed for planar, quadric, marched and same-domain
  configurations, volumes held to closed forms and to each other; a cut runs
  on an imported NIST part end to end. *Residual:* the cut-imported result's
  mesh does not yet weld everywhere, so its volume is structural rather than
  measured — see the deferred table.
- Healing survives the corpus. ✓ Closed as stated: all eleven NIST parts
  read, every shell closes, ring re-anchoring makes the smallest part
  measure a real volume. *Residual:* parts with edges on B-spline surfaces
  have no pcurves yet and stay unmeasurable — closed shells, honest refusals
  — see the deferred table.
- History carried through every boolean, operation by operation. ✓ Closed:
  every source face is recorded modified-into-pieces or deleted, asserted in
  the boolean tests and on the imported cut.

### M4 — a manufacturing-capable kernel · P4 · §10, §11 · **closed, residuals listed**

- Constant- and variable-radius fillets, chamfers, shelling, offsets, pipe
  sweeps and lofts, on a mechanical part end to end. ✓ Closed: one part
  through the whole vocabulary — a block shelled into an open tray, corner
  edges chamfered, filleted at constant radius and at a running radius, a
  lofted boss and a swept nozzle fused onto its floor — valid, closed, and
  its volume the exact sum of every operation's closed form. *Residuals:*
  the blends speak straight edges and circular rims between analytic faces
  (marching blends, vertex and face-face blends deferred); the solid offset
  speaks planes and cylindrical bands (general rebuild deferred); sweeps
  and lofts stop where sweep surfaces would begin — every entry in the
  deferred table.
- Parity with the field is the bar here, not perfection: these are the two most
  fragile areas of every kernel that has them. The stones paid for themselves
  in kernel currency along the way: three torus arms and coincident circles in
  the analytic intersector, curved partner charts in the same-domain
  resolution, seams said honestly in the sweeps, and "touching is not
  crossing" enforced at every level from surface pairs down to paves.

### M5 — a documents-and-drawings kernel · P5 · §13, §15, §17

- A multi-body STEP assembly imported with colours, names and PMI, modified,
  exported, reimported, and compared against what went in.
- A 2D drawing generated from a 3D model: hidden lines removed, visible and
  hidden edges classified, sections taken.

### M6 — a complete kernel · P6 · §14, §16, §18

- A 2D constraint solver that reports degrees of freedom and *names the
  conflicting constraints* when a sketch is over-constrained, rather than
  merely failing.
- Ray and rectangle picking with sub-shape granularity, and a stable mapping
  from a triangle back to the topology that produced it.
- Features recognized from raw topology.

---

## 1. Foundation · `og-core` · P1 · *mostly done*

Arenas and entity identity, errors as values, the tolerance model, geometric
predicates. Still to come: units, progress reporting and cancellation,
deterministic parallelism.

*Elsewhere:* `Standard`, `NCollection`, `TCollection`, `Precision`, `Message`,
`OSD`, `Quantity`, `Units`.

## 2. Mathematics · `og-math` · P1

| | *Elsewhere* |
|---|---|
| Points, vectors, directions, axes, planes, quadrics | `gp` |
| Rigid transforms with form classification; general affine transforms | `gp_Trsf`, `gp_GTrsf` |
| Quaternions, Euler sequences | `gp_Quaternion` |
| Dense and sparse linear algebra: LU, QR, SVD, eigen | `math` |
| Root finding: Newton, Brent, bisection, polynomial roots | `math_FunctionRoots` |
| Optimization: BFGS, Powell, Newton minimization | `math_BFGS`, `math_Powell` |
| Numerical integration, Gauss quadrature | `math_GaussSingleIntegration` |
| B-spline and Bézier basis: de Boor, knot insertion and removal, degree elevation and reduction, subdivision | `BSplCLib`, `BSplSLib`, `PLib` |
| Elementary curve and surface parameterization and inversion | `ElCLib`, `ElSLib` |
| Approximation and fitting: least squares, constrained, smoothing | `AppDef`, `AppParCurves`, `AdvApprox` |
| Interval arithmetic for filtered predicates | — |

## 3. Geometry · `og-geom` · P1

Curve and surface traits first — the adaptor abstraction of `DATA_MODEL.md` §10 —
then the full type vocabulary in 2D and 3D.

| | *Elsewhere* |
|---|---|
| Analytic curves: line, circle, ellipse, hyperbola, parabola | `Geom`, `Geom2d` |
| Free-form curves: Bézier, B-spline, **rational (NURBS)** | `Geom_BSplineCurve` |
| Derived curves: trimmed, offset, curve-on-surface | `Geom_OffsetCurve` |
| Analytic surfaces: plane, cylinder, cone, sphere, torus | `Geom_ElementarySurface` |
| Free-form surfaces: Bézier, B-spline, **rational (NURBS)** | `Geom_BSplineSurface` |
| Procedural surfaces: revolution, linear extrusion, ruled | `Geom_SurfaceOfRevolution` |
| Derived surfaces: trimmed, offset | `Geom_OffsetSurface` |
| Local properties: derivatives, tangent, normal, curvature, torsion, principal curvatures | `GeomLProp`, `LProp` |
| Conversion: analytic ⇄ NURBS, B-spline ⇄ Bézier, knot and degree manipulation | `GeomConvert` |
| Projection onto curves and surfaces | `ProjLib`, `GeomProjLib` |
| **2D geometric constructions** — circles tangent to three entities, lines at prescribed angles, the classical straightedge-and-compass repertoire | `Gcc`, `Geom2dGcc` |
| Bisectors of curves | `Bisector` |
| Fair / minimum-energy curves | `FairCurve` |
| Plate and filling surfaces from boundary and point constraints | `GeomPlate`, `Plate`, `FEmTool` |
| Hatching in parametric space | `Hatch`, `Geom2dHatch` |

Rational (weighted) NURBS are load-bearing, not an optional extra: exact circles,
cylinders, cones and spheres in free-form representation all require them, and
every exchange format assumes they exist.

## 4. Topology · `og-topo` · P1

The B-rep model of `DATA_MODEL.md`: shared topology nodes, location chains,
orientation composition, per-entity tolerances, multi-representation edges.

| | *Elsewhere* |
|---|---|
| Vertex, edge, wire, face, shell, solid, compsolid, compound | `TopoDS` |
| Traversal composing location and orientation correctly | `TopExp`, `TopoDS_Iterator` |
| Same / equal / partner identity with matching hashers | `TopTools` |
| Geometry attachment; the builder as sole mutation path | `BRep_Builder`, `BRep_Tool` |
| Ancestor and adjacency maps | `TopExp::MapShapesAndAncestors` |
| **Non-manifold topology** — edges bounding more than two faces, mixed-dimension compounds | — |
| Cached triangulation and polygon representations | `Poly` |

## 5. Construction · `og-algo` · P2

| | *Elsewhere* |
|---|---|
| Vertex, edge, wire, face, shell, solid builders | `BRepBuilderAPI` |
| Primitives: box, cylinder, cone, sphere, torus, wedge | `BRepPrimAPI` |
| Prism, revolution, generalized sweep | `BRepSweep` |
| Sewing free faces into shells | `BRepBuilderAPI_Sewing` |
| Wire ordering, closing and repair during construction | `ShapeAnalysis_WireOrder` |
| Face from wires with outer/inner classification | `BRepBuilderAPI_MakeFace` |
| Copy, rigid transform, general affine transform | `BRepBuilderAPI_Transform` |
| Whole-shape conversion to NURBS | `BRepBuilderAPI_NurbsConvert` |
| `same_parameter` repair | `BRepLib::SameParameter` |
| Text and font outlines to wires | — |

Every operation emits history from the commit that introduces it
(`DATA_MODEL.md` §7).

## 6. Analysis and measurement · `og-algo` · P2

| | *Elsewhere* | |
|---|---|---|
| Bounding volumes: axis-aligned, oriented, per-sub-shape | `Bnd`, `BRepBndLib` | done |
| Point projection onto curves and surfaces | `GeomAPI_ProjectPointOn*` | done |
| Mass properties: length, area, volume, centroid, inertia, principal axes | `GProp`, `BRepGProp` | done |
| Point-in-face classification | `BRepTopAdaptor_FClass2d` | done |
| Point-in-solid classification | `BRepClass3d` | done, both: tessellated with the deflection band stated, and exact by ray casting against the true surfaces |
| Minimum distance and proximity between shapes | `BRepExtrema` | done: every element pair swept, interior approaches from the geometry-level extrema, boundary candidates owned by the lower-dimensional pairs |
| Extrema between curves and surfaces | `Extrema` | done: stationary approaches by Newton on the squared distance, constant-distance loci reported as a family rather than a guessed point |
| Validity checking against the model invariants | `BRepCheck` | |
| Curvature analysis; inputs for zebra, draft and thickness analysis | `BRepLProp` | |
| Self-intersection detection | `BOPAlgo_CheckerSI` | after §8 |
| Arc length, arc-length parameterization, deflection-based sampling | `GCPnts` | |
| Interpolation and fitting: points → curve, points → surface | `GeomAPI_PointsToBSpline` | |

## 7. Intersection · `og-intersect` · P3

The hardest area in the project, and the one with no existing Rust
implementation.

| | *Elsewhere* | |
|---|---|---|
| Curve/curve, 2D and 3D | `IntCurve` | done — analytic lines/circles with overlap detection, general via seeded Newton; 3D crossings carry their gap |
| Curve/surface | `IntCurveSurface`, `IntAna` | done — line/plane and line/quadric in closed form, general via the well-posed 3×3 Newton; lying-in-plane reported as an overlap |
| **Surface/surface** — analytic special cases first, then general marching with an approximation stage | `IntPatch`, `GeomInt`, `IntWalk` | done — `intersect_surfaces` is the one call: analytic exact with same-parameter pcurves, marched-and-fitted otherwise, tolerance as a sum of stated parts |
| Polyhedral pre-filtering | `IntPolyh` | inside the seeding: both surfaces sampled to triangles, pairs tested, candidates Newton-corrected |
| Contour and silhouette computation | `Contap` | with §13, whose consumer it is |
| Degenerate handling: tangential contact, coincident surfaces, poles, seams | — | tangency and coincidence answered by the analytic layer and refused honestly by the marcher; seams unwrapped in pcurve fitting; branch-point stitching still owed |

**Gated.** Surface/surface quality is measured against analytic ground truth and
published benchmark datasets before the boolean pipeline is committed to. If it
does not clear the bar, the project ships as a geometry library rather than
sinking years into a boolean built on an intersector that cannot carry it.

## 8. Booleans · `og-bool` · P3

One algorithm — general fuse — plus selection predicates over its result.

The **pipeline is in, one path for any analytic face**: sections from §7's
`intersect_surfaces` — exact with same-parameter pcurves where the projection
has a closed form, marched and fitted where not — paves from exact curve/curve
intersection, every face split in its own parameter space by an arrangement of
tagged strands (polyline scaffolding naming exact sub-curves), pieces
classified by the exact ray classifier, rebuilt with pcurves both sides, sewn,
closure demanded, shells nested into solids and voids, history recorded.
Fuse, common, cut and section all select from the one general result. A box is
drilled by a cylinder — seam, chart folding, holes in planar faces and a
flipped curved wall, all of it — the planar cases run through the same general
path, and crossed cylinders run through *marched* sections: fitted curves with
pcurves both sides, held to each other by the volume identities
fuse = A + B − common and cut = A − common. Same-domain and tangential contact
are refused with instructions — see the deferred table.

| | *Elsewhere* |
|---|---|
| Interference data structure, pave blocks, common blocks | `BOPDS` |
| Pave filler: dimension-ordered intersection, section edges, pcurve generation | `BOPAlgo_PaveFiller` |
| Builder: 2D face splitting, same-domain unification, solid rebuilding, tolerance repair | `BOPAlgo_Builder` |
| Fuse, common, cut, section, split | `BRepAlgoAPI` |
| Volume maker from an unordered face set | `BOPAlgo_MakerVolume` |
| Cells builder for arbitrary set expressions | `BOPAlgo_CellsBuilder` |
| Defeaturing / feature removal | `BOPAlgo_RemoveFeatures` |
| Periodic models, connected assemblies | `BOPAlgo_MakePeriodic` |
| Fuzzy tolerance, glue mode, non-destructive mode, parallel mode | `BOPAlgo_Options` |

## 9. Healing and simplification · `og-heal` · P3

| | *Elsewhere* |
|---|---|
| Analysis: wire order and orientation, edge/face consistency, tolerance survey, free bounds | `ShapeAnalysis` |
| Fixing: wires, faces, shells, solids, small edges and faces, missing pcurves, seams | `ShapeFix` |
| Upgrade: same-domain unification, face and edge splitting, canonical conversion | `ShapeUpgrade` |
| Reshape and substitution framework | `ShapeBuild_ReShape` |
| **Canonical recognition** — detecting that a NURBS patch is really a cylinder | `ShapeAnalysis_CanonicalRecognition` |
| Simplification: face merging, edge merging, tolerance reduction | — |

Real files are broken. This is not a postscript to the boolean work; it is what
makes the boolean work usable on anything that came from outside.

## 10. Blending · `og-fillet` · P4

Constant- and variable-radius edge fillets, vertex blends, face-face blends,
chamfers (distance-distance and distance-angle), 2D fillets and chamfers, blend
surface analysis.

**Opened**: the symmetric straight-edge chamfer between planar faces is in —
a wedge of five explicit planar faces subtracted through the boolean, whose
same-domain resolution is what melts the wedge's coplanar legs into the
solid's own faces. Exact by construction, volume pinned against the closed
form, history reading as the cut it is.

**The rolling-ball fillet** followed on the same scaffolding: for a straight
convex edge between planes the ball's envelope is exactly a cylinder, so the
blend is the chamfer's wedge with the bevel exchanged for the tangent
cylinder — legs to the tangency lines, caps closed by arcs, every pcurve
closed-form. The wedge's tangencies are what made it expensive: the blend
arc is both a boundary edge and a section curve, and until the 3D
intersector learned that two coincident circles overlap (the circle
counterpart of collinear lines, closed-form since), the arrangement carried
the same arc twice at two samplings and made sliver pieces; and the exact
classifier's `On` band now carries the ring polylines' own sag, without
which a probe between a tangent chord and its arc reads `Out` at a
resolution that cannot distinguish them. Volume pinned against the closed
form; both faces refuse concave edges — see the deferred table.

**The toroidal rim fillet** extended the family to its first revolved seat:
around the rim where a cylindrical wall meets a perpendicular cap, the ball's
envelope is a torus, and the wedge is the straight case's prism revolved — a
wall band to the tangency parallel, a cap annulus to the tangency circle, and
the quarter-tube between them, the bands built by `make_revolution_band` with
their seams. It paid its way in kernel currency: the analytic intersector
gained the axis-normal plane/torus, coaxial cylinder/torus and coaxial
torus/torus arms — tangencies reported as the circles they are, the way a
tangent plane reports its line on a cylinder — and the boolean's same-domain
resolution learned to invert probes into *curved* partner charts through the
closed-form `elementary` inversions, where before it could only speak plane.
Volume pinned against Pappus; the hole's rim — the concave case — refused
and recorded.

**The variable-radius fillet** closed §10's milestone criteria: for a linear
law on a straight edge between planes, the rolling ball's envelope is an
*exact* rational B-spline surface — degree one along the edge, a rational
quadratic across it, the control net affine in the radius — so nothing is
fitted: every section of the surface is that section's exact blend arc, the
tangency lines degenerate to straight rails, and the wedge subtracts through
the boolean like its constant-radius siblings. The kernel currency this
stone paid: the plane pcurve arm learned B-splines (affine invariance:
project the control net, keep knots and weights — which also gives imported
spline edges on planar faces exact pcurves); the marcher learned that
touching is not crossing at its level too, dropping branches along which the
two surfaces share their normal instead of fitting the tangential valley's
noise into phantom edges; curve-level tangential paves are skipped by the
same principle; sections that hug a boundary edge of either face are
excluded by measurement, not just by support recognition; and the boolean's
junction vertices widen their tolerance to the crossing residual they
actually carry, which is what vertex tolerance is for.

**Concave and mirrored blends** closed the sign ledger: convexity is read
from the face's own trim (the leg construction chooses its side and so
cannot answer it), a concave edge's wedge fuses where the convex one cut —
its legs opposed to the faces they melt against, which is exactly what a
fuse cancels — and the four revolved seats (external rim, hole rim, boss
base, blind floor) are one construction over two signs: wall-outward and
wall-side. Getting there fixed three kernel defects that reversed and
overlapped boundaries had been hiding: a contact overlapping a target's
boundary edge now paves the *target* too, so faces across the overlap split
where their new neighbours do; a partner chart's trim test now sees a
seam's both columns; and a seam occurrence picks its chart side by
*continuity with the ring being walked* rather than by orientation flags —
which a reversed face flips while the columns stay put.

**The chamfer's other spellings, and the sketch plane**: the asymmetric
distance-distance chamfer names a face for its first distance; the
distance-angle form derives the second distance where the bevel, leaving the
named face at the given angle, meets the other — and refuses angles that
never do. Both end in the one wedge the symmetric form already built. In 2D,
`fillet_corner_2d` and `chamfer_corner_2d` replace a wire's corner between
straight edges with a tangent arc or a cut at set distances, trimming the
edges on their own curves and rebuilding the wire with history — corner
deleted, edges modified, connector generated. Corners with curved sides are
the 2D tangency problem proper — see the deferred table.

*Elsewhere:* `ChFi2d`, `ChFi3d`, `Blend`, `BRepBlend`, `BRepFilletAPI`,
`FilletSurf`.

## 11. Offsetting, sweeping, features · `og-offset` · P4

**Opened** with the 2D wire offset: each straight or circular edge offset on
its own support, corners deciding the rest — gaps closed by an arc about the
old corner or by extension to the sharp meeting, overlaps trimmed to the
pieces' intersection, all against the wire's own winding. Areas pin to
Minkowski's arithmetic (perimeter times offset plus `πw²` for rounded
corners), offsets compose, and the history reads edge-into-offset-edge,
corner-into-join. The refusals are named: free-form edges, open wires,
offsets that consume an edge or an arc's radius, and self-intersecting
results — the arrangement that resolves a collapsed offset into its valid
loops is in the deferred table.

**Offset solids and shells** followed: the topology-preserving offset moves
every face's surface along its own outward normal — planes translate,
cylinder radii grow or shrink — and rebuilds the topology one-for-one on the
moved supports. Vertices re-solve where their three planes now meet (normal
equations, regularized along the edge for two), edges re-derive with their
directions and parameterizations preserved so orientation flags carry over
unchanged, and band faces rebuild through `make_revolution_band` so seams
stay seams. Corners stay sharp — the parallel solid, not the Minkowski body.
Shelling is that offset pointed inward and the boolean pointed at the
result: the cavity keeps the *removed* faces exactly where they were, so it
reaches the boundary at the openings, and the cut's same-domain resolution
melts the flush faces away — which is what opens the shell. A box shells to
its open tray and a cylinder to its cup, volumes exact against the closed
forms; an inside-out rebuild is caught by its own measured volume. The
vocabulary is planes and full cylindrical bands, straight edges and
axis-normal rings — the general rebuild is in the deferred table.

**Pipes and lofts** close the sweeps whose surfaces were already in the
vocabulary. A circular profile along a straight spine is the cylinder
primitive placed; around a full circle it is the torus primitive; along an
arc it is a torus segment built from two half-tube patches whose tube
circles are framed so their own parameter *is* the tube angle — every pcurve
a straight line in the chart, and the outer equator an honest seam between
the halves, said through `attach_seam`. Ruled lofts between parallel
sections: coaxial circles delegate to the cone and cylinder primitives,
matched polygons build planar walls with shared rails — and a twisted loft
is refused as the skew it is. Volumes pin to Pappus and the frustum
formulae. Free-form spines, oblique and mixed sections, and smoothed
skinning through many sections are the sweep-surface machinery, deferred by
name.

| | *Elsewhere* |
|---|---|
| Offset shape with selectable join type | `BRepOffsetAPI_MakeOffsetShape` |
| Shelling and thickening, with face removal | `BRepOffsetAPI_MakeThickSolid` |
| 2D wire offset | `BRepOffsetAPI_MakeOffset` |
| Pipe and pipe-shell sweeps with law-driven scaling and orientation | `BRepOffsetAPI_MakePipeShell`, `GeomFill` |
| Lofting through sections, ruled and smoothed | `BRepOffsetAPI_ThruSections` |
| Draft angle application | `BRepOffsetAPI_DraftAngle` |
| Evolved shapes | `BRepOffsetAPI_MakeEvolved` |
| Surface filling from constraints | `BRepOffsetAPI_MakeFilling` |
| Normal projection onto shapes | `BRepOffsetAPI_NormalProjection` |
| Form features: prism, revol, rib, slot, pocket, glue | `BRepFeat`, `LocOpe` |
| Bi-tangent construction | `BiTgte` |

## 12. Tessellation and meshing · `og-mesh` · P2

| | *Elsewhere* |
|---|---|
| Edge discretization to linear and angular deflection | `BRepMesh` |
| Constrained Delaunay triangulation per face in parametric space, with refinement | `BRepMesh_Delaun` |
| Deflection control: absolute, relative, adaptive | `IMeshTools_Parameters` |
| Watertight output with vertices shared across face boundaries | — |
| Normals, UVs and edge polygons on the triangulation | `Poly_Triangulation` |
| Mesh simplification and decimation | — |
| **Mesh to B-rep**: plane, cylinder, cone and sphere fitting; region growing; surface reconstruction | — |

## 13. Drawing generation · `og-hlr` · P5

Hidden line removal, exact and polygonal; visible and hidden edge
classification; section and broken-section views; silhouette and isoparametric
curve extraction; reflect lines.

*Elsewhere:* `HLRAlgo`, `HLRBRep`, `HLRTopoBRep`, `HLRAppli`.

Producing a 2D drawing from a 3D model is a core CAD capability. A kernel that
omits it is not a complete alternative, regardless of how few applications reach
for it directly.

## 14. Sketching and constraints · `og-sketch` · P6

A 2D geometric constraint solver: coincidence, distance, angle, parallel,
perpendicular, tangent, symmetry, equality, radius, horizontal and vertical;
construction geometry; driving and driven dimensions; degree-of-freedom
analysis; under- and over-constraint diagnosis that *names the conflicting
constraints* rather than merely failing.

*Elsewhere:* **nothing.** No mainstream B-rep kernel ships a constraint solver,
so every application built on one supplies its own. A kernel that intends to be a
complete foundation should not push this back onto every consumer.

## 15. Document, assembly and product structure · `og-doc` · P5

| | *Elsewhere* |
|---|---|
| Product structure: parts, assemblies, instances with placements | `XCAFDoc_ShapeTool` |
| Appearance: colours per shape and per sub-shape, materials, textures, layers | `XCAFDoc_ColorTool` |
| Names and user-defined properties | `TDataStd_Name` |
| **PMI and GD&T**: dimensions, tolerances, datums, annotations | `XCAFDoc_DimTolTool` |
| Validation properties (mass, area, centroid check values) | `XCAFDoc_ValidationProps` |
| Undo and redo over a transactional model | `TDF_Delta` |
| Native persistence: versioned, forward-compatible, lossless | — |

Required as a *capability*. Deliberately not a reimplementation of OCAF's
label-and-attribute tree or its file formats — that design is widely considered
more machinery than the job needs, and format compatibility is a non-goal.

## 16. Selection and display support · `og-select` · P6

BVH construction over tessellated and analytic geometry; ray picking with
sub-shape granularity; rectangle and polygon selection; depth sorting; level of
detail; a stable mapping from triangles back to the topology that produced them.

*Elsewhere:* `SelectMgr`, `Select3D`, `BVH`.

Not a renderer. `tools/ogview` consumes this; it is not part of it.

## 17. Data exchange · `og-io` · P5

| Format | Direction | Notes |
|---|---|---|
| Native `.og` | r/w | **Done.** Text, versioned, lossless: topology, geometry, placements, per-entity tolerances, provenance and the cached tessellation, with arena handles preserved rather than renumbered. Writing what was read gives the same bytes, so `diff` is the comparison tool. Anything it cannot write it refuses |
| Conventional kernel's native B-rep text format | r/w | Its own serialization, implemented from the published format description. Interoperation, on the same footing as STEP — it is what lets a corpus written by an existing application be read here, and it is the cheapest possible bridge for an application mid-migration |
| **STEP** AP203 / AP214 / AP242 | r/w | Assemblies, colours, names, AP242 PMI. The most valuable target — no production-grade Rust implementation exists |
| IGES | r/w | Legacy, still shipped by industry |
| glTF 2.0 | r/w | Tessellated, with materials |
| STL | r/w | ASCII and binary |
| OBJ, PLY, 3MF | r/w | |
| VRML | w | |
| DXF | r/w | 2D interchange |
| Parasolid X_T, ACIS SAT, JT | r | Later, where the formats are documented |

Parsing is the easy fifth. The rest is semantic mapping onto our topology, unit
and assembly-transform handling, and surviving files that violate their own
specifications — which most real files do.

## 18. Feature recognition · `og-recognize` · later

Recognizing holes, pockets, slots, bosses, fillets and chamfers from raw
topology; deriving feature trees from imported dumb solids; manufacturing
feature extraction.

Orthogonal to everything above and sequenced after the modeling stack is
trustworthy.

---

## Sequencing discovered during implementation

Ordering constraints that were not obvious from the capability list, recorded as
they were found. None of these is a reduction in scope — everything below is
still in — but each moved because building it earlier would have meant shipping
something that answers wrongly for inputs it appears to accept.

**Mass properties depend on tessellation (§12), not just on geometry.** Area and
volume are integrals over the *trimmed* region of a face, and the trimming is
what makes them hard: there is no closed form for the area of an arbitrary
region of a NURBS patch. The two honest routes are tessellation or 2D quadrature
with point-in-face classification, and tessellation is what production kernels
use. Writing an exact-for-planar-faces version first was the tempting
alternative and is a trap: it would give correct answers for a box and silently
wrong ones for a cylinder, with nothing in the signature to say which.

**Point classification lands with tessellation (§12); the *exact* version waits
for intersection (§7).** This was recorded the other way round and the ordering
turned out to be wrong. Ray casting against the tessellation needs no
ray/surface hits at all, and it is what makes classification available to
everything downstream years before the intersector exists.

What it costs is exactness, and the answer says so: a point nearer the boundary
than the deflection is reported `On` rather than assigned a side. That band is a
real limitation, not a rounding detail — but it is a *stated* one, which was the
original objection. The trap the earlier note was guarding against is answering
`In` or `Out` for a point the method cannot actually place, and returning `On`
is exactly not doing that.

The exact classifier still wants ray/surface intersection, and replaces this one
when the intersector lands. Point-in-*face* was already unblocked either way: it
is a winding count over the face's pcurves in parameter space.

**Distance and extrema between shapes depend on intersection (§7)** for the same
reason: the minimum distance between two curved faces is a constrained
minimization whose machinery is the intersector's.

**Self-intersection detection depends on booleans (§8).** It is the boolean
engine's interference stage run against a single argument.

---

## Verification

**There is no external oracle, and there will not be one.** Comparing results
against another kernel would mean vendoring one, which `CONTRIBUTING.md`
forbids — and it would make that kernel's bugs the definition of correct.

Geometry code fails quietly, so correctness is established from five independent
directions. The point of having five is that a defect which hides from one is
unlikely to hide from all; each of the first three has already caught a real one.

1. **Analytic ground truth.** Closed-form volume, area, centroid and length for
   every shape that has them, compared from the side the approximation must fall
   on — an inscribed mesh cannot exceed the surface it inscribes, and a test
   that only checks "close" would not notice a mesh that is wrong in the
   direction it cannot be.

2. **Metamorphic agreement.** The same solid built two ways must answer the same
   way. A rectangle with one side on the axis, revolved a full turn, is a
   cylinder — so it must have the face count and the volume `make_cylinder`
   gives, at the same deflection. This is the strongest substitute for an oracle
   available without one, because the two constructions share no code path.
   *It found the revolution's face orientation inverted.*

3. **Self-consistency.** A shape carries two descriptions of itself — its
   topology and its tessellation — and they can disagree. `check_tessellation`
   asks whether they do. The failure it names is invisible to every other check:
   face counts look right, the shell closes, each face triangulates without
   error, and the solid still has a slit down it. *It found the prism defect on
   its first run, and a second one the same day.*

4. **Round-trip identity.** Write, read, write: the same bytes. Not "close" —
   the same. A format that drifted a little on every save would make a real
   disagreement indistinguishable from noise. *Establishing this exposed two
   places where exactness was being thrown away in `og-math`.*

5. **Property tests over the laws.** Composition, inversion, antisymmetry,
   containment — stated as laws and tested over generated inputs, repeated
   across seeds by `tools/check.sh` so a case that only some seeds reach does
   not slip through.

**The corpus.** Generated permutations of primitives and operations today.
Published benchmark datasets for the intersection gate (§7), where they exist —
that is the one place a shared, external standard is available. Real exchange
files with §17, which are input to the kernel rather than another kernel's
output about a shape, and so carry no dependency.

**The harness.** `./tools/check.sh` is it, for now: format, lints, docs, and the
suite repeated across seeds, exit-code driven. A separate corpus runner arrives
with §17, when there are outside files to run and a per-operation pass rate is a
number that means something. Building one earlier would be a dashboard reporting
on inputs we generated ourselves.

---

## Deferred implementation details

Places where something *is* built and works, but a narrower or exact version is
still owed. Each is recorded with the section that will pay it off, so it is a
scheduled debt rather than a surprise found by whoever hits it.

The rule these follow: **an approximation that says so is allowed; one that does
not is not.** Every entry below either reports its own accuracy or refuses the
input it cannot handle. Nothing here answers confidently and wrongly.

| Owed | Where it lands | What exists now |
|---|---|---|
| Analytic mass properties for analytic surfaces | §6 | Tessellation-based, with the chord deflection reported on every result (`MassProperties::deflection`). Correct within a stated band, and exact already for planar faces and straight edges, since their tessellation is exact. A closed form for a trimmed cylinder or sphere would be faster and exact; it is an optimization, not a correction. |
| `same_parameter` repair | §5 | The flag is carried and set false whenever a representation is added, which is honest but pessimistic: every primitive's edges claim disagreement they do not have. The repair routine has to verify agreement, and where it fails, widen the edge's tolerance until the claim is true. Until it exists, nothing may *rely* on the flag being true. |
| `PolygonOnTriangulation` edge representation (`DATA_MODEL.md` §6) | §12 | A face carries a `Triangulation` and an edge carries a `Polyline` with its parameters, which is what makes the stored tessellation watertight. The missing piece is the edge's path through a *specific* face's triangulation as node indices — what a renderer wants to draw a shared edge without hunting for coincident vertices. |
| Predicates routed through the `Predicates` trait, everywhere | §7 | The seam is real where it decides something *combinatorial*: `inside_boundary` — which decides whether a triangulated face keeps a triangle or drops it — goes through `orient2d`, and `inside_boundary_with` names the implementation. Pinned by a test that finds cases where the two implementations genuinely disagree, and a second showing the exact one is right down to one unit in the last place. What is *not* routed: the ray/triangle test in classification, and the tolerance comparisons in mass properties. Those are measurements against a tolerance rather than sign decisions, and forcing them through a predicate would be ceremony. The marching intersector is written against the trait from the start, which is why this moved here from §9.
| Collapsed wedges | §5 | `make_wedge` refuses a top extent of zero. A wedge tapering to a ridge has five faces and one tapering to a point has four; both are different topologies, not this one with a zero in it, and building them through the box path would give them a face with no area. They belong to the prism sweep over a triangular profile. |
| Compsolid builder | §5 | `make_compound` now exists alongside the other builders, history included. A *compsolid* — solids gluing along shared faces — remains: nothing builds one and no consumer asks yet. |
| Non-manifold topology (§4) | §4 | The model permits an edge bounding more than two faces — nothing rejects it — but nothing has been built that way and no test pins the behaviour. `is_shell_closed` counts edge *uses* and asks for an even number, which is the right rule for non-manifold input, but that is a design choice not yet exercised. |
| Deflection in parameter units | §12 | `discretize_planar` measures its chord tolerance in parameter units, not in space, because it has no surface to convert through. Callers wanting a spatial tolerance convert themselves. The face boundary already avoids this by discretizing the edge's 3D curve and sampling the pcurve at those parameters — the fallback path only runs for an edge with no 3D curve at all. |
| ~~Smoothness across a fitted closed loop's seam~~ **done** | §6 | `fit_points_closed` makes a loop's join C1, the constraint eliminated exactly inside the least-squares solve; opt-in, because it spends shape freedom `fit_points` still promises. |
| ~~Points to a B-spline *surface*~~ **done** | §6 | `fit_surface_grid`: rows then columns on shared knot vectors at fixed parameters, the error measured surface-against-every-input. Scattered data still needs a grid first. |
| General affine transforms over a whole *shape* | §5 | Done at the geometry level: `Curve::general_transformed` converts to the exact B-spline form and moves the control points, which an affine map does exactly, so a sheared circle lands on the ellipse it should rather than near it. `transformed` still takes only a `Transform`, so a placement keeps the analytic type and only this does not. What is owed is the shape-level operation, and it is owed for the same reason whole-shape NURBS conversion is: the parameterization moves, so every edge's range and every pcurve has to be re-derived. |
| Whole-shape NURBS conversion | §5 | Exact conversions exist for every curve — line, circle, ellipse, hyperbola, parabola, spline, trimmed — and for the ruled surfaces: plane, cylinder, cone, extrusion. Each is exact rather than fitted, and the tests measure it *implicitly*, by how far the patch strays from the analytic surface, because comparing parameters would report the reparameterization instead of the error. **What is owed is the whole-shape operation.** Conversion necessarily reparameterizes — a circle's parameter is its angle and a rational quadratic's is not, and no reparameterization of a rational quadratic makes it one — so converting a shape means restating every edge's range and re-deriving every pcurve against a surface whose parameterization has also moved. Re-deriving a pcurve is a fit, which the adaptive fitting core now provides; what is missing is the operator that walks the shape and restates every range against the moved parameterizations. |
| ~~Exact patches for revolved surfaces~~ **done** | §3 | `revolved_patch` revolves the exact rational profile by the exact rational unit circle, weights multiplied; sphere, torus and a revolved line verify to 1e-12 against independent authorities. Trimmed still refuses. |
| ~~Sweep surfaces proper~~ **done** | §11 | `make_loft_skinned` and `make_pipe_skinned` build skins over the grid fit — seams closed exactly through pinned row ends, pcurves as iso lines — with volumes checked against the frustum and Pappus. Closed free-form spines remain unspoken. |
| ~~The general offset rebuild~~ **done** | §11 | Vertices solve their faces' planar displacement constraints and Newton-polish onto the moved surfaces; edges re-derive from the moved pairs, circles re-framed on old axes so parameters carry; all five analytic surfaces offset. Free-form faces remain refused. |
| ~~Collapsed 2D offsets~~ **done** | §11 | The raw offset splits at its self-crossings, sub-pieces stand or fall by distance to the source, survivors reconnect into loops; open wires offset into capped thick-path outlines. |
| ~~2D corner blends with curved sides~~ **done** | §10 | Tangent-circle centres found on the sides' offset loci, candidates qualified GccAna-style; chamfer distances on an arc are arc lengths. Free-form sides remain refused. |
| ~~Concave blends~~ **done** | §10 | Blends dispatch on convexity read from the face's own trim; the additive wedge fuses where the subtractive one cut; the four revolved rim seats are one parameterization over two signs. |
| ~~Canonical surfaces from a sweep~~ **done** | §9 | `make_prism` builds straight walls as planes on the extrusion's own chart, so the boolean's plane-to-plane resolution sees them; a chamfer melts into a prism wall. |
| ~~Canonical surfaces from a revolution~~ **done** | §9 | `make_revolution` builds a radial line's sweep on the plane it turns in — no seam, no degenerate centre; a revolved rectangle matches `make_cylinder` count for count. |
| Asymmetric conic domains in the native format | §17 | A hyperbola or parabola is written as the symmetric `[-extent, extent]` its constructor produces, and the reader refuses anything else rather than silently re-centring it. Unreachable today — nothing builds an asymmetric one — and it is a refusal rather than a wrong answer, but the format has no way to express one if something later does. |
| A consumer for the half-space | §8 | `make_half_space` builds the solid on one side of a face, orienting the boundary away from the material. Nothing reads it yet, and nothing can: its shell is one face with free edges all round, so it is not closed, and every query needing an inside refuses rather than guessing — which is the correct answer for a boundary that does not close. A half space is an argument to a boolean, and the boolean is §8. Its bound comes back *empty* for the same honest reason: `surface_bounds` declines to bound an unbounded plane rather than quoting the declared extent. |
| Curved regions in mesh → B-rep | §9 | Planar recovery is done and exact: coplanar triangles are grown into regions, their boundaries chained into loops, and faces rebuilt on shared edges so the shell closes — a box round-trips through a mesh and measures the same. A *curved* region is refused. Fitting a cylinder to a band of triangles is easy; deciding that the band **is** a cylinder rather than a smooth patch resembling one is the difficulty, and a wrong answer gives a solid that looks right, measures nearly right, and has the wrong surface under every later operation. The crease angle that separates a model edge from a tessellation seam is exposed as a parameter with a stated default, because no fixed angle is right for every mesh — a coarse tessellation of a large cylinder turns further than a fine one of a small cylinder. Deciding it from the mesh's own deflection is canonical recognition. |
| Branch points in the marching intersector | §7 | Where the two normals are parallel the intersection has no single direction, and `trace` stops and says `Stalled` rather than marching through — pushing on is how a tracer changes branch without noticing, and a wrong branch is a plausible answer to a different question. What is owed is *stitching*: detecting the branch point, and joining the fragments that meet there into the curves they belong to. Until then a configuration with a genuine singularity comes back as fragments that each say why they stopped. |
| Boolean coverage beyond the analytic pairs | §8 | The pipeline runs on any faces §7's intersectors answer for, planar and quadric alike, through one general path — the drilled box is the pinned proof. What remains named: same-domain contact (coplanar or co-surface faces) refuses rather than unifies; tangential contact refuses where the touch lies within both faces' reach; scaled placements refuse because a scale changes a surface's parameterization out from under its pcurves; and a marched pair that resolves no branch refuses rather than guesses. The refusal gates are generous by design — a face's reach is its boundary's bound plus a bulge allowance, and a curved face's allowance is most of its own diagonal, so an innocent configuration can be refused but a degenerate one cannot slip through. |
| Same-domain faces in the boolean | §8 | Handled for planar pairs: two faces on one plane split each other by the *other* face's boundary edges — coincident surfaces offer no section curve, so the contact's edges are the splitting curves, projected exactly into the shared chart — and each On piece resolves by comparing outward normals with its coincident partner. Opposed materials vanish from a fuse and survive a cut on the first argument's side; aligned materials keep one copy for fuse and common and none for cut. Stacked boxes, partial stacks and flush walls are pinned. Still refused by name: curved same-domain pairs whose edges have no closed-form projection into the shared chart, and contact reaching only an edge or a vertex, where there is no coincident partner *face* to compare sides against. |
| Tangential contact along a curve | §7 | A plane resting on top of a torus touches it along a whole circle, and there is no transversal curve there to march: the crossing angle is zero along the entire contact. The analytic layer names this contact where it has the closed form (a plane tangent along a cylinder's ruling comes back `Along`), but for the pairs the marcher owns, seeds near the contact converge — the two surfaces sit within the correction's acceptance of each other — wander briefly, and stall. What comes back is fragments hugging the contact curve: on both surfaces to rounding, describing nothing, each saying `Stalled`. The pinned test `tangency_along_a_circle_produces_fragments_not_a_curve` holds that behaviour still. The honest answer needs tangential contact traced as its own kind of curve — following the valley of the gap function rather than a crossing — which is a different walker than the one built. |
| Crossings shallower than a microradian | §7 | The marcher refuses a crossing whose angle is below `SHALLOWEST` (sine 1e-6), and the number is set by the correction rather than by taste: the Newton step may leave a residual up to the confusion tolerance, and on near-coincident surfaces that residual masquerades as an angle between the two computed normals. A gate below that floor reads the correction's own noise as a direction — identical spheres came back as six confident curves that existed nowhere but in rounding. A genuine sub-microradian crossing is also one the correction could not have followed, so the refusal costs nothing that was available. Following such crossings needs the higher-precision arithmetic the `Predicates` seam exists to admit. |
| Smoothness of a fitted intersection across its own seam | §7 | The approximation stage is done: a traced branch becomes a 3D spline and one pcurve per surface, each fitted to a stated tolerance, with the error reached reported and the total against the true intersection stated as trace chord plus fit error. Pcurves crossing a periodic surface's seam are unwrapped before fitting, so they run continuously and may leave the stated domain — which is what crossing a seam is. The fitting core now offers `fit_points_closed` for a C1 join; what remains is the approximation stage adopting it for closed branches — the joint seven-dimensional fit needs the same constraint before an intersection loop stops creasing. |
| Exact pcurves for every analytic section | §7 | `intersect_surfaces` derives exact same-parameter pcurves where the projection has a closed form: any conic on a plane, axis lines and full circles on a cylinder, parallels of latitude on a sphere. Where it does not — an oblique ellipse on a cylinder runs as a cosine in parameter space — the pcurve is `None` rather than a fit, because an exact curve with a fitted pcurve would be a curve whose descriptions disagree by an amount nothing on it records. The consumer that needs those is the boolean, which marches the pair instead. |
| Restriction of closed exact sections to surface extents | §8 | Section *lines* are clipped to both surfaces' parameter extents through their exact pcurves, and closed curves wholly outside an extent are dropped — but a closed curve partially outside is kept whole. Cutting it into arcs is the restriction problem, and the restriction that matters is the face's trim, which is the boolean's 2D splitting stage; the surface extent is only the parameterization window. |
| ~~Pcurves for edges on B-spline surfaces~~ **done** | §9/§17 | Projection fitting at the curve's own parameters serves every corpus face; slop up to `confusion*1e6` is accepted and *recorded* on the edge's tolerance. Parts still refusing do so for the face-level reasons below. |
| ~~Welding the cut-imported result~~ **done** | §8/§9 | Edge ends anchor to their vertices within the recorded tolerance, chart points fold onto the ring-continuing periodic branch, and a border-only weld-and-stitch pass runs at the recorded tolerance, floored at a tenth of the chord. The imported cut measures. |
| Corpus faces that still refuse to triangulate | §9 | Nine corpus parts stay unmeasured for reasons the scan bench isolated per part. Three classes remain. (1) A face bounded by *one* ring that winds a periodic direction — a torus band, a cone around its apex — has a uv boundary that is an open chain once folded, so the trim culls incoherently; it needs the seam-and-split synthesis the two-ring bands already get (ftc_06/07/08/10 cones and tori, ctc_01's apex cones). (2) A plane face whose boundary insertion reports `TooSmall` (ctc_03). (3) ctc_02/04/05 close to within a few hundred border segments of measurable, all of class 1. The reader-side pcurve and weld machinery is done; what remains is per-face boundary synthesis. |
| Spindle-torus signed distance | §2 | `Torus::distance_to` handles every torus, including the folded branch of a spindle. `signed_distance_to` is documented as ring and horn tori only: a spindle torus encloses two regions and "inside" does not name one of them. |

---

## Explicitly out of scope

| | Why |
|---|---|
| Rendering | A kernel supplies geometry and picking, not a render loop. `tools/ogview` is a *software* renderer to an image, and it is a verification tool rather than a feature: it exists so a wrong result is visible, it has no dependency beyond the kernel, and `check.sh` can therefore run it. A real-time windowed viewer would wrap it and add interaction, not correctness. |
| FEA meshing — tetrahedral, boundary-layer | A separate discipline with its own literature. The kernel supplies the geometry and surface meshes it consumes. |
| Simulation, CAM toolpath generation | Applications, not kernel. |
| A parametric *recompute engine* | Application-level. The kernel supplies stable identity and history (`DATA_MODEL.md` §7, §8), which is what such an engine needs from it. |
| Binary ABI compatibility with another kernel | Not expressible from Rust without a C++ façade carrying the other kernel's class layout, and any host API that hands a raw shape pointer to a third-party binding runtime needs that layout to match exactly. Firmly downstream — see `docs/INTEGRATION.md`. |

---

## Appendix: usage evidence

`tools/apisurf` measures how an application exercises a kernel of this shape.
Run against one large, mature CAD application it reported 574 headers, 591 types
and 4,126 include sites, with per-member call counts and module attribution.

**The generated profile is not committed.** It is derived from another project's
headers, it is regenerable in seconds by anyone with both trees checked out
locally, and nothing in this repository depends on it. Generate your own:

```sh
python3 tools/apisurf/apisurf.py --reference /path/to/reference/headers \
                                 --consumer  /path/to/application \
                                 --out docs/api_surface.json
```

Read it as a **profile**, not a specification.

Good for:

- **Sequencing.** Shape topology at 5,398 references and primitive geometry at
  4,444 sit on every hot path; they earn the most design attention and the
  earliest benchmarks.
- **Ergonomics.** The members called a thousand times are the ones whose API
  shape matters. The core shape type is touched through 22 of its 38 members.
- **Realism.** It shows which capabilities carry weight in a real application
  rather than only on a feature list.

Not good for deciding **what to build**. Reference counts systematically
understate anything reached through a narrow façade — booleans show 38
references and represent perhaps a fifth of the implementation effort — and they
say nothing whatsoever about capabilities the sampled application implements
itself, such as constraint solving, or reaches for rarely, such as drawing
generation.

---

## Appendix: dependencies considered

A kernel that pulls in the wrong dependency inherits its representation, and a
representation is the one thing that cannot be swapped out later. So the ones
turned down are recorded here with the reason, not only the ones taken —
otherwise the question gets re-litigated every time someone notices a crate that
looks like it would help.

Every entry states a decision and why. Nothing here asserts a license:
`deny.toml` enforces the licence rule at build time, and it is checked at
adoption rather than remembered from a note.

**Taken.** Rationale lives beside each in the workspace `Cargo.toml`.

| | For |
|---|---|
| `nalgebra` | The arithmetic substrate. Generic over `RealField`, so extended-precision or interval scalars can be swapped in later without rewriting algorithms |
| `robust` | Shewchuk adaptive-exact predicates: `orient2d`, `orient3d`, `incircle`, `insphere`. Fast floating-point filter escalating to exact only when the error bound leaves the sign undecided |
| `spade` | Constrained Delaunay with refinement, over the same predicates |
| `hashbrown`, `indexmap`, `smallvec`, `thiserror` | Collections and error derivation |

**Not taken.**

| | Why not |
|---|---|
| `glam` | `f32` and graphics-oriented. A kernel is `f64` end to end |
| `cgmath` | Unmaintained |
| `csgrs` | Mesh CSG. It destroys the analytic surfaces that are the entire point of a B-rep — the result of a boolean between two cylinders has to still *be* cylindrical, not a triangulation that resembles one |
| `earcutr` | Ear clipping: no quality guarantee, and no way to insert an interior point. A curved face needs both, which is why `spade` was taken instead |
| `rapier3d` | Physics. Nothing here needs it |
| `parry3d` | Proposed for BVH broad-phase only. A BVH over our own triangulations is a small thing to own, and taking a physics crate's spatial types to get one would put a second geometry representation in the tree. Revisit at §16 if picking turns out to want more than a BVH |
| `inari` | Interval arithmetic, for filtered predicates. Genuinely wanted eventually (§2), and deferred until the filtered predicates it would serve are actually built — the `Predicates` trait is the seam it plugs into, and that seam already exists |
| `curvo` | NURBS evaluation, curve/curve intersection, 2D region booleans and trimming. The 2D booleans are close to what §8's face-splitting stage needs, so this is the one on the list worth genuinely re-examining — **at the §7 gate**, not before. Taking it earlier would mean two NURBS representations in the tree while `og-geom` is still being shaped |
| `truck-geometry` / `truck-topology` | Evaluated as a *reference design* rather than a dependency. Taking a whole topology crate means taking its data model, and `DATA_MODEL.md` diverges from the conventional one deliberately in two places (§8, §9). Those divergences are the point |

**A correction.** An earlier note in this project's planning claimed that
`robust` was 2D-only and that a separate crate would be needed for `orient3d`
and `insphere`. That is not true of the version in use: it provides all four
predicates, and `og-core`'s exact implementation calls them directly. Recorded
so the claim is not re-derived from the old note.
