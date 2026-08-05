# Scope

What a complete CAD kernel has to do, where each capability lands, and in what
order.

**Scope is set by the capability, not by any host application's demand for it.**
An earlier draft of this document tiered work by how heavily one large
application exercised each area. That was the wrong criterion and it produced
wrong answers: it marked hidden-line removal "deferrable" because only two of
that application's modules called it, and product structure "probably never"
because that application has its own document model. Both are things a CAD
kernel is expected to have. ogeom is meant to be a complete, independent
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

### M5 — a documents-and-drawings kernel · P5 · §13, §15, §17 · **closed, residuals listed**

- A multi-body STEP assembly imported with colours, names and PMI, modified,
  exported, reimported, and compared against what went in. ✓ Closed: the
  bolted-plate assembly arrives as three named products with per-product and
  per-face colours and its semantic PMI (a toleranced dimension, a flatness,
  a datum); the plate is drilled through the boolean, a third bolt instanced,
  a colour and a dimension added; the written file reimports with the same
  products, four occurrences at their placements, the assembly volume equal
  to the arithmetic of the modification, and every annotation intact.
  *Residuals:* PMI beyond the semantic core and §15's transactional
  machinery, in the deferred table — modifiers and composite datum
  references have since been paid; presentation polylines, datum targets
  and undo/redo remain, and land with M11.
- A 2D drawing generated from a 3D model: hidden lines removed, visible and
  hidden edges classified, sections taken. ✓ Closed: a bored block in
  three-quarter view classifies every curve, draws the bore by silhouette
  and keeps the far side dashed rather than gone; the section through the
  bore reveals the rectangle minus the chord-width slot to arithmetic.
  *Residuals:* the HLR is polygonal (exact interference deferred), and a
  section exactly through a bore's axis meets the boolean's tangential
  refusal — both in the deferred table.

### M6 — a complete kernel · P6 · §14, §16, §18

- A 2D constraint solver that reports degrees of freedom and *names the
  conflicting constraints* when a sketch is over-constrained, rather than
  merely failing. ✓ Closed: the full constraint vocabulary over damped
  Gauss–Newton with SVD structure-reading — a dimensioned bracket solves and
  reports itself exactly constrained; adding a disagreeing width makes the
  solve refuse and name the two distances that fight; the unanchored
  rectangle names its sliding motions and every movable point; a repeated
  dimension is called redundant, not conflicting.
- Ray and rectangle picking with sub-shape granularity, and a stable mapping
  from a triangle back to the topology that produced it. ✓ Closed: a
  median-split BVH over per-face tessellation, hits depth-sorted, an
  aperture resolving vertex before edge before face, rectangle and even-odd
  polygon marquees in inside and crossing modes, and triangle indices
  assigned in traversal order — the ray cast into a drilled part strikes a
  triangle whose stable owner is a face the recognized hole claims.
- Features recognized from raw topology. ✓ Closed: holes (through, blind,
  counterbored across their shoulders, countersunk), fillets by two-sided
  tangency (partial cylinders and torus bands), chamfers, pockets and
  slots, tested on built parts and on the corpus — where the recognizer
  stays honest: ftc_11's rim tori are tangent on one side only, and it
  refuses to call them fillets. The residuals paid since: bosses, partial
  rounds as their own category, and `feature_tree()` — all against a parity
  oracle that asks the solid itself.

### The conclusion — M7 through M12

M6 closed the last of the phase milestones while the section tables still
hold unbuilt rows and the deferred table still holds milestone-scale debt.
These six milestones are the remainder, all of it: when M12 closes, every
row of §1–§18 is either done or refused by name with the refusal pinned by
a test and recorded as design rather than debt. Ordering is by dependency,
as ever.

### M7 — exact-geometry completeness · §1, §2, §3, §6, §7

- Predicates: the ray/triangle classification test routed through the
  trait; an interval type with outward rounding deciding "genuine tangency
  or correction noise" at the marcher's shallow-crossing gate; sparse
  linear algebra serving the sketch solver's normal equations.
- Progress reporting and cancellation through every long operation;
  deterministic parallel tessellation whose output is bit-identical at any
  thread count; a benchmark harness with a checked-in baseline.
- The §3 repertoire: offset curves and surfaces and curve-on-surface as
  first-class types; the Gcc tangency constructions with qualified
  solutions; bisectors; fair curves; plate/filling surfaces from boundary
  constraints; hatching; text to wires; a helical curve with exact
  derivatives and closed-form arc length, swept into a thread.
- §6 closed: draft and thickness analysis with their sampling stated;
  self-intersection detection as the interference stage against one
  argument.
- §7's deferred entries paid: sinusoidal exact pcurves where the closed
  form exists, branch-point stitching, tangential contact traced as its
  own curve, and the sub-microradian gate decided by intervals rather
  than taste.

### M8 — boolean and healing completeness · §8, §9

- Curved same-domain pairs unified through the elementary inversions;
  edge- and vertex-only contact handled by name; partially-outside closed
  sections restricted to arcs.
- The option set: fuzzy tolerance, glue mode, non-destructive pinned;
  a consumer for the half-space; volumes from an unordered face soup;
  the cells builder over the general-fuse result; scaled placements baked
  through the whole-shape conversion or refused with instructions.
- Canonical recognition — a NURBS patch confessed as the cylinder it is,
  within a stated tolerance — and, on it, curved regions in mesh→B-rep;
  same-domain face unification, edge-chain merging and tolerance
  reduction as standalone upgrades; the reshape substitution framework.

### M9 — blend and feature completeness · §10, §11

- The restored blend vocabulary: face-face blends, marching blends over
  offset-surface spines, the setback vertex blend, and blend analysis
  that reports tangency error instead of asserting it.
- Draft angle about a neutral plane; the form-feature vocabulary (prism,
  revol, rib, slot, pocket, glue) as history-carrying compositions;
  normal projection of wires onto shapes; evolved shapes over planar
  spines; bi-tangent construction where those need it.

### M10 — exact drawings · §13

- Exact silhouettes on the analytic surfaces, projected edges and
  silhouettes split at their exact 2D crossings, visibility decided by
  the exact classifier — the polygonal pipeline kept as the measured
  fallback. Isoparametric extraction and reflect lines. Broken sections;
  the on-axis half-section once M8's contact handling lands.

### M11 — interaction and document completeness · §14, §15, §16, §18

- Sketch: analytic Jacobians pinned against the numeric oracle,
  construction geometry, driven dimensions, snapping.
- Selection: one pick structure serving several deflections.
- Documents: undo/redo over transactions, textures, datum targets,
  presentation PMI polylines both directions.
- Recognition: the feature tree mapped to manufacturing operations.

### M12 — exchange completeness · §17

- The conventional kernel's B-rep text format, both directions, from the
  published description. IGES both directions, staged like STEP was.
  Readers for the mesh formats already written (glTF, OBJ, PLY), VRML
  and 3MF writing and reading, DXF reading. The documented proprietary
  tail — SAT, X_T, JT — read exactly as far as public documentation
  carries, refusals named past that line.

---

## 1. Foundation · `ogeom-core` · P1 · *mostly done*

Arenas and entity identity, errors as values, the tolerance model, geometric
predicates, units carried by the model and honoured by the formats. Progress
and cancellation: a scope-installed `Watch` whose `checkpoint()` every long
loop offers — tessellation, the boolean pipeline, the marcher, the STEP
reader — cancellation landing as an error, never a partial result.
Parallelism: `parallel::map_ordered`, whose rule is that the answer is
bit-identical at any thread count; face tessellation runs through it, its
compute phase immutable and its attachment sequential in face order, pinned
by a byte-equality test on the native serialization.

*Elsewhere:* `Standard`, `NCollection`, `TCollection`, `Precision`, `Message`,
`OSD`, `Quantity`, `Units`.

## 2. Mathematics · `ogeom-math` · P1

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

## 3. Geometry · `ogeom-geom` · P1

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
| Helical curves, with exact derivatives and closed-form arc length | — (the conventional kernel approximates helices; threads deserve better) |
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

## 4. Topology · `ogeom-topo` · P1

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

## 5. Construction · `ogeom-algo` · P2

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

## 6. Analysis and measurement · `ogeom-algo` · P2

| | *Elsewhere* | |
|---|---|---|
| Bounding volumes: axis-aligned, oriented, per-sub-shape | `Bnd`, `BRepBndLib` | done |
| Point projection onto curves and surfaces | `GeomAPI_ProjectPointOn*` | done |
| Mass properties: length, area, volume, centroid, inertia, principal axes | `GProp`, `BRepGProp` | done |
| Point-in-face classification | `BRepTopAdaptor_FClass2d` | done |
| Point-in-solid classification | `BRepClass3d` | done, both: tessellated with the deflection band stated, and exact by ray casting against the true surfaces |
| Minimum distance and proximity between shapes | `BRepExtrema` | done: every element pair swept, interior approaches from the geometry-level extrema, boundary candidates owned by the lower-dimensional pairs |
| Extrema between curves and surfaces | `Extrema` | done: stationary approaches by Newton on the squared distance, constant-distance loci reported as a family rather than a guessed point |
| Validity checking against the model invariants | `BRepCheck` | done: `check` grades every finding by severity, `check_tessellation` holds the mesh against the topology |
| Curvature analysis; inputs for zebra, draft and thickness analysis | `BRepLProp` | local properties done at the geometry level; draft and thickness analysis land with M7 |
| Self-intersection detection | `BOPAlgo_CheckerSI` | lands with M7: the interference stage against one argument |
| Arc length, arc-length parameterization, deflection-based sampling | `GCPnts` | |
| Interpolation and fitting: points → curve, points → surface | `GeomAPI_PointsToBSpline` | |

## 7. Intersection · `ogeom-intersect` · P3

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

## 8. Booleans · `ogeom-bool` · P3

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

## 9. Healing and simplification · `ogeom-heal` · P3

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

## 10. Blending · `ogeom-fillet` · P4

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

## 11. Offsetting, sweeping, features · `ogeom-offset` · P4

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

## 12. Tessellation and meshing · `ogeom-mesh` · P2

| | *Elsewhere* |
|---|---|
| Edge discretization to linear and angular deflection | `BRepMesh` |
| Constrained Delaunay triangulation per face in parametric space, with refinement | `BRepMesh_Delaun` |
| Deflection control: absolute, relative, adaptive | `IMeshTools_Parameters` |
| Watertight output with vertices shared across face boundaries | — |
| Normals, UVs and edge polygons on the triangulation | `Poly_Triangulation` |
| Mesh simplification and decimation | — |
| **Mesh to B-rep**: plane, cylinder, cone and sphere fitting; region growing; surface reconstruction | — |

## 13. Drawing generation · `ogeom-hlr` · P5

Hidden line removal, exact and polygonal; visible and hidden edge
classification; section and broken-section views; silhouette and isoparametric
curve extraction; reflect lines.

*Elsewhere:* `HLRAlgo`, `HLRBRep`, `HLRTopoBRep`, `HLRAppli`.

Producing a 2D drawing from a 3D model is a core CAD capability. A kernel that
omits it is not a complete alternative, regardless of how few applications reach
for it directly.

## 14. Sketching and constraints · `ogeom-sketch` · P6

A 2D geometric constraint solver: coincidence, distance, angle, parallel,
perpendicular, tangent, symmetry, equality, radius, horizontal and vertical;
construction geometry; driving and driven dimensions; degree-of-freedom
analysis; under- and over-constraint diagnosis that *names the conflicting
constraints* rather than merely failing.

*Elsewhere:* **nothing.** No mainstream B-rep kernel ships a constraint solver,
so every application built on one supplies its own. A kernel that intends to be a
complete foundation should not push this back onto every consumer.

## 15. Document, assembly and product structure · `ogeom-doc` · P5

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

## 16. Selection and display support · `ogeom-select` · P6

BVH construction over tessellated and analytic geometry; ray picking with
sub-shape granularity; rectangle and polygon selection; depth sorting; level of
detail; a stable mapping from triangles back to the topology that produced them.

*Elsewhere:* `SelectMgr`, `Select3D`, `BVH`.

Not a renderer. `tools/ogeom-view` consumes this; it is not part of it.

## 17. Data exchange · `ogeom-io` · P5

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

## 18. Feature recognition · `ogeom-recognize` · P6

Recognizing holes, pockets, slots, bosses, fillets and chamfers from raw
topology; deriving feature trees from imported dumb solids; manufacturing
feature extraction.

Holes, pockets, slots, fillets, chamfers, bosses and partial rounds
recognize today, and `feature_tree()` orders them into ancestry — all
against a parity oracle that asks the solid itself. Manufacturing feature
extraction — the recognized tree mapped to machining operations — remains,
and lands with M11.

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
   places where exactness was being thrown away in `ogeom-math`.*

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
| ~~Analytic mass properties for analytic surfaces~~ **done** | §6 | Volume, area, centroid and inertia integrate on the *exact* surfaces where the trims allow — chart rectangles on any analytic surface, discs on planes — by panelled Gauss quadrature over the chart, exact to rounding for the trigonometric integrands, `deflection` reported as zero. Every primitive lands on its closed form (a sphere's 2/5·m·r² included); scaled placements and free-form trims fall back to the tessellation with its stated chord. |
| ~~`same_parameter` repair~~ **done** | §5 | `repair_same_parameter` samples every pcurve against the edge's own curve at matched parameters — the triangulator's linear range mapping, seam sides included — and either confirms the claim or widens the edge's tolerance to what was measured, so afterwards the flag means something. A native primitive verifies with nothing widened; the imported NIST part widens where its fitted pcurves carry the file's slop, and every claim holds after. |
| ~~`PolygonOnTriangulation` edge representation~~ **done** | §12 | `tessellate` attaches each edge's path through every adjacent face's stored triangulation as node indices, chosen so consecutive indices are triangle edges — position alone cannot decide at a seam, where two chart columns lift to the same points. Serialized in the native format alongside the meshes it indexes, replaced together with them on retessellation. |
| Predicates routed through the `Predicates` trait, everywhere | §7 · M7 | The seam is real where it decides something *combinatorial*: `inside_boundary` — which decides whether a triangulated face keeps a triangle or drops it — goes through `orient2d` with `inside_boundary_with` naming the implementation, and the boolean arrangement's even-odd containment now goes through the same predicate, the crossing question answered by `orient2d` instead of a computed intersection abscissa. Pinned by a test that finds cases where the two implementations genuinely disagree, and a second showing the exact one is right down to one unit in the last place. What is *not* routed, by argument: the ray/triangle test in classification and the tolerance comparisons in mass properties are measurements against a tolerance rather than sign decisions, and forcing them through a predicate would be ceremony. What remains for M7: the marcher's shallow-crossing gate decided by the interval filter rather than a fixed floor. |
| ~~Collapsed wedges~~ **done** | §5 | A zero top extent builds the topology it is: a ridge wedge or a pyramid, five explicit planar faces over shared edges through the `faceted_solid` generalization of the box path, volumes checked against the prismatoid closed form. |
| ~~Compsolid builder~~ **done** | §5 | `make_compsolid` builds it and enforces what the word means: the members must glue into one connected whole through shared *face nodes* — the sharing a boolean or a sew produces. Solids that merely touch, faces coincident but distinct, are a compound's arrangement and are refused as one. History included, like every builder. |
| Non-manifold topology (§4) | §4 | The model permits an edge bounding more than two faces — nothing rejects it — but nothing has been built that way and no test pins the behaviour. `is_shell_closed` counts edge *uses* and asks for an even number, which is the right rule for non-manifold input, but that is a design choice not yet exercised. |
| ~~Deflection in parameter units~~ **done** | §12 | `discretize_on_surface` measures the chord in space through the surface — the lifted midpoint against the lifted chord — so one tolerance means one thing whatever the chart's scale, and the triangulator's no-3D-curve fallback now runs through it. `discretize_planar` keeps its parameter-unit meaning for callers that genuinely work in the chart. |
| ~~Smoothness across a fitted closed loop's seam~~ **done** | §6 | `fit_points_closed` makes a loop's join C1, the constraint eliminated exactly inside the least-squares solve; opt-in, because it spends shape freedom `fit_points` still promises. |
| ~~Points to a B-spline *surface*~~ **done** | §6 | `fit_surface_grid`: rows then columns on shared knot vectors at fixed parameters, the error measured surface-against-every-input. Scattered data still needs a grid first. |
| ~~General affine transforms over a whole *shape*~~ **done** | §5 | `general_transformed_shape` rides the whole-shape conversion: convert exactly, then move control points — which an affine map does exactly, and which leaves every parameterization and every derived pcurve untouched. A sheared box keeps its unit-determinant volume; a sheared cylinder — the case placements cannot express — still encloses exactly its volume. |
| ~~Whole-shape NURBS conversion~~ **done** | §5 | `to_nurbs` walks a solid and restates everything against the moved parameterizations: surfaces bounded to each face's own chart region before converting (a plane's declared domain would map the face to a dot), edges converted exactly with their ranges moved to the new domains, and pcurves re-derived — *exactly* where the edge ran along a chart direction, because the boundary conversion and the patch direction share one parameterization by construction, and by projection fitting at the edge's own parameters otherwise. Seams become the clamped chart's two edge columns. A box measures itself exactly; a cylinder converts seams, rims and all. |
| ~~Exact patches for revolved surfaces~~ **done** | §3 | `revolved_patch` revolves the exact rational profile by the exact rational unit circle, weights multiplied; sphere, torus and a revolved line verify to 1e-12 against independent authorities. Trimmed still refuses. |
| ~~Sweep surfaces proper~~ **done** | §11 | `make_loft_skinned` and `make_pipe_skinned` build skins over the grid fit — seams closed exactly through pinned row ends, pcurves as iso lines — with volumes checked against the frustum and Pappus. Closed free-form spines remain unspoken. |
| ~~The general offset rebuild~~ **done** | §11 | Vertices solve their faces' planar displacement constraints and Newton-polish onto the moved surfaces; edges re-derive from the moved pairs, circles re-framed on old axes so parameters carry; all five analytic surfaces offset. Free-form faces remain refused. |
| ~~Collapsed 2D offsets~~ **done** | §11 | The raw offset splits at its self-crossings, sub-pieces stand or fall by distance to the source, survivors reconnect into loops; open wires offset into capped thick-path outlines. |
| ~~2D corner blends with curved sides~~ **done** | §10 | Tangent-circle centres found on the sides' offset loci, candidates qualified GccAna-style; chamfer distances on an arc are arc lengths. Free-form sides remain refused. |
| ~~Concave blends~~ **done** | §10 | Blends dispatch on convexity read from the face's own trim; the additive wedge fuses where the subtractive one cut; the four revolved rim seats are one parameterization over two signs. |
| ~~Canonical surfaces from a sweep~~ **done** | §9 | `make_prism` builds straight walls as planes on the extrusion's own chart, so the boolean's plane-to-plane resolution sees them; a chamfer melts into a prism wall. |
| ~~Canonical surfaces from a revolution~~ **done** | §9 | `make_revolution` builds a radial line's sweep on the plane it turns in — no seam, no degenerate centre; a revolved rectangle matches `make_cylinder` count for count. |
| ~~Asymmetric conic domains in the native format~~ **done** | §17 | `HyperbolaCurve::over` and `ParabolaCurve::over` take arbitrary increasing domains, and the native reader accepts whatever range the file states instead of refusing asymmetry; write–read–write reproduces the bytes. |
| A consumer for the half-space | §8 | `make_half_space` builds the solid on one side of a face, orienting the boundary away from the material. Nothing reads it yet, and nothing can: its shell is one face with free edges all round, so it is not closed, and every query needing an inside refuses rather than guessing — which is the correct answer for a boundary that does not close. A half space is an argument to a boolean, and the boolean is §8. Its bound comes back *empty* for the same honest reason: `surface_bounds` declines to bound an unbounded plane rather than quoting the declared extent. |
| Curved regions in mesh → B-rep | §9 | Planar recovery is done and exact: coplanar triangles are grown into regions, their boundaries chained into loops, and faces rebuilt on shared edges so the shell closes — a box round-trips through a mesh and measures the same. A *curved* region is refused. Fitting a cylinder to a band of triangles is easy; deciding that the band **is** a cylinder rather than a smooth patch resembling one is the difficulty, and a wrong answer gives a solid that looks right, measures nearly right, and has the wrong surface under every later operation. The crease angle that separates a model edge from a tessellation seam is exposed as a parameter with a stated default, because no fixed angle is right for every mesh — a coarse tessellation of a large cylinder turns further than a fine one of a small cylinder. Deciding it from the mesh's own deflection is canonical recognition. |
| Branch points in the marching intersector | §7 | Where the two normals are parallel the intersection has no single direction, and `trace` stops and says `Stalled` rather than marching through — pushing on is how a tracer changes branch without noticing, and a wrong branch is a plausible answer to a different question. What is owed is *stitching*: detecting the branch point, and joining the fragments that meet there into the curves they belong to. Until then a configuration with a genuine singularity comes back as fragments that each say why they stopped. |
| Boolean coverage beyond the analytic pairs | §8 | The pipeline runs on any faces §7's intersectors answer for, planar and quadric alike, through one general path — the drilled box is the pinned proof. What remains named: same-domain contact (coplanar or co-surface faces) refuses rather than unifies; tangential contact refuses where the touch lies within both faces' reach; scaled placements refuse because a scale changes a surface's parameterization out from under its pcurves; and a marched pair that resolves no branch refuses rather than guesses. The refusal gates are generous by design — a face's reach is its boundary's bound plus a bulge allowance, and a curved face's allowance is most of its own diagonal, so an innocent configuration can be refused but a degenerate one cannot slip through. |
| Same-domain faces in the boolean | §8 | Handled for planar pairs: two faces on one plane split each other by the *other* face's boundary edges — coincident surfaces offer no section curve, so the contact's edges are the splitting curves, projected exactly into the shared chart — and each On piece resolves by comparing outward normals with its coincident partner. Opposed materials vanish from a fuse and survive a cut on the first argument's side; aligned materials keep one copy for fuse and common and none for cut. Stacked boxes, partial stacks and flush walls are pinned. Still refused by name: curved same-domain pairs whose edges have no closed-form projection into the shared chart, and contact reaching only an edge or a vertex, where there is no coincident partner *face* to compare sides against. |
| Tangential contact along a curve | §7 | A plane resting on top of a torus touches it along a whole circle, and there is no transversal curve there to march: the crossing angle is zero along the entire contact. The analytic layer names this contact where it has the closed form (a plane tangent along a cylinder's ruling comes back `Along`), but for the pairs the marcher owns, seeds near the contact converge — the two surfaces sit within the correction's acceptance of each other — wander briefly, and stall. What comes back is fragments hugging the contact curve: on both surfaces to rounding, describing nothing, each saying `Stalled`. The pinned test `tangency_along_a_circle_produces_fragments_not_a_curve` holds that behaviour still. The honest answer needs tangential contact traced as its own kind of curve — following the valley of the gap function rather than a crossing — which is a different walker than the one built. |
| Crossings shallower than a microradian | §7 | The marcher refuses a crossing whose angle is below `SHALLOWEST` (sine 1e-6), and the number is set by the correction rather than by taste: the Newton step may leave a residual up to the confusion tolerance, and on near-coincident surfaces that residual masquerades as an angle between the two computed normals. A gate below that floor reads the correction's own noise as a direction — identical spheres came back as six confident curves that existed nowhere but in rounding. A genuine sub-microradian crossing is also one the correction could not have followed, so the refusal costs nothing that was available. Following such crossings needs the higher-precision arithmetic the `Predicates` seam exists to admit. |
| ~~Smoothness of a fitted intersection across its own seam~~ **done** | §7 | `fit_points_joint_closed` carries the C1 join constraint through the shared seven-dimensional solve, and the approximation stage uses it for every closed branch: the section curve and both pcurves cross their own seam with matching tangents, still same-parameter by construction. |
| Exact pcurves for every analytic section | §7 | `intersect_surfaces` derives exact same-parameter pcurves where the projection has a closed form: any conic on a plane, axis lines and full circles on a cylinder, parallels of latitude on a sphere. Where it does not — an oblique ellipse on a cylinder runs as a cosine in parameter space — the pcurve is `None` rather than a fit, because an exact curve with a fitted pcurve would be a curve whose descriptions disagree by an amount nothing on it records. The consumer that needs those is the boolean, which marches the pair instead. |
| Restriction of closed exact sections to surface extents | §8 | Section *lines* are clipped to both surfaces' parameter extents through their exact pcurves, and closed curves wholly outside an extent are dropped — but a closed curve partially outside is kept whole. Cutting it into arcs is the restriction problem, and the restriction that matters is the face's trim, which is the boolean's 2D splitting stage; the surface extent is only the parameterization window. |
| ~~Pcurves for edges on B-spline surfaces~~ **done** | §9/§17 | Projection fitting at the curve's own parameters serves every corpus face; slop up to `confusion*1e6` is accepted and *recorded* on the edge's tolerance. Parts still refusing do so for the face-level reasons below. |
| ~~Welding the cut-imported result~~ **done** | §8/§9 | Edge ends anchor to their vertices within the recorded tolerance, chart points fold onto the ring-continuing periodic branch, and a border-only weld-and-stitch pass runs at the recorded tolerance, floored at a tenth of the chord. The imported cut measures. |
| ~~Corpus faces that still refuse to triangulate~~ **done** | §9 | All eleven corpus parts triangulate, weld fully shut and measure, deterministically. What it took, in order: bands accept a degenerate apex/pole ring and a one-ring cone face synthesises its apex; a ruling's chart angle is measured over the used range, not at a stated location that may be the apex; half-period fold ties break toward not moving, a mis-wound ring unwinds from its last undecided tie, and two opposite-wound rings merge into one band; subnormal chart coordinates round to the zero they mean and chart-degenerate hair triangles are dropped, not refined; a sphere parallel carries a winding; chart angles across a collapsed direction are interpolated from sound neighbours — except through a pole, where the jump is real; face meshes are re-oriented by what they actually walk, majority-vote parity union-find, so inconsistent file flags cannot poison the flux; and the unit context cited by the geometry's own shape representation wins over whichever one a hash map served first. |
| Vertex blends, face-face blends, marching blends, blend analysis | §10 · M9 | The blend vocabulary speaks edges whose envelope is analytic: straight and circular edges, constant and linear-law radii, both signs. Not yet spoken: the corner patch where three fillets meet (the setback vertex blend — a sphere octant in the equal-radius planar case), the blend between two faces that share no edge, edges whose envelope has no closed form (the spine is the intersection of the two offset surfaces, marched, with circular sections skinned — every part of which now exists separately), and a report of a blend's own tangency error. An earlier cleanup dropped this row while its work was still owed; restored, and scheduled under M9. |
| ~~Sequential blends on one body~~ **done** | §11/§7 | The failure was the marcher's: a second blend's tools march against the first blend's faces — a bevel plane against a blend cylinder is an ellipse with no closed-form pcurve — and an empty marching result was reported as a refusal even though the extended surfaces genuinely never meet. No branch found now asks which story it is in: a conservative projected-grid measurement, and a pair that measurably keeps its distance over its stated extents gets the empty section that is the true answer, while a close pair that resolves nothing keeps the honest refusal. A chamfer and a fillet share one block and the volume pins to the arithmetic. |
| ~~Boss recognition, feature trees, partial rounds~~ **done** | §18 | Convexity is now measured against the solid itself — two parity probes standing just off the edge, both-in-material concave, both-in-air convex — with outward normals calibrated per face from the surface's own foot and the parity ray run in a deliberately generic direction, because rebuilt boolean faces wind their wires however they like, orientation flags lie, and axis-aligned rays graze axis-aligned mesh edges. On that oracle: bosses (an all-convex top whose every wall lands concave on the base), partial rounds as their own category — ftc_11's rim tori come back as two 1.5 mm partial rounds and no false fillets — and `feature_tree()` ordering the flat list into ancestry by face adjacency. |
| A shared multi-resolution pick hierarchy | §16 | Exact picking landed: `pick_refined` lets the tessellation find and order the hits, then intersects the ray with the struck face's *exact* surface over a segment bracketing the mesh answer — a coarsely meshed drum's half-millimetre sag refines onto the true radius within a nanometre — and `Hit::refined` says which kind of number the caller is holding when refinement resolves nothing. What remains is level of detail as a structure: one hierarchy serving several deflections, instead of independent scenes each built at one. |
| Sketch conveniences: analytic derivatives, snapping | §14 | Interactive scale landed in two moves. Differencing exploits structural sparsity — each constraint's rows are formed over only the parameters of the entities it names, a fact the model knows exactly, pinned entry-for-entry against the dense reference — so the Jacobian's cost scales with the constraint count. And `drag()` is the per-frame primitive: the pointer enters as a soft objective every real constraint outweighs, the solve runs warm-started, and a release polish re-solves without the pull so even the whisper leaves no residue. Analytic per-constraint derivatives and snapping remain application-facing conveniences the crate does not offer. |
| Exact HLR | §13 | The drawing pipeline is polygonal: model edges and mesh silhouettes classified by occlusion sampling against the tessellation, as fine as the chord and the sampling. Exact HLR — curve/surface interference resolved analytically, silhouettes as exact curves on the surfaces — is the other half of `HLRBRep`, and it rests on §7's intersectors plus a curve/surface interference walk that does not exist yet. Isoparametric extraction and reflect lines land with it. |
| Sections through tangential contact | §13/§8 | A section plane exactly through a bore's axis meets the cylinder wall along its rulings, and the boolean refuses the configuration ("kept pieces did not close") — the same tangential-contact class already recorded for §7. The textbook half-section is exactly this cut, so resolving that class buys the drawing too. Off-axis sections through the same bore measure to arithmetic today. |
| Documents: textures, undo/redo | §15 | The document holds products, instances, colours, names and semantic PMI, round-trips through STEP, and persists natively byte-stable. The attribute layer landed: user-defined properties (text, number, flag — accumulating by name with replacement), materials with density and an optional colour assigned per shape by id, layers with visibility a shape may sit on several of, and validation properties whose `agrees_with()` compares another measurement at a relative tolerance so a receiver can verify a translation lost nothing — computing the numbers stays with the code that owns the geometry. Still owed from §15's own table: textures, and undo/redo over a transactional model. |
| PMI beyond the semantic core | §15/§17 | Read and written: dimensional characteristics with values and plus/minus bounds — sizes and locations, linear and angular by their own spellings — geometric tolerances with magnitudes, datum precedence, material-condition modifiers (the complex tolerance form, both directions) and composite datum references (a hyphen-joined `A-B` entry, bound into a compartment on write and resolved back out of one on read), datum letters. Locations keep their two feature groups through the round trip, and the deeper aspect walk resolves the NIST angle to its faces. Still not carried: presentation PMI (the annotation polylines and planes) and datum targets. |
| Exchange formats beyond the written set | §17 | DXF writes: R12 polylines on VISIBLE and HIDDEN layers with their linetypes, taking bare polylines so drawings, sections and sketches all serve without the io crate knowing them. glTF (as GLB: one buffer, one node per mesh, POSITION bounds, a metallic-roughness material where a colour was given), OBJ (the file-wide one-based index space) and PLY (ASCII, normals always, per-vertex colour when given) write on the same philosophy, taking bare tessellations. The §17 table's remainder — the conventional kernel's B-rep text format, IGES, 3MF, VRML, *reading* any of the mesh and drawing formats, and the documented proprietary tail — is not started; it is the whole of M12. STEP runs both directions with assemblies, colours, names and semantic PMI; STL both directions; the native format both directions, document layer included. |
| ~~Spindle-torus signed distance~~ **done** | §2 | The sign names the apple — the region the torus-as-a-solid occupies, bounded by the outer sheet, containing the lemon — and the magnitude is measured to whichever sheet is nearer, the folded inner branch included, so a point on either sheet reads zero. For a ring torus the folded branch is never the nearer one and the function reduces to the classical tube distance exactly, pinned point by point, so every consumer keeps the values it had. |

---

## Explicitly out of scope

| | Why |
|---|---|
| Rendering | A kernel supplies geometry and picking, not a render loop. `tools/ogeom-view` is a *software* renderer to an image, and it is a verification tool rather than a feature: it exists so a wrong result is visible, it has no dependency beyond the kernel, and `check.sh` can therefore run it. A real-time windowed viewer would wrap it and add interaction, not correctness. |
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
| `curvo` | NURBS evaluation, curve/curve intersection, 2D region booleans and trimming. The 2D booleans are close to what §8's face-splitting stage needs, so this is the one on the list worth genuinely re-examining — **at the §7 gate**, not before. Taking it earlier would mean two NURBS representations in the tree while `ogeom-geom` is still being shaped |
| `truck-geometry` / `truck-topology` | Evaluated as a *reference design* rather than a dependency. Taking a whole topology crate means taking its data model, and `DATA_MODEL.md` diverges from the conventional one deliberately in two places (§8, §9). Those divergences are the point |

**A correction.** An earlier note in this project's planning claimed that
`robust` was 2D-only and that a separate crate would be needed for `orient3d`
and `insphere`. That is not true of the version in use: it provides all four
predicates, and `ogeom-core`'s exact implementation calls them directly. Recorded
so the claim is not re-derived from the old note.
