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

Still owed inside P2, none of it blocking the above: NURBS conversion, general
affine transforms, mesh simplification, mesh → B-rep.

### M3 — a CAD kernel · P3 · §7–9

- **The gate comes first.** Surface/surface intersection is measured against
  analytic ground truth and published benchmark datasets *before* any boolean
  work is committed to (§7). If it does not clear the bar, the project ships at
  M2 as a geometry library rather than spending years on a boolean built over
  an intersector that cannot carry one.
- General fuse over the corpus, with fuse, common, cut, section and split as
  selection predicates over one result rather than five algorithms.
- Healing survives a corpus of real imported files — which is what makes the
  boolean usable on anything that came from outside.
- History and provenance carried through every boolean and compared operation
  by operation, not just on the final shape.

### M4 — a manufacturing-capable kernel · P4 · §10, §11

- Constant- and variable-radius fillets, chamfers, shelling, offsets, pipe
  sweeps and lofts, on a mechanical part end to end.
- Parity with the field is the bar here, not perfection: these are the two most
  fragile areas of every kernel that has them.

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
| Point-in-solid classification | `BRepClass3d` | done; exact after §7 |
| Minimum distance and proximity between shapes | `BRepExtrema` | after §7 |
| Extrema between curves and surfaces | `Extrema` | after §7 |
| Validity checking against the model invariants | `BRepCheck` | |
| Curvature analysis; inputs for zebra, draft and thickness analysis | `BRepLProp` | |
| Self-intersection detection | `BOPAlgo_CheckerSI` | after §8 |
| Arc length, arc-length parameterization, deflection-based sampling | `GCPnts` | |
| Interpolation and fitting: points → curve, points → surface | `GeomAPI_PointsToBSpline` | |

## 7. Intersection · `og-intersect` · P3

The hardest area in the project, and the one with no existing Rust
implementation.

| | *Elsewhere* |
|---|---|
| Curve/curve, 2D and 3D | `IntCurve` |
| Curve/surface | `IntCurveSurface`, `IntAna` |
| **Surface/surface** — analytic special cases first, then general marching with an approximation stage | `IntPatch`, `GeomInt`, `IntWalk` |
| Polyhedral pre-filtering | `IntPolyh` |
| Contour and silhouette computation | `Contap` |
| Degenerate handling: tangential contact, coincident surfaces, poles, seams | — |

**Gated.** Surface/surface quality is measured against analytic ground truth and
published benchmark datasets before the boolean pipeline is committed to. If it
does not clear the bar, the project ships as a geometry library rather than
sinking years into a boolean built on an intersector that cannot carry it.

## 8. Booleans · `og-bool` · P3

One algorithm — general fuse — plus selection predicates over its result.

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

*Elsewhere:* `ChFi2d`, `ChFi3d`, `Blend`, `BRepBlend`, `BRepFilletAPI`,
`FilletSurf`.

## 11. Offsetting, sweeping, features · `og-offset` · P4

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
| Exact point-in-solid | §7 | Ray casting against the tessellation. Returns `On` for any point within the deflection of the boundary rather than assigning it a side. |
| `same_parameter` repair | §5 | The flag is carried and set false whenever a representation is added, which is honest but pessimistic: every primitive's edges claim disagreement they do not have. The repair routine has to verify agreement, and where it fails, widen the edge's tolerance until the claim is true. Until it exists, nothing may *rely* on the flag being true. |
| `PolygonOnTriangulation` edge representation (`DATA_MODEL.md` §6) | §12 | A face carries a `Triangulation` and an edge carries a `Polyline` with its parameters, which is what makes the stored tessellation watertight. The missing piece is the edge's path through a *specific* face's triangulation as node indices — what a renderer wants to draw a shared edge without hunting for coincident vertices. |
| Predicates routed through the `Predicates` trait | §7 | The trait and both implementations exist in `og-core` and nothing calls them yet. Tessellation, classification and mass properties all use plain `f64` comparisons against tolerances. Retrofitting is mechanical *because* the trait was designed in first — that was the point — but it has not been done, and until it is, "predicates are swappable" is a claim about the design and not about the code. **Scheduled at §7, not later.** The first place swapping the numerics strategy actually pays is the intersector and its polyhedral pre-filter, and the gate that decides whether the boolean pipeline is attempted at all is measured *there* — so the seam has to be real before the measurement, or the gate reports on one strategy and calls it the answer. Tessellation does not need it: the constrained Delaunay already runs on adaptive-exact predicates through `spade`. |
| Collapsed wedges | §5 | `make_wedge` refuses a top extent of zero. A wedge tapering to a ridge has five faces and one tapering to a point has four; both are different topologies, not this one with a zero in it, and building them through the box path would give them a face with no area. They belong to the prism sweep over a triangular profile. |
| Compound and compsolid builders | §5 | `Model::add_compound` exists at the raw level. There is no `make_compound` alongside the other builders, so compounds get no history and no roles. |
| Non-manifold topology (§4) | §4 | The model permits an edge bounding more than two faces — nothing rejects it — but nothing has been built that way and no test pins the behaviour. `is_shell_closed` counts edge *uses* and asks for an even number, which is the right rule for non-manifold input, but that is a design choice not yet exercised. |
| Deflection in parameter units | §12 | `discretize_planar` measures its chord tolerance in parameter units, not in space, because it has no surface to convert through. Callers wanting a spatial tolerance convert themselves. The face boundary already avoids this by discretizing the edge's 3D curve and sampling the pcurve at those parameters — the fallback path only runs for an edge with no 3D curve at all. |
| Adaptive fitting to an error target | §6 | `interpolate` and `approximate` are done. Neither chooses the number of control points for you, and neither stops when the fit is good enough — the caller names the count. Doing it properly means measuring the deviation, inserting a knot where it is worst, and repeating to a tolerance. Left out rather than guessed at: a fit that silently picks its own resolution and reports success is the shape of answer that gets trusted and should not be. |
| Points to a B-spline *surface* | §6 | Curve fitting is done; surface fitting is not. It is the same solve applied twice — rows then columns — over a grid of points, and needs the grid to be a grid, which arbitrary scattered data is not. |
| General affine transforms | §5 | `transformed` covers every rigid motion and uniform scale, and does it by adding a placement rather than touching geometry. A shear or a non-uniform scale is a different type — `GeneralTransform` — which `transformed` does not accept, so the mistake cannot be made. Applying one means rebuilding every curve and surface in its image, including the conversions that make an ellipse out of a transformed circle. Not written; the type refuses the input rather than a stub accepting it and failing. |
| Canonical surfaces from a revolution | §9 | `make_revolution` builds every face it makes as a surface of revolution, including the ones that are really planes: revolving a rectangle with one side on the axis gives a cylinder whose *caps* are polar coordinates on a plane rather than planes. They are seamed accordingly, so the solid has the same face count as `make_cylinder`'s but more edges — a seam and a degenerate centre edge per cap. Nothing is wrong with the result; it is the same solid, validly described. Noticing that a revolved line is a plane is canonical recognition, which belongs to healing. Reproduction: revolve the profile from `a_rectangle_with_a_side_on_the_axis_revolves_into_a_cylinder_face_for_face` and count edges — seven against `make_cylinder`'s three. |
| Exact axis clearance for a revolution | §7 | `make_revolution` refuses a profile that crosses the axis, because revolving one sweeps a surface through itself and the resulting solid has a finite, plausible volume that counts part of space twice. The test is a sampled one: thirty-two places along each edge, looking for the radial direction reversing rather than for the distance to the axis reaching zero, since a crossing almost never lands on a sample. It catches every crossing coarser than its step. Deciding exactly where a curve meets a line is the intersector's work. |
| Unit scale in the native format | §1 | The `.og` file records no unit scale, because a `Model` does not hold one — tolerances are passed to each call rather than carried by the document (`DATA_MODEL.md` §9 says the scale is explicit, and it is, at every call site). So a document authored in metres reads back correctly but is *validated* on the way in against millimetre thresholds, and geometry that is legitimate at one scale could be refused at another. Reproduction: nothing yet fails, because every constructor the reader calls validates against `confusion`, which is scale-relative in `Tolerances` but fixed at millimetres in `native::read`. The fix is a unit scale on the model, written and read with it. |
| Asymmetric conic domains in the native format | §17 | A hyperbola or parabola is written as the symmetric `[-extent, extent]` its constructor produces, and the reader refuses anything else rather than silently re-centring it. Unreachable today — nothing builds an asymmetric one — and it is a refusal rather than a wrong answer, but the format has no way to express one if something later does. |
| A consumer for the half-space | §8 | `make_half_space` builds the solid on one side of a face, orienting the boundary away from the material. Nothing reads it yet, and nothing can: its shell is one face with free edges all round, so it is not closed, and every query needing an inside refuses rather than guessing — which is the correct answer for a boundary that does not close. A half space is an argument to a boolean, and the boolean is §8. Its bound comes back *empty* for the same honest reason: `surface_bounds` declines to bound an unbounded plane rather than quoting the declared extent. |
| Spindle-torus signed distance | §2 | `Torus::distance_to` handles every torus, including the folded branch of a spindle. `signed_distance_to` is documented as ring and horn tori only: a spindle torus encloses two regions and "inside" does not name one of them. |

---

## Explicitly out of scope

| | Why |
|---|---|
| Rendering | `tools/ogview` uses wgpu directly. A kernel supplies geometry and picking, not a render loop. |
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
