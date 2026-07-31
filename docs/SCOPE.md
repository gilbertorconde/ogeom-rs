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
| Mass properties: length, area, volume, centroid, inertia, principal axes | `GProp`, `BRepGProp` | after §12 |
| Point-in-face classification | `BRepTopAdaptor_FClass2d` | after §12 |
| Point-in-solid classification | `BRepClass3d` | after §7 |
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
| Native `.og` | r/w | Versioned, forward-compatible, lossless. Lands in P1 as the test interchange format |
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

**Point-in-solid depends on intersection (§7).** Ray casting needs ray/surface
hits. Those are closed-form for the quadrics but not for a NURBS patch, so the
general case waits for the marching intersector. Point-in-*face* needs neither —
it is a winding count over the face's pcurves in parameter space — and lands
with tessellation.

**Distance and extrema between shapes depend on intersection (§7)** for the same
reason: the minimum distance between two curved faces is a constrained
minimization whose machinery is the intersector's.

**Self-intersection detection depends on booleans (§8).** It is the boolean
engine's interference stage run against a single argument.

---

## Explicitly out of scope

| | Why |
|---|---|
| Rendering | `tools/ogview` uses wgpu directly. A kernel supplies geometry and picking, not a render loop. |
| FEA meshing — tetrahedral, boundary-layer | A separate discipline with its own literature. The kernel supplies the geometry and surface meshes it consumes. |
| Simulation, CAM toolpath generation | Applications, not kernel. |
| A parametric *recompute engine* | Application-level. The kernel supplies stable identity and history (`DATA_MODEL.md` §7, §8), which is what such an engine needs from it. |
| ABI or file-format compatibility with any other kernel | Independence is the point. Interoperation happens through published exchange formats. |

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
