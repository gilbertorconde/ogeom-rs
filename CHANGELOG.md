# Changelog

The fifteen library crates are one kernel and move together: they share a
version, and an entry here covers all of them. `tools/` is not published and
is not recorded here.

Dates are the release date. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/), which before 1.0 means a minor
bump may break the API and a patch bump may not.

## [Unreleased]

### Fixed

- **Meshing a real assembly closed.** 1,499 of 1,563 solids meshed
  watertight before; 1,562 do now, and none refuse. Six causes, each with a
  fixture cut from the file that showed it:
  - `Triangulation::is_closed` counted an edge's uses and wanted exactly
    two. It now requires every edge to be crossed as often each way, which
    is what the divergence theorem needs — it catches a face wound inside
    out (which counting passed) and accepts four triangles round one edge
    (which counting refused).
  - A face narrower than the chord error its boundary is drawn with had a
    boundary that crossed itself; the triangulator returned fragments. The
    boundary is redrawn finer, keyed to the edge so both faces sharing it
    agree.
  - A run along a cone's apex row was dropped when its rulings stood
    exactly a quarter turn apart. Whether a gap is an apex run is decided
    by its width and by whether its chart midpoint lifts to the shared
    vertex, not by a fraction of the period.
  - A chart direction was called degenerate by ratio to the other
    direction, so a densely parameterised `v` made every `u` look weak and
    two half-millimetre edges fitted as a point. It is judged by what
    crossing the whole span moves, in length.
  - Ring folding across a chart's join engaged on periodicity only; a
    B-spline tube that closes without repeating never folded. Closure is
    the test now, folded rings are slid back into the domain, and
    `Surface` evaluation wraps a parameter past a closed join instead of
    refusing it.
  - An inner loop thinner than a micron — arcs out, fitted splines back —
    is a slit, not a hole, and is no longer handed to the triangulator.
- **`triangulate_face` drew every face twice** since the sliver refinement
  landed. It draws once, and again only where the first came up short.
  Face by face, the assembly above meshes in 9.4 s where it took 18.6.

### Changed

- `BSplineSurface` settles whether its net closes at construction, so
  evaluation past a closed join costs what evaluation inside costs.
- CI verifies the declared `rust-version` on every push, reading it from
  the manifest so the two cannot drift.

## [0.1.0] — 2026-09-18

First public release.

### Added

- **Geometry.** Parametric curves and surfaces in two and three dimensions —
  lines, conics, B-splines rational and not, extrusions, revolutions, offsets
  and trims — behind adaptor traits, on a B-spline substrate with knot
  insertion, splitting, degree elevation and Bézier decomposition.
- **Topology.** One shared B-rep model: geometry and topology in arenas, a
  shape a cheap handle into them, the same node placed, mirrored or instanced
  many times without being copied. Per-entity tolerances, location chains,
  orientation composed through the tree, and edges carrying a list of
  representations rather than one curve.
- **Construction and measurement.** Primitives, polyhedra, sewing, shape
  validity, mass properties exactly where a closed form exists and from the
  mesh otherwise, bounds, projection and classification.
- **Booleans.** A general fuse with the filters over it — union, difference,
  intersection, section — plus defeaturing over the same machinery.
- **Blends.** Constant and variable radius fillets and chamfers on single
  edges, tangent chains and full rims, closed forms where they exist and a
  marched rolling ball where they do not, and the corner where three blends
  meet.
- **Offsets and sweeps.** Offsetting, shelling, sweeping along spines open and
  closed, lofting through sections, and draft.
- **Tessellation.** Edge discretization to chord and angular tolerances, and
  constrained Delaunay triangulation per face in its own chart.
- **Healing.** Validity diagnosis, sewing, reanchoring periodic rings, and
  instructed fixes for trims a reader refused and boundaries sitting off the
  surface they bound.
- **Drawings.** Hidden line removal, sections and hatching.
- **Documents.** Assemblies and product structure, appearance, PMI, saved
  views, transactions and persistence.
- **Exchange.** STEP and IGES in both directions — STEP carrying assemblies,
  colours, semantic PMI and saved views — plus the native format, `.brep`,
  STL, DXF, glTF, OBJ, PLY, VRML and 3MF.

### Known restrictions

This release is measured against 97 named capabilities in
`docs/PARITY.md`, and the ledger is part of the build: 66 are covered, 16
carry a stated restriction, 9 diverge by design, 5 do not apply and 1 is
unreviewed. `tools/check.sh` fails if the audit and the code drift apart.
The `partial` rows say exactly what each one does not do; the largest are
blend hosts beyond planes and cylinders, the general N-edged setback vertex,
shape healing's wire reordering and small-feature removal, the medial axis
beyond convex polygons, and the IGES entities listed as refused by name.

Nothing here is a silent gap. A capability that is not implemented refuses
by name rather than returning an answer it cannot stand behind.

[Unreleased]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/gilbertorconde/ogeom-rs/releases/tag/v0.1.0
