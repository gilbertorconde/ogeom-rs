# Changelog

The fifteen library crates are one kernel and move together: they share a
version, and an entry here covers all of them. `tools/` is not published and
is not recorded here.

Dates are the release date. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/), which before 1.0 means a minor
bump may break the API and a patch bump may not.

## [Unreleased]

### Added

- **IGES reads more of the 1980s.** A conic arc (104) of any axis-aligned
  kind reads as its own curve — hyperbola and parabola alongside the
  ellipse, the arc's ends read off the file's points whichever way round
  they are listed. A ruled surface (118) reads as the degree-one patch
  between the two curves' exact spline forms, raised to one degree and
  refined to one knot vector, the second curve walked backward where the
  direction flag says. A constant offset curve (130) and an offset
  surface (140) read as this vocabulary's offsets, the file's direction
  conventions carried by the sign. Parametric spline surfaces (114) and
  rotated conics stay refused by name.
- **Gauss–Kronrod quadrature, and an integral over a rectangle of
  parameters.** The adaptive integral estimated its error by halving and
  comparing, thirty evaluations an interval; it now takes the seven-point
  Gauss and fifteen-point Kronrod pair, fifteen evaluations shared, the
  gap between the two the error. `integrate_2d` integrates over a chart
  rectangle the same way, in tensor form, halving a cell along whichever
  direction the pair finds rougher, so a crease running across the chart
  costs a line of cells rather than a field of them.
- **A stated-tolerance fit for what has no exact spline, and a degree
  limit.** `Curve::to_bspline` still refuses a helix, an offset curve or a
  surface curve rather than approximate silently; `Curve::fitted_bspline_over`
  is the explicit ask, fitting the curve at its own parameters to a named
  tolerance, measured between the samples as well as at them, the error
  reported as found. `BSplineCurve::restricted_to_degree` and
  `BSplineSurface::restricted_to_degree` bring a curve or a patch under a
  degree limit the same way, which is what an exchange format with limits
  needs written.
- **B-spline curves and patches extend.** A B-spline used to end where
  its net ended: an analytic surface's window widens over its carrier,
  but a patch had nothing past its last row. `BSplineCurve::extended`
  and `BSplineSurface::extended` continue a curve past either end, or a
  patch past any of its four sides, by about a length in space: the
  polynomial continuation of the curve's own end derivatives to the
  order asked, raised to the degree and joined on — a rational arc
  continued at order two stays on its circle, and a cylinder patch
  continued round its circle or along its axis stays on its cylinder.
  `widened_to_hold` continues a patch, side by side, as far as a point
  stands off the side its projection clamped to, which is what the
  neighbour-extension steps need on spline faces.
- **Marched blends on cone, sphere and torus hosts.** A seat with one
  host a cone, a sphere or a torus — a bore through a cone's wall, a hole
  drilled through a ring's tube, a ball drilled off its centre — was
  refused by name; the chart inversions those surfaces have in closed
  form are carried now, and a band whose rings start on different columns
  gets a fitted connector where only a cylinder had an exact one. A full
  circular rim whose hosts are not a planar cap and its coaxial wall — a
  bore straight down a ball's axis — takes the march too, where the
  revolved blend had refused it by name. Under it, two fixes any seat
  could hit: a looping seat whose join is a
  corner — the two arcs of a boolean's seam joined end to end — steers
  the march by a smooth refit of itself, since a march cannot cross a
  corner in the guide's derivatives; and the band's first station is
  re-solved exactly on the apex column, where a station a fraction of a
  stride off the host's own seam left a sliver the boolean could not
  hold.
- **Draft on any face.** A wall on a raw fitted patch, or a wall of
  revolution about a neutral plane tilted against its axis, used to be
  refused. Both draft now the way a mould-maker drafts: the face becomes
  the ruled surface through its crossing with the neutral plane, each
  ruling the pull direction turned by the draft angle about the
  crossing's tangent — the construction the planar and revolved drafts
  are the closed forms of. The crossing is read off the face's own mesh
  and corrected onto the surface, so it lies inside the face whatever
  the chart does; the support is the ruled surface itself, degree one
  along the ruling, between a cubic fitted through the crossing and one
  fitted through the rulings' tips at the same parameters. The rebuild
  underneath learned what a fitted support needs: a band on one is
  assembled wire by wire, its rings marched against the caps and
  re-seamed at the vertex the seam starts from, its seam the support's
  own iso-curve (`BSplineSurface::iso_u_curve`, new), and a seam's end
  vertex the crossing of that column with the other seat rather than
  the nearest point two seats agree on.
- **A skinned loft caps a non-planar end section.** A wavy rim, on which
  no plane can stand, is capped by a patch skinned from the rim down to a
  point inside it, sharing the wall's ring edge; the loft used to refuse
  it. A planar end keeps its plane.

### Fixed

- **Removing a tangent chain of blends.** A stadium's top rim rounded
  in one call, its four bands named for removal together, failed by
  name: at a tangent junction neither band's crease pierces the other's
  wall — the straight crease grazes the round wall and the round crease
  grazes the flat one — so no corner was placed and the straight crease
  was taken for a wrapping band. The junction is now where two creases
  touch, placed exactly by the cross-section edge the two bands share:
  the foot of that edge on either crease. And a band standing across a
  circular crease's seam is read as one run about its own centre rather
  than its complement, so an end wall grows back its outer half and not
  its inner.
- **A torus band between two parallels drew as two whole tori laid over
  each other.** Each wire of such a face is a single edge winding the
  chart once. Walked one wire at a time, each rim closed on its own
  translate a tube-period over — which is the whole torus cut along that
  rim — and two of those cancel where they overlap: the band drew to
  more area than its torus has. Wound rims are paired into one ring
  before any rim is closed alone, on either periodic axis, so a ball's
  belt between two latitude wires holds the same way. The parity check
  that guards a boundary had been flagging exactly these faces and
  drawing them again, eight times finer, every time; the benchmark part
  that carries two such bands tessellates in a twentieth of the time.
- **A face narrower than the chord shaded as a saw of fins.** Drawn at
  the caller's chord its boundary sags between points by more than the
  face is wide — a thread flank a tenth of a millimetre wide at a chord
  three times that — and every triangle across the width stands off the
  surface by the sag. A face narrower than four chords has its edges
  drawn to a quarter of its width, and the faces across those edges draw
  them the same. The whole-shape pass reads every face's width off its
  rings before drawing anything, so the big faces round a thousand small
  fillets are drawn once; on a real assembly the area standing more than
  45° off the surface falls from 540 to 80 mm², for 7% more triangles
  and a quarter more time.

- **A vertex at a cone's apex, a sphere's pole or a patch's collapsed
  corner had no normal.** The surface has none there, and the mesh carried
  a zero, which shaded the tip black and dragged every normal it was
  welded with towards nothing. The vertex now carries the limit normal
  along its own column — the apex of a cone seen up a ruling is that
  ruling's normal — found a step inside the domain. On a real assembly,
  1,199 such vertices; none now.

## [0.2.0] — 2026-09-22

A minor bump before 1.0, which is to say it breaks one thing: the default
angular deflection is half a radian. Everything else is additions and
fixes, most of them to meshing, measured on a community printer assembly
of 1,563 solids: every one of them meshes watertight at the default
deflection, with no triangle standing off its surface, in a third of the
time 0.1.0 took.

### Added

- **`heal::fix_shape`**, one pass over a shape nobody promised was
  well-formed: diagnose; put a wire's edges end to end where an order
  exists; collapse an edge shorter than its own vertices' tolerances;
  give a pcurve-less edge the trim projection can honestly fit; sew a
  compound of loose faces or an open shell; tighten tolerances; diagnose
  again. `Fixed` carries the result, its history, and a `FixReport` of
  what was done and what the checker still sees. Small faces and small
  solids are not removed — that is defeaturing — and the report says so
  by leaving them in `after`.

- **`Surface::curvature_at`**: the principal curvatures and directions at
  a point, signed against the surface's own normal, with the mean and
  Gaussian curvatures they combine to and whether the point is an
  umbilic — what curvature display and zebra analysis ask, from the two
  fundamental forms of the jet the trait already carried.
- **Ruled lofts take skew walls.** A wall between two segments that are
  not coplanar — a square lofted to the same square turned an eighth of a
  turn — is the bilinear patch through its four corners, the ruled
  surface between them and exact, where it was refused by name and
  routed to the skinned loft.

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
- **A long bore drew square between its cross holes.** A cylinder never
  sags along its axis, so its grid got one interior row, and the Delaunay
  triangulation bridged from each rim to it with triangles a quarter turn
  wide — under the three chords the repair pass fires at, so they stayed.
  Grid cells are held to a bounded aspect: rows close enough, measured
  through the surface, that no triangle reaches across more than a few
  columns. Fewer triangles on the bore than before, not more.
- **Faces shaded as quilts of creases.** Delaunay in the chart is not
  Delaunay on the surface when the chart's units differ by axis — a
  cylinder's `u` in radians against its `v` in millimetres, a fitted
  strip's `u` over a fiftieth of a unit against a `v` over one — and the
  slivers it makes across the narrow way lift folded, flat across a bend
  the surface takes in between, their normals pointing where none of
  their vertices' do. Three things, measured on a real assembly by the
  area of triangles standing more than 45° off their vertices' normals,
  61,000 mm² before and 540 after:
  - the triangulation runs in the chart scaled to the surface's own
    metric, the mean tangent length each way, so Delaunay sees distances
    as space does;
  - an interior grid point keeps a third of a cell clear of the boundary,
    where a point hugging a boundary chord makes a sliver that stands off
    the surface as a fin;
  - the angular deflection is not asked of a segment both shorter than
    the chord tolerance and a sixteenth of its edge — a fitted edge's
    end hook, a few microns long, which bisection chased down to the
    resolution of the parameter and handed the face a fan of hairs.
  The same assembly meshes in a third of the time with a fifth fewer
  triangles, the slivers the repair pass used to chase now never made.
- **A bore's inside was coarser than its rims.** Edges are discretized
  to the chord *and* the angular deflection; the interior grid was held
  to the chord alone. A viewer scaling its chord to a body's size gave a
  long extrusion's bore thirteen-sided rims and a seven-sided inside. Grid
  cells are held to the normal's turn as well now, so a surface's inside
  is as round as its boundary.
- **A face whose chart sat far from its origin drew as two triangles.**
  The scale a degenerate chart triangle was measured against was taken
  from the ring's coordinates rather than its span, so a cylinder whose
  axis point the file placed half a metre away had every honest cell
  called a hair. Fixture cut from the file that showed it.
- **At a coarser angular deflection, three more ways a face came apart**,
  each with a fixture cut from the file that showed it:
  - A boundary drawn coarsely enough to cross itself and enclose *nothing*
    — an annulus narrower than its rims' sag — was refused, where one that
    enclosed fragments was drawn again with finer edges. Empty is short
    too; a face is refused only when finer edges still enclose nothing.
  - The sag repair inserted a sliver's centre that was its own apex to the
    last bits, round after round, each a hair on the last; the degenerate
    filter dropped the hairs and left a hole. A repair point that lands on
    a vertex already there is not inserted.
  - A grid point that fell exactly on a boundary segment running
    diagonally across the chart was inserted and split that constraint on
    one face alone — a T-junction against the face across the edge. Grid
    points on the boundary stay out.
- **A crossing buried under interior points, a spike, and slop at the
  closing vertex.** Whether a face's boundary crossed itself was read
  from the count of triangles against boundary points, which over a face
  with hundreds of interior points misses a crossing that costs a
  handful; a slotted cylinder wall drew with two holes more than it has.
  The boundary is triangulated on its own first, where the count is
  exact and told by constraint parity rather than by where a hair's
  centre rounds to. What that exposed: a ring that runs out along an edge
  and straight back is a spike that bounds nothing and comes off; and a
  ring whose last point is its first a file's slop away — the two edges'
  own ends of the vertex they share — closed with a fold over its first
  segment, a crossing a fraction of a micron deep. Merged, every solid in
  the real assembly meshes watertight at the default deflection, 1,563 of
  1,563; at half a radian, 1,562.
- **A closed edge's seam is moved to its vertex.** A fitted loop written
  with its start wherever the fit began, the edge's one vertex millimetres
  along it, was held to the curve's own seam, and the vertex's tolerance
  widened to the miss — 2.18 mm in a real assembly, a reach a solid's
  border weld then used. The reader moves the seam to the vertex: the same
  curve, begun where the edge does. `BSplineCurve::reseamed_at` and
  `bspline::join` are new.
- **`triangulate_face` drew every face twice** since the sliver refinement
  landed. It draws once, and again only where the first came up short.
  Face by face, the assembly above meshes in 9.4 s where it took 18.6.

### Changed

- **The default angular deflection is half a radian**, twenty-eight
  degrees, a circle in thirteen segments — what B-rep kernels have long
  defaulted to — where it was 0.2, eleven degrees and thirty-two. With
  the interior of a face now held to the angular deflection as its
  edges are, the old default cost four times the triangles on every
  cylinder; a real assembly meshes to 4.0M triangles at the new default
  where the old gave 19.5M. Ask for `angular: 0.2` to have what the old
  default drew.
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

[Unreleased]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/gilbertorconde/ogeom-rs/releases/tag/v0.1.0
