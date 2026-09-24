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

- **`project_edge_onto_plane`** projects an edge orthogonally onto a plane
  as an exact curve in the plane's coordinates, for a sketch to constrain
  against. A line comes back a line or a point; a circle or ellipse a
  circle, an ellipse or a segment, its arc range counter-clockwise; a
  B-spline the spline on its projected control points. Any other curve
  is fitted and the fit's error returned.

### Fixed

- **A ball or a ring with holes in it stayed faceted.** Only a ball cut
  by parallel planes, or a ring cut square to its tube, was rebuilt. A
  ball or ring whose boundary only makes holes in it (a boss fused on a
  ball, a ball bored across twice, a ring pierced through its tube) now
  comes back as the whole surface with the holes as inner wires: its axis
  is turned so its poles stand farthest from every hole, and its seam into
  the widest angle the holes leave free.

- **Curved faces meeting along a curve that is no circle stayed
  faceted.** A bar drilled across, or a pipe with a branch, came back as
  hundreds of facets. Where two recognized surfaces meet along any curve,
  the meeting is solved onto both and interpolated once, its image on
  each surface made at the same parameters. A face round its axis between
  two rims of any shape, holed or not, is built with a seam of its own
  placed clear of its holes. Both parts convert exactly at a fine mesh and
  a coarse one.
- **A sliver face left a hole in the mesh of its solid.** A face narrower
  than a millionth of its length was drawn as nothing, and its
  neighbours' meshes did not meet across it. It is drawn as a fan across
  its own boundary points, which its neighbours share.

- **Coarse meshes lost their fillets.** At a printer's usual export
  settings a small fillet is three or four facets across, and it came
  back as flat faces. The first fit now tries a dozen vertices when they
  already turn through a bend, a fit needs half as many samples again as
  its surface has parameters, a failed sample that spanned a bend retires
  only its seed, and a facet spanning two rows of a tight fillet joins it
  when the surface accounts for its lean. A fully filleted box meshed at
  the default deflection comes back as its 26 faces.

- **A cut along a cylinder's or torus's seam refused to close.** A box
  face lying on the seam split the plane along the seam edge in pieces
  the curved face did not share. An open section is no longer mistaken
  for a closed loop and split at its middle, and a plane through a
  torus's axis meets it in its two exact tube circles instead of fitted
  ones.
- **Concave faces sent volumes to the mesh.** The check that faces agree
  about which way is out misread faces such as a three-quarter disc, and
  faces cut along their seam. It now reads each face's side from its
  walked boundary, so these shapes are integrated exactly.

- **Mass properties were a part in a hundred off on elliptic walls at the
  default deflection.** Only faces bounded by a chart rectangle or a full
  circle on an analytic surface were integrated exactly; anything else
  sent the shape to its mesh. A face on an analytic surface, or on a line
  or conic swept or revolved, is now integrated round its own trim, so an
  elliptic pad measures to rounding whatever deflection is asked for.
  Shapes with spline faces are still measured on their mesh, and say so
  in `deflection`. The default deflection no longer claims a part in a
  thousand.

## [0.3.1] - 2026-09-23

A patch release with no API change; `cargo semver-checks` finds nothing
against 0.3.0. Converted meshes rebuild whole spheres and tori, sphere caps
and zones, and bent tubes exactly, and drills through converted parts are
faster and close where they refused. The docs are rewritten shorter and
plainer.

### Added

- **`solid_from_mesh` rebuilds spheres and tori that close on themselves.**
  A whole sphere or torus comes back as one face. A sphere cut by one
  plane comes back as a cap, cut by two parallel planes as a zone, and a
  torus cut across its tube as a bent tube, each bounded by full circles.
  Before, these stayed faceted.
  The thin triangles a mesh often has round a pole, whose normals lean
  far from the surface's, join the region that surrounds them.

### Fixed

- **A drill through a finely faceted part spent minutes pairing
  sections.** Every section on the drill's wall was tried for crossings
  against every other section on it, a curve-curve intersection each, and
  a wall met by a few thousand facets of a converted mesh carries a few
  thousand sections: seven minutes on one part. A crossing between two
  sections that share a face counts only inside the faces they do not
  share, so only sections whose other faces' bounds meet are tried: the
  same part pairs its sections in eighteen seconds.
- **A drill through a face of many holes spent seconds on that one
  face.** A section was tried against every boundary edge of its faces,
  and a flat face bounded by a few thousand edges (a panel with four
  hundred vents) took them all; and choosing a split piece's nine probe
  points compared each of hundreds of thousands of candidates with every
  column already taken. Each boundary edge carries its own bound, and
  one outside the face the section must also lie in is skipped; the probe
  search stops at nine columns. A drill through such a part takes under a
  second where it took twenty.
- **A drill through near-coplanar facets refused to close.** A converted
  mesh merges facets that lie within its coplanar distance of one plane
  into one planar face, whose boundary then stands off that plane by up
  to that distance. A drill's section, lying in the plane, met such a
  boundary only within its own reach and missed it, so a whole face's
  section was dropped; neighbouring sections ended on their shared edge
  microns apart, farther than the drill wall's weld reached; and a
  section's end and the edge's split point were built as two vertices.
  On a plane a section meets an edge within the edge's radius, the weld
  on a face reaches as far as the edges its sections end on, and every
  crossing of a tolerant edge is one junction for both: the part that
  refused cuts to a valid solid.

## [0.3.0] - 2026-09-23

A minor release, for one change a caller can see: `FixReport` gains a
field, and a struct literal that built one no longer compiles. That is the only
break `cargo semver-checks` finds against 0.2.1. Meshes come in: 3MF
packages read, deflated, multi-part and ZIP64 alike, and a triangle mesh
becomes a solid whose planar regions merge into faces and whose curved
regions are rebuilt on the cylinders, cones, spheres and tori they lie on,
a scope addition `docs/SCOPE.md` argues. A boolean on such a solid of
thousands of faces takes a second where it took minutes. Imports hold the
tolerance rules they state, IGES round-trips real parts, and a face drawn
out to its whole plane, or the long way round an ellipse, meshes on
itself.

### Changed

- **`FixReport` gains `tolerances_widened`.** A public field on a struct
  callers can build: a literal that built one names the new field too.

### Added

- **`solid_from_mesh`**, a B-rep from a triangle mesh, with
  `MeshSolidOptions`, `MeshSolid` and `MeshSolidReport`. The topology is
  built from the mesh's own connectivity once its repeated vertices are
  welded, with no sewing search, so a 200 000-triangle mesh converts in
  about a second. Coplanar triangles merge into planar faces bounded by
  their outer loops and holes, and collinear boundary runs into single
  edges: an STL cube becomes six faces and twelve edges. Windings are
  made to agree and face outward, a closed piece inside another becomes a
  void, and a mesh that does not close comes back as open shells, with
  holes, non-manifold edges, dropped triangles and flipped windings
  counted in the report.
- **`solid_from_mesh` recognizes curved regions.** Triangles across which
  the surface turns smoothly grow into regions for as long as their
  vertices lie on one cylinder, cone, sphere or torus, verified at the
  coplanar distance, and each is rebuilt on that surface; its boundary
  with each neighbour is placed on the surface exactly, as a parallel
  circle or a ruling, and a band all the way round its axis gets a seam.
  A meshed bore comes back as a cylinder between two circles, and a
  rounded box as six planes, twelve cylinders and eight spheres. A region
  whose boundary is no such curve, or whose surface is round both ways,
  stays faceted and is counted (`curved_faces`, `curved_faceted`);
  `MeshSolidOptions::recognize` turns it off and `crease` sets the angle
  that counts as an edge. `recognize_points` exposes the recognition
  itself, with its measured deviation as the certificate. The mesh
  conversion is in scope by `docs/SCOPE.md`, which says why.
- **`read_3mf`**, reading a 3MF package into one placed mesh per build
  item, with `ThreeMfImport`, `ThreeMfObject` and `ObjectType`.
  Components are flattened, the production extension's multi-part
  packages followed (every slicer writes one object per part), the
  model's unit scaled to millimetres, a mirroring transform's triangles
  rewound, and a uniform object colour kept; per-triangle colours,
  textures, supports and unknown required extensions come back as
  warnings. Deflated entries are inflated by a decoder in the crate, and
  the archive is read through its central directory, so entries streamed
  with their sizes after the data read too. `read_package` reads deflated
  entries as well, checks every entry's checksum, and reads ZIP64
  archives (which streaming writers emit whatever a package's size, and
  which half the 3MF files downloaded from model sites are), refusing
  encrypted entries and archives spanning several disks by name.
- **`restore_containment`**, the pass that establishes the tolerance
  containment rule the checker enforces: every edge widened to at least
  the faces it bounds, every vertex to at least the edges it bounds,
  only ever growing, returning how many entities grew.

### Fixed

- **A sliver face was drawn as its whole plane.** A planar face narrower
  than a millionth of its length (a strip 163 mm long and 40 nm wide, as
  one CAD export leaves along a mirrored panel) has its two long sides
  merged into one line by the mesher, and no boundary ring survived; the
  mesher then took it for a face without wires and drew the plane's whole
  window, a hundred metres each way, so an imported printer's skirts
  rendered as bars across the scene. A face with wires whose rings all
  collapse encloses nothing, and meshes to nothing.
- **A boolean on a solid of many small faces took minutes.** Rebuilding
  the result compared every strand end with every vertex minted so far,
  and `sew` compared every edge with every other (projecting one's
  middle onto the other's curve for each pair), and every vertex and
  face likewise, so the cost grew with the square of the face count even
  when the tool touched a handful of faces. A thin drill through a
  converted mesh of 6 000 faces took two minutes. Vertices, junctions
  and edge ends are binned by position and compared only with those
  near them, in the order the full scan visited them, so the results are
  the same: the same drill takes under a second, and 12 000 faces two.
- **An edge on a very eccentric ellipse could run the long way round.**
  The STEP and IGES readers placed a vertex on an ellipse by its
  eccentric anomaly, which is exact only for a point on the curve. A
  vertex two microns off an ellipse 6.5 m by 1.8 mm read as nine
  millimetres along it, the end parameter fell before the start, and the
  edge took almost the whole ellipse: `triangulate_face` drew its faces
  out to ten metres. The readers refine the closed form to the nearest
  point on the curve.
- **A STEP read broke the tolerance containment rule, and `fix_shape`
  did not restore it.** The reader widens an edge to how far its pcurves
  sit from its curve, and left the edge's vertices at the confusion
  tolerance: `check` called every such vertex broken, hundreds on an
  ordinary part and thousands on an assembly, and `fix_shape` only ever
  tightened, so it never cleared them. The STEP and IGES readers run the
  containment pass on every body they build, and `fix_shape` runs it
  last, after the reduction.
- **A reshape reversed a wire's walk once too often.** A substitution
  rebuilds everything above the substituted entity, and read each
  occurrence's children with that occurrence's orientation composed in,
  then oriented the new node as the occurrence again. The occurrence being
  rebuilt came out right, but the new node it stored was inverted, and
  every other occurrence of it (a wire shared by a face reached reversed)
  walked its edges tail to head. Children are read through the forward
  occurrence and the result oriented once.
- **Collapsing a degenerate edge opened a gap.** `fix_shape` collapses an
  edge shorter than its vertices' tolerances onto one vertex, and the
  surviving vertex kept its own tolerance, so the neighbouring edges'
  curves stopped the collapsed length short of it. The survivor widens to
  reach every vertex it absorbs.
- **`write_iges` overflowed the eighty-column record, and its own reader
  refused the file.** A near-zero real went out in positional notation:
  a coefficient of 1.5e-51 spelt in sixty-nine characters where a record
  holds sixty-four. Every real takes the shorter of its positional and
  exponent spellings, and a parameter longer than a whole record fills
  every record it crosses, so the reader rejoins it with nothing inserted.
- **`read_iges` refused a solid whose curves missed their vertices by
  nanometres.** IGES states no tolerances, so every vertex was built at the
  confusion tolerance, and a writer's last digit (half a nanometre)
  refused the whole solid; a vertex a hair past a bounded curve's end
  also asked for a parameter outside the curve. A vertex's tolerance now
  grows to state its curves' miss, up to the millimetre past which a
  boundary is not that curve's at all, and the range is held to the
  curve's domain. Every solid in the exchange corpus survives the
  kernel's own IGES write and read, at the volume it went out with.

## [0.2.1] - 2026-09-23

A patch release: no public signature changed in any crate, checked against
0.2.0 by `cargo semver-checks`. It adds, and it fixes output that was wrong.
Fillets reach every host the march can seat on, fitted patches included,
and the corner tool rounds any convex planar vertex. Edges asked together
meet: `chamfer_edges` mitres, and `fillet_edges` closes a shared corner
with the rolling ball's patch rather than leaving the bands' caps. This is a
visible change for any caller that rounds a box corner in one call. A
viewer meshing face by face can now agree with the whole-shape mesh along
every edge, and STEP files from plainer writers and mesh converters read.

### Added

- **The marched blend takes fitted hosts as far as the melt.** A
  B-spline patch, or a swept or revolved surface, is a host the march
  had refused by name. It marches now: the chart inverts by projection
  warm-started from the last station, a patch that meets itself at its
  seam wraps as a periodic one does, the patch is continued past its
  face by a few radii so the ball can run out through a wall, and a
  seat the boolean split into arcs at the hosts' seams is closed back
  into one loop through the neighbours the two hosts share, arc by arc
  where each continues the last tangentially, and marched as the loop it
  is. The wedge is built on the host's own patch, widened in place, its
  rail loop slid into the host's own window, and it melts: a converted
  box's edge, the crease round a converted post, a converted cone's rim
  and a converted drum's rim all round to valid solids. A melt the
  boolean still cannot resolve refuses by name. Under it, `to_nurbs` now carries the
  old vertices' and edges' recorded tolerances into the converted solid,
  which it had dropped, and the boolean inverts a probe on a patch by
  projection and treats two identical patches as one chart.
- **Chamfers bevel a chain of edges as one operation, mitred.**
  `chamfer_edges` and `chamfer_edges_with` (a `Chamfer` per edge, in any
  of the three spellings) build every wedge on the solid as it stands
  before the call and then apply them all, so where two bevels meet at a
  vertex each still reaches the corner and the bevel planes meet along
  their own line. Beveling a box's four top edges one call at a time left
  a small tetrahedron and two extra faces at every corner; as one
  operation it is ten faces and the exact volume. Three bevels at a
  convex vertex meet at one point.
- **The corner tool rounds a vertex no single ball touches.** A
  rectangular pyramid's apex, an irregular pentagonal one, or any convex
  planar vertex whose faces share no tangent ball was refused by name.
  The ball's centre may sit anywhere a radius in from every host plane,
  and at such a vertex that region's tip is a few points joined by short
  ridges; the rounded corner is the exact envelope of the rolling ball, a
  sphere at each tip vertex and a cylinder along each ridge, each cut
  with its own compartment (the sphere by the one-ball tool on its three
  planes, the ridge by the flush fillet of a virtual crease), the
  compartments meeting cap to cap on the planes square to the ridges. A
  ridge seven microns long stays the sliver of cylinder it is. The
  sphere's pole stands along a ridge where there is one, so the ridge
  fillet's cap meets it on a circle the chart images exactly. The flush
  fillets follow at a sharp rectangular apex; where a sphere clears a
  fourth plane by a few hundredths of a millimetre the envelope keeps a
  sliver of that plane, and the flush fillet meeting the sliver still
  dies in the cut.
- **IGES reads more of the 1980s.** A conic arc (104) of any axis-aligned
  kind reads as its own curve: hyperbola and parabola alongside the
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
  order asked, raised to the degree and joined on. A rational arc
  continued at order two stays on its circle, and a cylinder patch
  continued round its circle or along its axis stays on its cylinder.
  `widened_to_hold` continues a patch, side by side, as far as a point
  stands off the side its projection clamped to, which is what the
  neighbour-extension steps need on spline faces.
- **Marched blends on cone, sphere and torus hosts.** A seat with one
  host a cone, a sphere or a torus (a bore through a cone's wall, a hole
  drilled through a ring's tube, a ball drilled off its centre) was
  refused by name; the chart inversions those surfaces have in closed
  form are carried now, and a band whose rings start on different columns
  gets a fitted connector where only a cylinder had an exact one. A full
  circular rim whose hosts are not a planar cap and its coaxial wall (a
  bore straight down a ball's axis) takes the march too, where the
  revolved blend had refused it by name. Under it, two fixes any seat
  could hit: a looping seat whose join is a
  corner (the two arcs of a boolean's seam joined end to end) steers
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
  crossing's tangent. The planar and revolved drafts are the closed
  forms of this construction. The crossing is read off the face's own mesh
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

- **A prism swept from a clockwise profile was inside out.** A closed
  wire on a plane bounds one region however it is walked, but the prism
  read each wall's side off its edge's direction, so a square walked
  clockwise about the travel swept its walls facing into the material
  while the caps faced out, and the volume was refused as wound inward.
  Each ring's walls now follow the ring's winding about the travel, the
  outer ring turning positively and every hole the other way, whichever
  way the caller walked them.
- **Fillets asked together at a corner close it with the ball's patch.**
  Three fillets at a box corner through `fillet_edges` left the bands'
  flat caps standing where a ball rolls round the corner. Where three or
  more edges of one call meet at a vertex, the corner tool now rounds
  the vertex first and the bands stop flush against its patch: an
  octant of a sphere at a box corner, the envelope of spheres and
  cylinders at an apex no single ball touches. A corner the tool does
  not speak keeps the caps. Under it, the corner tool read a vertex's
  point where its node was built rather than where it stands, and at a
  prism's far end it rounded the near corner instead.
- **A band meeting two bands at an apex left its triple point as two vertices.** Three
  descriptions of one point arrived a tenth of a micron apart in turn,
  each within the last's reach and none within the first's, and the
  junctions the arrangement adds after the first merge were never merged
  with each other; they are, an end welds into a vertex within the
  vertex's own tolerance and the weld, and a crossing at an edge's end is
  the end on the edge as well as on the section. What is still owed is
  named in the ledger: a band's section touching an edge describes the
  touch point twice, a tenth of a micron apart, and the third flush fillet
  at a sharp apex still dies on it.
- **STEP faces on a surface of linear extrusion were skipped, and the
  body drew open.** A curve swept along a vector is how some writers
  spell a drum's wall, and the reader skipped every such face by name,
  so a part with a slot's rounded ends on them showed its inside through
  the gaps. A circle swept along its own axis now reads as the cylinder
  it is, a line swept as the plane it is, and anything else as the swept
  surface itself over a window the face's own edges size, its trims
  fitted by projection.
- **An open edge whose vertices stand on its curve short of the ends.**
  A file that writes the whole spline and lets the vertices say where
  the edge stops, millimetres in, had the edge held to the whole curve,
  overshooting its neighbours, and the faces it bounded drew as nothing.
  The window between the vertices' own feet on the curve is taken, and
  the report tallies it as `vertex-window`.
- **A marched section between two patches wandered, and a loop cut at
  a seam stayed open.** The boolean fits a marched section jointly with
  its two chart images, and parameterised the samples centripetally: a
  walk that starts on a window's rim halves its step to a micron and
  grows it back only twofold a point, so most of its samples crowd one
  end, and a cubic through them could not be the straight line they lay
  on: a section came back twelve hundred millimetres off it, and the
  boolean refused the wedge for running beyond a turn. A joint fit that
  misses its target centripetally is fitted again by chord length, and
  the closer of the two stands. A loop
  walked round a converted drum reached the seam from both sides and was
  flagged as having left the domain; its ends coincide and it is closed,
  its chart image unwrapped across the seam as across a period.
- **A fused post converted to patches drew open, and the conversion
  failed the checker.** The second half of a rim split by a fuse was
  written as the chart segment between its projected endpoints, the seam
  vertex projecting to either column, and it drew the front half
  mirrored; the wall's ring lost its far side. A seam endpoint's column
  is decided by the edge's own interior. The fit's slop widened each
  converted edge past the vertices bounding it; the vertices widen with
  the edge.
- **A fitted pcurve hooked off its curve between the fitter's samples.**
  The reader fits an edge's pcurve by projecting samples spaced evenly in
  the curve's parameter and measuring the fit at those samples alone. A
  curve is free to run far faster at one end than the other (a fan
  blade's root meets its hub through most of its bend inside the first
  sample interval), and there the cubic hooked four tenths of a
  millimetre past the curve while every sample sat within a hundredth:
  the ring crossed the face's neighbouring ring in the chart and the hub
  drew with overlapping triangles and open edges at every blade. The fit
  is now measured between its samples as well, through the surface, and
  where it leaves the curve by more than a micron it is asked again at
  twice the samples, a few rounds at most, the closest round kept; the
  reported error and `fit-short` count what was found between the
  samples too.
- **A mesh assembled face by face cracked along every narrow face's
  edges.** A face narrower than a few chords draws its edges finer than
  asked, as does one whose boundary crosses itself at the chord, and
  `triangulate` tells the faces across those edges to draw them the
  same. `triangulate_face` could not be told, so a viewer keeping one
  mesh per face (and the kernel's own stored tessellation, built the
  same way) drew the two sides of such an edge to different points: a
  seam of cracks round every fillet and along every thin plate.
  `edge_chords_for` answers the agreement once per shape and
  `triangulate_face_with` draws a face to it; `tessellate` stores its
  faces and its edge polylines by the same answer.
- **A bore crossed by a hole at its wall drew with a ring missing.** The
  hole's two rims wind the bore's chart between them, and pairing them
  along the first rim's last column could run into the hole's own
  opening, or leave a rim's tail on the far side of the seam from the
  ring it was folded into. The pairing column is chosen in the widest
  gap free of the other rings, either rim is cut whichever way it runs,
  and the rings the band passes are slid into the same stretch of the
  period; a cross hole through a bore wall draws closed.
- **STEP files spelt in older or plainer vocabulary.** A face written as
  `FACE_SURFACE` rather than `ADVANCED_FACE`, a shape carried by a
  `SHAPE_REPRESENTATION` or a manifold surface or faceted representation
  rather than an advanced B-rep one, and a product whose formation is
  written with its source read as they were meant; the source-spelt
  formation's product is taken from its third slot whether it stands
  alone or inside a complex instance, and a product whose name slot is
  blank (a mesh converter fills the id and leaves the name) reads by
  its id rather than its entity number. An edge whose
  vertices stand at descending parameters on its curve whatever its
  sense flag says (every edge written forward and the line left to run
  the other way) is reversed to run with the curve, tallied as
  `edge-against-curve`.
- **Removing a tangent chain of blends.** A stadium's top rim rounded
  in one call, its four bands named for removal together, failed by
  name: at a tangent junction neither band's crease pierces the other's
  wall (the straight crease grazes the round wall and the round crease
  grazes the flat one), so no corner was placed and the straight crease
  was taken for a wrapping band. The junction is now where two creases
  touch, placed exactly by the cross-section edge the two bands share:
  the foot of that edge on either crease. And a band standing across a
  circular crease's seam is read as one run about its own centre rather
  than its complement, so an end wall grows back its outer half and not
  its inner.
- **A torus band between two parallels drew as two whole tori laid over
  each other.** Each wire of such a face is a single edge winding the
  chart once. Walked one wire at a time, each rim closed on its own
  translate a tube-period over, which is the whole torus cut along that
  rim, and two of those cancel where they overlap: the band drew to
  more area than its torus has. Wound rims are paired into one ring
  before any rim is closed alone, on either periodic axis, so a ball's
  belt between two latitude wires holds the same way. The parity check
  that guards a boundary had been flagging exactly these faces and
  drawing them again, eight times finer, every time; the benchmark part
  that carries two such bands tessellates in a twentieth of the time.
- **A face narrower than the chord shaded as a saw of fins.** Drawn at
  the caller's chord its boundary sags between points by more than the
  face is wide (a thread flank a tenth of a millimetre wide at a chord
  three times that), and every triangle across the width stands off the
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
  along its own column, found a step inside the domain: the apex of a
  cone seen up a ruling has that ruling's normal. On a real assembly,
  1,199 such vertices; none now.

## [0.2.0] - 2026-09-22

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
  solids are not removed (that is defeaturing), and the report says so
  by leaving them in `after`.

- **`Surface::curvature_at`**: the principal curvatures and directions at
  a point, signed against the surface's own normal, with the mean and
  Gaussian curvatures they combine to and whether the point is an
  umbilic: what curvature display and zebra analysis ask, from the two
  fundamental forms of the jet the trait already carried.
- **Ruled lofts take skew walls.** A wall between two segments that are
  not coplanar (a square lofted to the same square turned an eighth of a
  turn) is the bilinear patch through its four corners, the ruled
  surface between them and exact, where it was refused by name and
  routed to the skinned loft.

### Fixed

- **Meshing a real assembly closed.** 1,499 of 1,563 solids meshed
  watertight before; 1,562 do now, and none refuse. Six causes, each with a
  fixture cut from the file that showed it:
  - `Triangulation::is_closed` counted an edge's uses and wanted exactly
    two. It now requires every edge to be crossed as often each way, which
    is what the divergence theorem needs: it catches a face wound inside
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
  - An inner loop thinner than a micron (arcs out, fitted splines back)
    is a slit, not a hole, and is no longer handed to the triangulator.
- **A long bore drew square between its cross holes.** A cylinder never
  sags along its axis, so its grid got one interior row, and the Delaunay
  triangulation bridged from each rim to it with triangles a quarter turn
  wide, under the three chords the repair pass fires at, so they stayed.
  Grid cells are held to a bounded aspect: rows close enough, measured
  through the surface, that no triangle reaches across more than a few
  columns. Fewer triangles on the bore than before, not more.
- **Faces shaded as quilts of creases.** Delaunay in the chart is not
  Delaunay on the surface when the chart's units differ by axis (a
  cylinder's `u` in radians against its `v` in millimetres, a fitted
  strip's `u` over a fiftieth of a unit against a `v` over one), and the
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
    the chord tolerance and a sixteenth of its edge, such as a fitted edge's
    end hook a few microns long, which bisection chased down to the
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
    (an annulus narrower than its rims' sag) was refused, where one that
    enclosed fragments was drawn again with finer edges. Empty is short
    too; a face is refused only when finer edges still enclose nothing.
  - The sag repair inserted a sliver's centre that was its own apex to the
    last bits, round after round, each a hair on the last; the degenerate
    filter dropped the hairs and left a hole. A repair point that lands on
    a vertex already there is not inserted.
  - A grid point that fell exactly on a boundary segment running
    diagonally across the chart was inserted and split that constraint on
    one face alone: a T-junction against the face across the edge. Grid
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
  ring whose last point is its first a file's slop away (the two edges'
  own ends of the vertex they share) closed with a fold over its first
  segment, a crossing a fraction of a micron deep. Merged, every solid in
  the real assembly meshes watertight at the default deflection, 1,563 of
  1,563; at half a radian, 1,562.
- **A closed edge's seam is moved to its vertex.** A fitted loop written
  with its start wherever the fit began, the edge's one vertex millimetres
  along it, was held to the curve's own seam, and the vertex's tolerance
  widened to the miss: 2.18 mm in a real assembly, a reach a solid's
  border weld then used. The reader moves the seam to the vertex: the same
  curve, begun where the edge does. `BSplineCurve::reseamed_at` and
  `bspline::join` are new.
- **`triangulate_face` drew every face twice** since the sliver refinement
  landed. It draws once, and again only where the first came up short.
  Face by face, the assembly above meshes in 9.4 s where it took 18.6.

### Changed

- **The default angular deflection is half a radian**, twenty-eight
  degrees, a circle in thirteen segments (what B-rep kernels have long
  defaulted to), where it was 0.2, eleven degrees and thirty-two. With
  the interior of a face now held to the angular deflection as its
  edges are, the old default cost four times the triangles on every
  cylinder; a real assembly meshes to 4.0M triangles at the new default
  where the old gave 19.5M. Ask for `angular: 0.2` to have what the old
  default drew.
- `BSplineSurface` settles whether its net closes at construction, so
  evaluation past a closed join costs what evaluation inside costs.
- CI verifies the declared `rust-version` on every push, reading it from
  the manifest so the two cannot drift.

## [0.1.0] - 2026-09-18

First public release.

### Added

- **Geometry.** Parametric curves and surfaces in two and three dimensions
  (lines, conics, B-splines rational and not, extrusions, revolutions, offsets
  and trims) behind adaptor traits, on a B-spline substrate with knot
  insertion, splitting, degree elevation and Bézier decomposition.
- **Topology.** One shared B-rep model: geometry and topology in arenas, a
  shape a cheap handle into them, the same node placed, mirrored or instanced
  many times without being copied. Per-entity tolerances, location chains,
  orientation composed through the tree, and edges carrying a list of
  representations rather than one curve.
- **Construction and measurement.** Primitives, polyhedra, sewing, shape
  validity, mass properties exactly where a closed form exists and from the
  mesh otherwise, bounds, projection and classification.
- **Booleans.** A general fuse with the filters over it (union, difference,
  intersection, section) plus defeaturing over the same machinery.
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
- **Exchange.** STEP and IGES in both directions (STEP carrying assemblies,
  colours, semantic PMI and saved views) plus the native format, `.brep`,
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

[Unreleased]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.1...HEAD
[0.3.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/gilbertorconde/ogeom-rs/releases/tag/v0.1.0
