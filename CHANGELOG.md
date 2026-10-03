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

- **Fine meshes of free-form surfaces convert to one patch.** Where
  recognition cut a smooth free-form area of a fine mesh into small
  spheres, cones, cylinders and tori, `solid_from_mesh` fits the area and
  the pieces it encloses as one B-spline patch, kept only if the patch
  verifies. A bumped plate meshed at 80 or 120 cells a side comes back as
  six faces instead of twelve to twenty-four thousand; a canonical face
  tangent to anything outside the smooth area is never taken.

### Fixed

- **Pads fuse into slabs drafted under a thousandth.** A planar face
  leaning off another by under a thousandth of a radian, the two sharing
  an edge to their points' rounding, crosses it on that edge rather than
  on the solved line a sliver off it, where the band in which they stand
  within the weld distance is under a tenth of either face. A pad pushed
  into a converted slab drafted 0.0001 to 0.0007 fuses, cuts and is taken
  in common instead of being refused.

- **A round's cap flush with a side face of the solid.** A marched round
  whose cap stands in the plane of a side face (an edge round or face blend
  of a leaning drum's foot running off the block's edge) melts into the
  solid: the band's end row lies in its cap plane and carries its measured
  tolerance. Arcs of one crease split at a seam meet in one cap plane
  instead of overlapping or leaving a sliver, which also raises fills that
  came out short at the seam. Fills tend to the upright closed form within
  3e-6.

## [0.8.0] - 2026-10-03

A minor release, since the API changes: `MeshSolidReport` gains a field,
and the mesh converter's steps are a public API (`MeshRegions`, with
`merge`, `split`, `fit` and `build`). The converter keeps an open mesh's
free boundary on its surfaces, seats a holed torus's seam clear of its
holes, tells a finely drawn curve from scatter and runs 3 to 5 times
faster on large meshes. Three silent wrong answers are fixed: a ball
drilled through its poles, empty boolean results taken as valid, and an
edge round over a drum's foot run off the block. Marched sections state
their measured error, pads fuse into slabs drafted less than their chords
sag, and `intersect_surfaces` answers `Same` for coincident patches.

### Added

- **Mesh conversion's steps as an API.** `MeshRegions::find` returns a
  mesh's regions (triangles, surface, deviation, neighbours); `merge`,
  `split` along a vertex path and `fit` of a `SurfaceKind` with
  `FitConstraints` (a fixed axis or radius) correct them, each verified at
  every vertex, and `build` makes the solid through the same pipeline as
  `solid_from_mesh`, which is now `find` then `build` with unchanged
  results. A step that cannot be taken returns a named `RegionRefusal` and
  leaves the regions as they were.

- **Tori drilled across their equators convert whole.** `solid_from_mesh`
  places a torus's seam round the axis on a parallel no hole crosses, as
  its seam round the tube already was, so a hole across the outer equator
  is an inner wire of one toroidal face instead of leaving the tube in
  facets.

- **Open meshes keep their curved faces' free boundary.** `solid_from_mesh`
  cuts a run of boundary with one face and nothing across it where it
  turns by the crease angle, and places it on the face's surface: an exact
  parallel circle or ruling where one holds, otherwise a curve fitted
  through its vertices' feet, with a pcurve and a tolerance measured
  against the mesh. Open tubes, half cylinder sheets, domes and free-form
  sheets come back as one face each instead of facets.
  `MeshSolidReport::free_edges_fitted` counts the fitted edges.

- **Face blends end where the seat runs off a face.** `blend_faces`
  rounds two curved faces of a solid, marched between them, where the ball
  touches both over only part of its seat: a crease that runs off the
  solid, or a face holding part of the loop. The round ends where the
  ball's line of contact leaves a face, closed by a flat cap in the ball's
  section, and is cut off or fused on as the closed round is. Slanted-drum
  fills tend to the torus ring's share within 2e-5.

- **Mesh conversion fits B-spline patches.** `solid_from_mesh` rebuilds a
  smooth region no plane, canonical surface or sweep fits as one B-spline
  patch, where the region is one disk and the patch verifies: charted by a
  nearly fitting canonical surface or by the mean-value map onto a square,
  fitted with thin-plate fairing and parameter correction, and checked at
  every vertex and inside every triangle, Jacobian and normals included.
  `MeshSolidOptions::patches` (on by default) controls it, and the report
  counts patch faces and the regions refused (`patches_not_disk`,
  `patches_narrow`, `patches_unverified`). `fit::fit_surface_scattered_at`
  fits a scattered surface at the caller's own parameters and knots, with
  thin-plate fairing.

- **Rounds between curved faces with no shared edge or axis.**
  `blend_faces` and `fillet_faces` round curved faces that share no edge
  and whose surfaces share no direction or axis (a B-spline face, a
  cylinder at a slant to a plane). The ball is marched along where the
  surfaces cross, or, between separate faces whose surfaces do not cross,
  along where their offsets by the radius cross; the round is a B-spline
  band fitted through the ball's arcs, tangent within a tenth of a degree.
  On a solid the corner is cut off or filled where the ball touches both
  faces all the way round a closed seat; between separate faces the round
  spans where both faces reach, and with `trim` the faces and the round
  sew into one shell. A seat a solid's faces hold over part of its length
  only, a marched round closing on itself between separate faces, and a
  line of contact crossing its face more than once are refused by name.

- **A coincident patch's edge that weaves across the other's window is
  carried stretch by stretch.** An edge of a B-spline patch that leaves a
  coincident patch's window and comes back onto it was refused. It is
  carried into that patch's chart as one contact per stretch on the
  window, each with its own pcurve; a sheared bar across a sheared drum on
  a shared plane fuses, cuts and intersects to the exact volumes.

- **`blend_faces` blends curved faces.** Faces that meet along edges of
  the solid are rounded along those edges: exact cylinders and tori where
  the seat has a closed form, marched B-spline bands otherwise, B-spline
  faces included. Faces that share no edge are rounded with an exact
  cylinder or torus where their surfaces share a direction or an axis,
  the corner cut off or filled as the solid shows. Faces sharing no edge
  and no direction or axis are refused by name. Closed marched fillet
  bands sample each ball arc more finely, and stay tangent within a tenth
  of a degree.

- **Tangential contact traced as its own curve.** Two surfaces touching
  along a curve no closed form names (a torus resting in a free-form
  cradle) came back as fragments. The trace predicts along the null
  direction of the difference of the two surfaces' second fundamental
  forms and corrects onto the locus where they meet with parallel
  normals, landing on a patch's edges and closing across a seam: one
  curve with pcurves on both surfaces, marked tangential. A seed at an
  isolated touch, or where the surfaces come near without touching, is
  refused by name.

### Changed

- **`intersect_surfaces` answers `Same` for coincident patches.** Two
  patches with no closed form that lie on one another within the requested
  tolerance are reported as coincident for every caller, where fillets,
  offsets and sections used to march noise along them. Surface projection
  (`project_on_surface`, `project_on_surface_from`, `SurfaceSeeds`,
  `SurfaceProjection`) lives in `ogeom-geom`; the `ogeom-algo` paths still
  work.

- **Mesh conversion is 3 to 5 times faster on large meshes.** Each seam
  curve is solved once per conversion and the per-face checks (areas,
  probes, the exact volume) run in parallel, with identical results: a
  23k-triangle part converts in 9 s instead of 48 s.
  `ogeom_mesh::triangulate_with_chords` returns a shape's triangulation
  together with the edge chords its faces agreed on.

### Fixed

- **Marched sections state the error their curves have.** A marched
  intersection's tolerance is measured in space (the curve against its
  trace, each pcurve on its surface against the curve, and the surfaces'
  gap over the sine of their angle where they cross shallowly) instead of
  the fit's residual times the surface's stretch, which on a wide drum
  stated hundreds of times the real error and widened the boolean's
  junctions as far. A branch whose fit misses by a hundred times its
  tolerance is fitted in pieces, and a section reaches an edge within its
  own reach plus the edge's tolerance. Two more drills in the stress run
  succeed.

- **Booleans with a ball whose poles lie on the other solid's wall.** A
  drill whose wall passes through both poles of a ball cut nothing away or
  everything (an empty or whole result `check` called valid), or took the
  common a little short. A whole ball meeting a curved wall at a pole is
  charted about another axis of its frame, clear of the wall, before the
  boolean; cut and common add up to the ball and match an independent
  integral within 1e-6 over drills in six directions and three radii.
  Empty results are checked too: a cut leaves nothing only where no part of
  the solid lies outside the tool, a common only where the solids share no
  volume, and a fuse never; otherwise each is refused by name.

- **Mesh conversion tells a finely drawn curve from scatter.** The
  estimate behind the default coplanar distance weighs each nearly flat
  edge's turn against its neighbours': a steady turn from one edge to the
  next is a curve, not noise. Slabs drafted by fractions of a degree keep
  each row a face of its own at the default distance (pads into them fuse,
  cut and are taken in common at every row count tried), and a free-form
  sheet with no flat stretch comes back one fitted patch instead of
  thousands of facets.

- **Exact area and volume on a face holed by a closed spline.** A face's
  loops are told apart by their areas, integrated span by span; one rule
  over a whole closed spline could read a hole's direction wrong and add
  the hole to the face.

- **Edge rounds stop where a drum's foot runs off the block.**
  `fillet_edges` on an open elliptic crease (a leaning drum whose foot runs
  off the block) rounded the whole ellipse and left material past the
  block's side; the round now covers the arc only, or is refused by name
  where the boolean cannot merge its end with the face the crease runs out
  on. Fillets and chamfers on a circular foot that runs off the block took
  the whole turn, read the wrong side of the wall and removed material;
  they now round only the arc, with exact torus or cone fills.

- **Pads fuse into slabs whose rounded corners lean in by less than
  their chords sag.** Where a pad wall crosses an edge between two facets
  lying almost in the wall, the two facets' sections reach the edge up to
  a few ten-thousandths of a millimetre apart; they now join as one point,
  and a stub of section between two such crossings no longer leaves a face
  piece with no inside. A slab drafted 0.001 fuses, cuts and is taken in
  common at every row count tried, while its rows stay faces of their own.

- **A boolean cuts periodic faces whose outlines sit turns away from
  their chart, and fits sharp-tipped section loops in pieces.** A face on
  a periodic surface whose outline was stated one or two turns outside the
  chart's range read every section piece as outside it, so the section was
  dropped and the cut did not close; its pcurves are moved home by whole
  periods when it is gathered. A closed section loop that turns sharply
  where two drums all but touch, fitted whole, missed its trace by
  millimetres and swelled the junctions it fed; it is fitted in halves.
  Two oblique drills through a NIST part that were refused give valid
  solids whose cut and common add up to the part.

- **A sphere face bounded by a great circle through its poles is
  healed into half meridians.** Its stored pcurve was a fit that strayed
  far from the meridian, so a boolean's section ended away from where the
  boundary was paved and a drill across the hemisphere was refused.
  Healing splits such a circle at the poles into exact half meridians,
  joins pieces meeting on a pole's row by a pole edge, and rebuilds the
  neighbouring face on the pieces; the repair is refused where the rebuilt
  face's area would change. The part's volume then integrates exactly.

## [0.7.0] - 2026-10-02

A minor release, since the API changes: `Canonical` is no longer `Copy`,
`MeshSolid` and `BlendContact` gain fields. Its larger part is surface
modelling on sheets: sewing across a gap and a prism's lids, trimming and
splitting sheets with booleans, splitting a face along curves, growing a
face past an edge, offsetting and thickening sheets, lofts, ruled surfaces
and sweeps as sheets with guide curves, blend surfaces and curves, N-sided
fills on their own edges, fillets on sheets and between separate faces,
curvature continuity in blend analysis, and sheet export to STEP and IGES
(#104 to #115). The mesh converter rebuilds surfaces of revolution all the
way round, fillets and corner balls from the faces they blend, and every
converted part it is measured on now tessellates closed; a boolean keeps a
piece reaching across a periodic face's seam as one face and classifies a
piece touching the other solid only along lines.

### Added

- **Fillets on sheets and between separate faces.** `fillet_sheet_edges`
  rounds edges where two faces of a shell meet, and `fillet_faces` rounds
  the corner between two faces of separate shapes, the ball on the side
  their normals point to, optionally trimming both back into one shell
  with the round. The round is an exact cylinder where the supports share
  a direction, an exact torus where they share an axis (a plane square to
  it, cylinders, cones, spheres, tori), and on a sheet a marched B-spline
  band elsewhere. Each face is rebuilt along its line of contact and
  shares that edge with the round. A B-spline pair given to
  `fillet_faces`, and a closed sheet edge with no closed-form round, are
  refused by name.

- **Sheets offset and thicken.** `offset_sheet` moves a face or a shell
  along its normals, and `make_thick_sheet` thickens a face or an open
  shell into a solid, on one side or both, closing it with side faces
  along its free edges. Analytic faces move to their exact parallel;
  B-spline and other faces are fitted within the approximation tolerance,
  the measured deviation held in the face's tolerance. Trims carry over
  unchanged. Creases between faces, offsets past a radius of curvature,
  edges shared by more than two faces and layers that run into each other
  are refused by name. `fit_surface_grid_at` fits a tensor-product grid at
  given parameters, keeping the caller's chart.

- **N-sided fills on their own edges.** `make_filling_n` fills a hole
  bounded by any number of edges with one face bounded by those same
  edges, with pcurves, so `sew` joins it to its neighbours. Each side
  meets its support face G0, G1 or G2, and the surface passes through
  given interior vertices and edges. The gap, angle and curvature step
  reached are measured per side and returned in `Filled`. A fit that
  cannot meet the tolerance, and a hole that is not a height field over
  its boundary's plane, are refused by name.

- **Blend surfaces and blend curves bridge two edges.**
  `make_blend_surface` bridges two edges of different faces with a
  B-spline face that shares both edges and meets each face C0, G1 or G2:
  exactly against planes, and within a measured budget (1e-5 rad, a
  curvature step of 1e-4 over the gap's width) against curved faces.
  `make_blend_curve` bridges the ends of two edges G1 or G2, `End` naming
  an edge's start or end as traversed. Circle and ellipse edges, closed
  and seam edges, and C1 or C2 asked across unrelated surfaces are refused
  by name.

- **`split_face` cuts a face along curves on it or projected onto it.**
  Curves lie on the face, or are dropped onto it along its normals or
  along a direction (`Projection`); curves meeting end to end are
  followed as one, and a curve crossing a periodic face's seam is followed
  across it. Each stretch from boundary to boundary divides the face, and
  a closed curve inside it cuts out the region it encloses. The pieces
  share the new edges, neighbouring faces keep sharing the boundary edges
  a cut ends on, and history maps the face to its pieces. A curve that
  ends inside the face, misses it, or runs along its boundary is refused
  with a reason.

- **Booleans trim and split sheets.** `common` and `cut` accept a face, a
  shell or a compound of those against a solid or a half space, and return
  what of the sheet lies inside or outside it, sewn into a shell or a
  compound of shells. `split_sheet` splits a sheet where another sheet
  crosses it and keeps every piece. Pieces keep their surfaces, and
  history maps each to the face it came from. A sheet as the tool of
  `cut`, two sheets given to `common`, a sheet given to `fuse`, and a
  piece lying on the tool's boundary are refused by name.

- **Sheets from lofts, ruled surfaces and sweeps.** `make_ruled`,
  `make_loft_surface`, `make_sweep_surface` and `make_sweep_two_rails`
  build sheets from edges or wires, open or closed, planar or not. Lofts
  and ruled surfaces pass through their sections exactly (exact rational
  B-splines, degree one across when ruled, a plane between coplanar
  segments, optionally closed back to the first section). Sweeps follow a
  frame law or two rails, held to ten confusions of the swept profile and
  measured. A sheet has a face per section edge, and bounding section
  edges are the caller's own, so the sheet sews to them. A loft follows
  guide curves as a Gordon surface through single-edge sections and their
  guides: exact on the sections and within ten confusions on each guide,
  measured; a guide missing a section is refused by name.

- **Faces grow past a boundary edge with `extend_face`.** A face is
  rebuilt with the edge moved `length` out along the surface, square to
  it, and the face's other edges stay shared with its neighbours. Planes
  grow across any straight edge; cylinders, cones, spheres and tori grow
  along a parameter line on their own surface; B-spline faces grow past
  the side of their patch, naturally or linearly at G1 or G2
  (`Extension`). History maps the old face to the new one, and the cases
  it does not handle are refused by name.

- **Blend analysis measures curvature continuity, and faces give
  curvature samples.** `BlendContact` carries `curvature_error`, the
  largest step in normal curvature square to the shared edge over the
  stations; an edge that cannot be measured reports it infinite.
  `face_curvature_samples` returns the principal curvatures on a grid
  inside a face's trim, placed where the face stands and signed against
  its own normal, for curvature maps and zebra stripes.
  `SurfaceCurvature::normal_curvature` gives the curvature along any
  tangent.

- **STEP and IGES export write sheet bodies.** `write_step` writes shells
  no solid owns and faces no shell owns as `SHELL_BASED_SURFACE_MODEL`s in
  a `MANIFOLD_SURFACE_SHAPE_REPRESENTATION`, and a part holding solids and
  sheets as both representations in one part. `write_iges` writes sheet
  faces as trimmed surfaces (144) beside the manifold solid B-reps it
  writes for solids, keeping a reversed face's side. Both writers refuse
  wireframe (a free wire, edge or vertex) and empty parts by name instead
  of dropping them.

- **Sewing across a gap.** `sew_within` sews faces whose edges meet up
  to a given distance apart, as surfaces built separately (an imported
  sheet, a fitted fill beside an exact extrusion) do: two edges whose ends
  and middle lie within the gap are one, and the edge kept and its
  vertices widen their tolerances to reach the edge it replaces. An edge
  that runs along two or more shorter edges of other faces is first split
  where their vertices meet it, so each piece has a twin. Two unit squares
  a twentieth of a millimetre apart sew across a tenth into one shell and
  stay two across a hundredth; a 10 mm edge against two 5 mm edges sews
  into one shell.

- **A boolean takes a compound of disjoint solids as an argument.** A
  grid of drums cut from a plate at once is one boolean, not one per drum:
  four hundred bores take 0.2 s where four hundred cuts took over a
  minute.

- **Refine around the faces an operation made.**
  `unify_same_domain_around` merges the given faces with the faces on
  their carriers they reach across shared edges, and leaves every other
  face of the solid as it is, split or not: an application refining after
  a feature unifies what the feature made without undoing splits the part
  already held.

- **Mesh conversion rebuilds extrusions and surfaces of revolution.**
  With `MeshSolidOptions::sweeps`, a smooth region none of the canonical
  surfaces fits is tried as an extrusion of a fitted profile (its normals
  all square to one direction) and as a surface of revolution (its
  normals all meeting one axis), each held to the coplanar distance at
  every vertex. The axis read from the normals is brought close by a
  search on the samples' radii against their heights, then settled by
  Gauss-Newton on the samples themselves: each round fits the profile and
  solves for the axis's offsets and tilts with a correction to the
  profile, which separates a tilt from the profile's shape over a whole
  turn as over part of one. A recognized band the region runs into
  smoothly is offered with it, so a stretch of a wavy profile that passes
  for a cone joins the one sweep. A wall extruded from a wavy spline comes
  back one extrusion and five planes, not nineteen facets; a wavy profile
  turned half round, one surface of revolution, and turned all the way
  round, one surface of revolution closed on itself between two flat
  ends. On by default.

- **Mesh conversion builds fillets and corner balls from their
  neighbours.** A torus between a plane and a cylinder or cone on an axis
  square to it is put on that axis, its tube's centre a radius off both;
  rounds between planes that meet at corner balls take one radius, and
  each ball is centred a radius off the three planes its rounds run
  between, where their axes meet. Only the radius is fitted, and each
  surface is kept only where it holds every vertex of its region. From a
  single-precision mesh of a rounded block the corner balls now meet their
  fillets tangentially to rounding, where the fitted ones stood off their
  axes by a few tenths of a micron.

### Fixed

- **A boolean piece touching the other solid along lines or points is
  classified, not refused.** A face piece whose every probe landed on the
  other solid's boundary (a cylinder inscribed in a box touches its walls
  at the quarter turns the probes use) was refused for want of a
  coincident partner. It is asked again at points off the contact, which
  a contact along lines or points cannot cover. A hole touching its outer
  boundary at isolated points (a disc inscribed in a square) is placed by
  a vote of chord midpoints, not by a point that may be the one it
  touches at.

- **A boolean keeps a piece reaching across a periodic face's seam as
  one face.** The boolean splits each face in its own chart, whose edge
  is the seam, so a piece reaching across the seam came back as two faces
  sewn along it: a tube or drum halved square to its seam gave three
  faces. Pieces of a cylinder, cone, sphere, torus or full revolution
  that meet across its seam are joined into one face, on the surface
  turned so its seam falls where the face does not reach, in solid
  booleans, sheet booleans and `split_sheet` alike. `History` gains
  `without_repeated_images`.

- **A prism's far end edges sew.** `sew` told edges apart by their nodes,
  and a prism's far end edges are its profile's own edge nodes carried
  along it by a location: each was read where the profile stands, so a
  lid at the far end joined nothing. Faces that hold an edge node placed
  twice are baked first, every edge its own node where it stands, and an
  edge's curve is read through its own location as well as the edge's.

- **A converted solid that would mesh open gives up the curved faces that
  open it.** A face whose trim folds within the tolerance its own seams
  claim passes every check on the solid and still meshes over itself.
  Where the solid's tessellation does not close, each face is meshed on
  its own, and the curved faces using a mesh edge other than twice (or,
  where only planes do, with an edge drawn within the reach of it) fall
  back to facets. Every converted part of the corpus and the NIST set now
  tessellates closed: ftc_07 gives up 20 curved faces of 143, Body28 3,
  ctc_02 2, handle-pickup 1. A part that needs it converts slower (ctc_02
  from 27 s to 50 s); one that meshes closed pays one tessellation.

- **A converted torus is no longer built with its seam through a hole.**
  A hole's outline was unwrapped round the torus's axis but not round its
  tube, so a hole across the outer equator, where the tube's angle starts,
  read as clear of the seam there. The whole torus was built with the hole
  as an inner wire and its seam edge running through the hole: a face
  `check` passed and no tessellation could close. Such a ring's tube now
  stays facets; its bore still comes back a cylinder.

- **Faces meet at a vertex their edges end off.** An edge's curve may end
  off its vertex by up to the vertex's tolerance, and the next edge's off
  it on another side. Each face round the vertex closed its own gap
  between two of the ends, and between them the faces left a hole; on a
  curved face the boundary also ran from the curve's end to the vertex and
  back in its chart, a fold its triangles covered twice. Each edge is now
  drawn from its vertices, in space and in every face's chart, so the
  faces round a vertex meet at one point. Converted parts whose faces,
  meshed one by one, left edges unmatched: 77777_1 from 64 to none, Body28
  from 122 to 31, ctc_02 from 25 to 6, ctc_01 and ctc_03 to none. Every
  imported STEP solid in the corpus closes as before.

- **Facets left between curved faces join them.** Where curved regions
  meet, or a fit stops a row short, a mesh leaves planar facets of a
  triangle or two that no surface claims, bounded by seams as wide as
  themselves, whose trims fold. A facet whose corners all lie within the
  seams' reach of a curved neighbour joins it; a region that cannot be
  built with the facets it took gives them back before it is faceted, and
  a conversion that leaves a face turned into the material is made again
  without them. A coarse plate with filleted rims comes back its own 18
  faces instead of 56, the truth bench's faces over the exact count fall
  from a mean log ratio of 0.27 to 0.17, and large real exports lose up
  to a third of their faces with every curved surface kept.

- **A needle on a flat face is part of it.** A triangle no higher than
  the coplanar distance (three corners all but on a line) has no plane of
  its own at that distance, but its rounded normal could lean far enough
  that it became a face of almost no area, turned whichever way its
  rounding pointed, which `check` could then find facing into the
  material. Such a triangle joins a neighbouring plane by its corners'
  distance alone.

- **Mesh conversion raises a coplanar distance the mesh cannot meet.**
  With recognition, a distance asked below two and a half times the
  mesh's measured scatter (or twice its quantum) left regions in
  fragments whose seams opened: a single-precision rounded block asked to
  a tenth of a micron tessellated open while `check` passed. It is raised
  to that floor, `MeshSolidReport::coplanar_distance_raised` says so, and
  the new `MeshSolid::coplanar_distance` is the distance the faces were
  built to.

- **A converted face whose boundary folds over itself is withdrawn.** Each
  built face's boundary is drawn in its surface's chart, and two of its
  chords that cross are a fold. A fold deeper than the tolerance of the
  two edges that cross, outside every corner's tolerance, on edges of the
  face other than a seam, and still there when the boundary is drawn a
  hundred times finer, names its two seams, and the curved faces across
  them fall back to facets. `check` passes such a face; its triangles
  cover part of it twice.
- **A seam threaded along a tangency no longer hooks past its corner.**
  Where a fillet meets a face tangentially the mesh's boundary between
  them wanders, and can step past its corner and back; the curve threaded
  through it ran on past the corner and came back, folding the trims of
  both faces. Only the points that close in on the seam's end are
  threaded, and where a cubic through them still overshoots, the polyline
  through them is.

- **A point on a face of any repeating surface classifies on the right
  side.** Points were folded by a period toward a face's outline only on
  cylinders, cones, spheres and tori; on a whole turn of a revolution or a
  periodic patch, an outline kept at negative angles left the point a
  period away and the solid classified inside out. The surface now says
  whether and how it repeats.

- **Crossings along an edge each known to within its doubt meet where
  both doubts reach.** A pad pushed into a converted slab whose corner's
  rows hold different numbers of facets refused: the facets' diagonal
  edges run all but in the pad's walls, and the walls' sections on the
  facets either side of one reached it tens of microns apart, each within
  its own doubt but not within the larger one alone. Two crossings are now
  one point where they stand no further apart than their doubts together.

- **A coarse mesh comes back the same whatever order its triangles are
  listed in.** A fillet two or three facets wide, first sampled together
  with the corners beside it, fitted nothing and was left faceted, or not,
  by the luck of the order. Its triangles are now sampled again once the
  corners are claimed, and a first sample no longer reaches across a flat
  face into the fillets on its far side. A sphere or a torus fitted to a
  handful of noisy vertices is no longer taken. Several converted parts
  come back with far fewer faces; one of 6268 triangles, 1615 faces before,
  comes back in 13.

- **A converted part tessellates closed where a single facet meets a
  fillet.** A facet left beside a fillet, thinner than the curve its plane
  meets the fillet along bulges, had that curve as its seam: the facet's
  trim folded back across itself, its triangles covered the fold twice,
  and the whole part's tessellation opened there while every check passed.
  Such a seam is now threaded straight through the facet's corners, its
  tolerance how far it stands off the fillet, and the fillet keeps its
  surface; where the straight seam strays too far, the fillet is built
  faceted. A sliver between two curved faces whose seams cross in its
  plane, and any facet whose mesh does not meet its neighbours', is
  threaded the same way. A threaded build that leaves a face turned into
  the material or collapsed onto a line is set aside for the unthreaded
  one. Three converted corpus parts and the rough rounded box now
  tessellate closed.

- **A polynomial's roots are found whatever their pattern.** The real
  roots of a quartic or higher came from the companion matrix's
  eigenvalues, whose iteration had no limit; a quartic whose roots pair
  off by sign (a line through a torus's middle) defeated its shifts and
  the call never returned. The iteration is now capped, and where it does
  not settle the roots are found by Durand-Kerner instead.

- **Halves of a cone sharing their seam ruling mesh closed.** The ruling
  two halves on one surface share has one image in the chart, at the seam,
  and the half running from a full turn back to half a turn needs it a
  turn over. Walked into the apex the two columns are equally near, the
  ring wound onto the other half, and both halves meshed over one side: a
  NIST part read from AP242 tessellated open. The winding is now undone
  from the tie where the walk left it.

- **A converted solid no longer comes back with a facet facing into its
  material.** On a coarse mesh a fitted face could run past the facet it
  should have ended on and enclose a sliver outside it; the solid checked
  sound but for that facet's orientation. Once a build would be accepted,
  the faces round every recognized face are probed as the validity check
  probes orientation, and a recognized face that overlaps another is built
  faceted. The STEP corpus, meshed and converted back, is now a test.

- **A face moved or offset clear of a neighbour it was tangent to is
  refused.** Moving the side of a block rounded all round outward by half
  a millimetre left the face's surface there and its edges where they had
  been, on fillets it no longer met, and the solid came back valid with
  less volume than before. Every corner the edit re-solves is now held to
  lie on each surface it sits on, and the edit is refused by name where it
  does not.

### Changed

- **`recognize::Canonical` is no longer `Copy`, and has a `Swept`
  variant** holding a surface of revolution or an extrusion
  ([`SweptShape`]). Code matching on it needs the new arm; code copying
  it, a `clone`.

- **`tight_bounds` is hundreds of times faster on ruled and round
  parts.** Every face was meshed and searched for an interior extreme, and
  every edge searched numerically. A face on a plane, drum, cone or linear
  extrusion reaches no further than its boundary, a straight edge leads at
  its ends, and a circle at an angle its frame names: those are taken as
  they are. A plate with four hundred bores bounds in 1 ms, not 311 ms, to
  the same box.
- **A small cut into a large solid costs what it touches.** A face the
  tool leaves alone (no section on it, no edge of it split, no coincident
  partner) is neither split nor arranged: its rings are its own wires and
  its state is read once, at the middle of one of its edges. The other
  solid's trims are drawn only where a point or a ray comes near them, a
  crossed face's outline is drawn once, circles are sampled to their sag
  in closed form, outline ends are welded among their neighbours, a
  piece's scanlines come from a heap as the search asks, and an edge
  rebuilt whole is built once for both its faces. Sewing looks up eight
  grid cells per point, not sixty-four, compares a loose edge only with
  edges starting within its reach, and finds a middle's foot on a line or
  circle in closed form. A pocket into a plate with four hundred bores
  takes about 50 ms, from 160.
- **A line meets a cone or a torus in closed form.** Both took the
  general seeded path, a Newton polish per seed; a cone's crossings are
  now a quadratic's roots and a torus's a quartic's, tangencies included,
  each polished once onto the surface, the quartic solved about the
  line's nearest approach to the torus's centre in units of its size, so
  a torus far down the line is not lost to rounding. Every ray a
  classification casts
  through a part with tori is cheaper: converting a part with three tori
  takes 2.5 s, not 27.

- **A chain of edges rounds in a few booleans, not one per edge.** The
  chain is taken in rounds, no two edges of a round sharing a vertex
  (pieces of one circle count as one); each edge's blend is built against
  what the earlier rounds left, and a round's blends are applied together,
  one cut and one fuse. The checks before the blends mesh each face once
  and pass over edge pairs whose boxes stand apart. A hundred bore rims
  round in 1.4 s, from 10, and a box's twelve edges in half the time.

## [0.6.1] - 2026-10-01

A patch release: `cargo semver-checks` against 0.6.0 finds no breaking
change, and every crate packages and verifies in a publish dry run. Most of
it answers an adversarial review of the kernel: mirrors and negative scales
transform curves, surfaces, primitives and booleans the right way round;
spheres come back whole through STEP and IGES; unclamped and periodic
splines reverse, rescale and split unchanged; offsets and fillets refuse what
they cannot build instead of returning it inside out. The boolean settles
near-coincident faces by one rule, welding a sliver thinner than the weld
distance and keeping anything thicker as real geometry, and pads on a part
converted face for facet now fuse back into it. A pipe tee's junction loop
and a slab's full round fillet, and a thread groove cuts the same at any
start angle. A history now says which faces came through a boolean
unchanged (`History::copy_of`), and a boolean on a face with hundreds of
holes is several times faster.

Known limits, each with its diagnosis in `docs/PLAN.md` or
`docs/REVIEW.md`: a converted rounded block's whole-shape tessellation can
stay open along its corner spheres while `check` passes, and a pad into a
slab whose faceted walls lean in by less than their chords sag refuses at
some facet counts.

### Fixed

- **Pads on a part converted face for facet.** A pad pushed down from a
  top face of a converted part, its walls standing on facets a hundredth of
  a radian off them, refused ("the kept pieces did not close into a
  shell"). Two planes sharing an edge only to their points' rounding were
  crossed on a line the solve left a sliver under the edge; they now cross
  on the edge itself. And where a facet row's edge runs all but in a wall's
  plane, the sections on the facets either side reached it microns apart
  and split it twice; they are now one junction.

- **A slab rounded at half its thickness.** Rounding every edge of a box
  at half its thinnest side refused ("an earlier blend in the chain
  consumed this edge"): the balls of the two corners at each short edge's
  ends meet across it and leave nothing of it. An edge of the chain whose
  ends are both corners the chain rounded, and which they consumed, is now
  rounded by them; the result is the full round, its volume Steiner's.

- **A closed crease asked for as a chain of its arcs.** `fillet_edges` on
  every arc of a junction loop refused after the first ("an earlier blend
  in the chain consumed this edge"): the first arc's blend runs the whole
  loop. An arc of the chain lying between the same two faces as one
  already rounded is now done, and the history records it consumed.

- **A pipe tee's junction loop would not round.** Filleting the loop where
  a branch drum meets a main drum refused ("the kept pieces did not close
  into a shell"). The blend's leg on the branch is a band between the loop
  and the ball's contact ring, and the two rings' chart images started on
  one column a period apart, so the band's seam wound a full turn round the
  branch as a helix. The second ring's image is now moved into the first's
  period and the seam runs straight along the column.

- **Drums the weld distance apart in axis.** Two drums whose axes stood
  1e-5 to about 1.07e-5 apart refused every operation. The rim arc across
  the thin crescent one leaves on the other's cap was asked whether it lies
  inside the cap only at its middle, where it lies inside by less than the
  boundary's own doubt; it was dropped and the drum was not split. A piece
  is now asked at points either side of its middle as well.

- **A thread groove cut into a rod depended on its start angle.** A
  helical sweep cut into a body of revolution failed at some angles about
  the axis ("edge N ends where edge N+1 does not begin", or a shell that
  did not close) and passed at others. A marched section stopped a
  fraction of a step short of the edge of the patch it left, so it missed
  the next patch's section by more than either curve's tolerance. The
  marcher now ends such a section on the edge itself; a sphere's pole and
  other bounds that collapse to a point are not edges and are left alone.

- **Blocks tilted by a hair, and slivers at the weld distance.** Two
  blocks stacked and hinged, or standing in one another turned, by 1e-5 to
  1e-4 rad refused. Two planes all but parallel now cross on their true
  line (measured from a point on the first, not the world origin, whose
  offsets cancelled to a third of a micron); an edge of either lying on the
  other's plane splits it as a coincident pair's edge would; and a target
  edge splits where a contact ends on it. The welding retry takes a second,
  wider step, held to the same tolerance limit, and drums whose sliver
  lies at the weld distance itself are sectioned rather than lost between
  the two readings.

- **Slivers under the weld distance are welded.** Drums a few microns apart
  in radius, axis or end, and a box tilted a few microradians on another,
  refused. A boolean that refuses now runs once more at a confusion
  distance ten times wider, where such a sliver is one the boolean welds,
  and keeps that result only if nothing in it states a tolerance past the
  weld distance, a hundred confusion distances, so no thicker sliver is
  lost. Drums 3e-5 apart in axis, whose crossing lines were not sought
  closer than 1e-4, now meet in them.

- **Drums a hundredth of a micron to a tenth of a millimetre apart.** Two
  drums 2e-5 to 1e-4 apart in radius, or 1e-4 apart in axis, refused or
  (offset in axis) came out with a misbuilt face. Where a planar face has a
  curved edge, a section's place inside it and a piece's probe are now read
  against the exact boundary rather than its chords, a point within the
  weld distance of the boundary counting as on it; two circles in one plane
  cross where their closed form says, not along a sampled run millimetres
  long; a closed section no longer has a crossing pulled to its own
  arbitrary start; and an exact section ending at a junction is carried
  onto it rather than widening the face's weld to the junction's span.
  Volumes match the closed forms.

- **A cavity touching its wall is the solid's void.** Cutting a ball from
  a box it touches from inside left the box as one solid and the ball as an
  inside-out solid of its own beside it, where the ball reached a tenth of
  a micron past the wall: shells were nested only when one's box held the
  other's exactly. The container's box now takes the touch tolerance.

- **Faces within tolerance of the other solid's.** A drill whose wall
  pokes a tenth of a micron out of a box side, a ball resting a tenth of a
  micron into a box top, and a box stacked on another and tilted a
  microradian all refused: a piece lay on the other solid's boundary with
  no coincident partner face. Where that refusal stands after the nested
  fallback, the boolean now runs once more, settling each such piece by
  asking just off it on both sides: outside or inside on both, it is
  that; inside only on its own material side, the other's face runs along
  it backed the same way.

- **Small corner rounds and an L-bracket's every edge.** Rounding a
  block's three edges at a corner refused at any radius up to 0.1: the
  ball patch's rim leaves each sharp edge tangentially, and an arc stated
  as a trim of its circle was held to the loose overlap width meant for
  fitted curves, so the boolean read the rim as running along the edge for
  a stretch and split it where the band beside it was whole. A trimmed
  exact curve is now as exact as its basis. Every edge of an L-bracket
  failed in either order: the corner tool took the step's re-entrant
  corner, where a concave edge meets two convex ones, for a convex corner,
  because the planes through it taken whole look like one. It now probes
  beside each corner edge and leaves a re-entrant corner to the blends.

- **A STEP void reads as a cavity whichever way it is written (#102, #103).** A void
  shell whose flags turn it to add material (written inside out and then
  marked reversed as well, as some exporters do) read inside out: the
  cavity counted as extra material and `check` reported its faces broken.
  Each void is now oriented by its geometry as it is read, with a warning
  naming it, and keeps its faces and placement.

- **Drums a thousandth apart.** Two drums 1e-3 apart in radius or axis
  refused: the wider one's cap leaves an annulus or crescent that thin, and
  its outline's chords bowed by more than its width, so the point the
  piece was classified by fell outside it. A curved strand on a plane is
  now sampled until its chords stand within 1e-4 of it. Two drums 1e-4
  apart in axis, whose union came out short by 0.03, now refuse instead.

- **A cut by an identical solid, and a cap that all but meets a wall.**
  Cutting a solid by an identical copy described another way (two spheres
  on turned frames) failed in the arrangement; a solid lying within its
  tool is now cut away entirely. A drum a hair wider than the box it
  stands in refused at most turns of its seam: a planar face's bound was
  its rim's sampled points, and a disc's rim bows out between them by more
  than the drum's cap pokes past the box wall, so the pair was never
  intersected. A plane's bound now carries its rim's sampling slack, as a
  ruled face's does.

- **Three intersector edge cases answer.** A plane tangent to a cone
  along a ruling read `Apart`: a plane through a cone's apex now holds its
  rulings in closed form (one tangent, two crossing, or the apex alone),
  each stated over the cone's window. One circle written twice in the plane
  gave both whole domains as its overlap, whatever the phase and winding;
  `on_b` now names where the second stands at the first's ends, as in
  space. Collinear segments meeting end to end, in the plane or in space,
  returned nothing; they share their end point.

- **A periodic pattern traces its faces and may have gaps.**
  `make_periodic` recorded history only for the whole solid; every face
  and edge of the cell now traces to its own image and generates its
  copies' images through each fuse. A pattern whose copies stand apart
  failed at the second copy, the fuse refusing the compound the first one
  left; the copies are now fused only into the pieces they join.

- **Every outline is drawn once.** HLR `project` read its silhouettes off
  a mesh whose faces were meshed apart, so every face border counted as an
  outline and every outline edge was drawn a second time over the model
  edge. The mesh is welded first, and a silhouette along a drawn edge is
  left to the edge.

- **A bare or placed profile tapers.** `make_prism_tapered` refused a
  face straight from `make_face`, whose edges carry no pcurves yet, and a
  placed profile, whose material side it probed in the wrong frame. The
  bare profile gains its exact pcurves, and a placed one is swept where it
  was built and the prism placed where it stands.

- **A straight run of corners is no obstacle to a face.**
  `make_polyhedron` took a face's normal from its first three corners and
  refused a ring starting on a straight run (and would have turned a
  concave ring starting on a reflex corner inside out); the whole ring's
  winding decides now. `medial_axis` refused a polygon with a corner on a
  straight side; such a corner is dropped, and the axis is the polygon's
  without it.

- **Touching circles have their third tangent line.**
  `lines_tangent_to_two_circles` dropped a pair of tangents that had
  closed into one, so circles touching from outside got two lines instead
  of three and circles touching from inside none instead of one.

- **Polynomial roots keep their double roots.** A tangency's double root
  sits where the discriminant is zero only up to rounding, and the quadratic
  and cubic solvers dropped it when rounding left the discriminant a hair
  below or above; above the cubic, the companion matrix split it into a
  complex pair and it was dropped there too. All three now find it. The doc
  no longer claims a closed-form quartic.

- **Degree elevation keeps the continuity at every knot.**
  `elevate_degree` left each interior knot at full multiplicity, a C0
  joint where the curve was C2, against its doc; the knots it introduced
  now come out exactly, each interior knot one higher than it was.

- **Conics and helices answer every derivative order.** A circle's,
  ellipse's, hyperbola's and helix's derivatives past the second (the
  helix's past the third) came back zero, as did a 2D circle's and
  ellipse's; they are now the closed forms at any order. A parabola's and a
  line's higher orders are zero, as they should be.

- **A near half turn and a small chart keep their digits.**
  `Quaternion::between` refused directions within some 5e-4 rad of a half
  turn (a length tolerance held against a cosine) and lost digits short of
  that; it now reads the angle off the cross product and refuses only within
  the angular tolerance. `SurfacePoint::normal` held the tangents' cross
  product, an area, against a length and refused a healthy chart a tenth of a
  millimetre across.

- **A large box measures at once and exactly.** The exact surface
  integral cut every chart rectangle into quarter-turn panels in both
  parameters, so a plane's millimetre parameters gave a metre cube some
  four hundred thousand panels: seconds per measure, longer at three metres,
  and drifting in the tenth digit. A length parameter now takes one panel
  and an angular one at most sixty-four.

- **`make_wedge` documents its zero top extents.** The doc said a zero
  top extent is refused; it builds a ridge or a point top, five faces either
  way, and the doc now says so.

- **A mesh of positions alone welds.** `Triangulation::welded` indexed
  the normals and parameters of every vertex and panicked on a mesh read
  without them; the optional arrays now follow only where they are full.

- **A converted torus measures exactly and divides into spans.** A torus
  stated as a spline far from the origin had its closing seam sampled back
  onto the chart's start, so its seam sides were fitted with ringing and the
  volume came out off by a few tenths of a percent; a spline seam that is a
  chart line now takes the exact measure, also when it is cut into pieces.
  `divide_by_continuity` walked a torus's ring of seams with its first seam
  on the wrong chart side and never cut it, so `to_bezier` left one
  multi-span face; it now cuts into sixteen, each patch held to its own span.

- **A drum shorter than a unit divides.** `divide_by_angle` and
  `divide_by_area` read a straight cut line's direction off the surface at
  parameters zero and one, past the window of a drum or cone less than a
  unit high, and refused. The line is read inside its own span.

- **A mesh reports the chord it meets.** On a surface curved both ways a
  grid cell's middle can stand off the surface by up to twice the chord,
  while `deflection_met` said the chord was met. Each triangle's middle is
  now measured against the surface, and the flag is false where any
  stands off by more than the chord.

- **A scaled body meshes to the chord where it stands.** A face is drawn
  in its own frame and placed after, so a body scaled by its placement
  was meshed to the caller's chord before the scale stretched it: ten
  times the chord at a tenfold scale, and mass properties two percent off.
  The chord is now taken where the face stands.

- **Tight bounds reach a sphere's and a torus's extremes.** In a general
  frame `tight_bounds` could come out up to 45 microns small on a sphere or
  a torus, its descent stopping short of the extreme. Their extremes are
  closed form and are now taken wherever they lie inside the face.

- **A placed face's medial axis moves with it.** `medial_axis` read a
  moved face's boundary in its surface's unplaced frame and lifted the
  axis back through it, giving an axis neither where the face stands nor
  where it was written. The frame is now placed with the face.

- **Lengths and points along a curve turn its kinks.** `curve_length`
  could fail across a polyline spline's corner, and `parameter_at_length`
  (with `points_by_count` and `points_by_spacing` on it) took the failure
  for a zero length and returned wrong points silently. Lengths are
  integrated span by span between a spline's knots, and a failure is
  reported.

- **A periodic spline is closed.** `is_closed` on a B-spline compared its
  end control points, which are its ends only for clamped knots: a
  periodic spline, and a surface extruded from one, read as open. A
  periodic spline is closed, and any other is closed where its ends meet.

- **A conical helix measures its own length.** `HelixCurve::arc_length`
  ignored the taper and reported a cylindrical helix's length; a conical
  helix's speed grows with its radius, and its length now integrates it
  exactly.

- **Scaled curves name the images of their points.** Under a scaling a
  planar trim of anything but a bare line (a trim of a trim, of an offset)
  kept its old window, a parabola its old domain, and a curve on a surface
  its old chart curve while the surface's length directions rescaled; each
  then named the wrong points. All three now follow the rescaled
  parameter, and the whole-shape conversion carries edge ranges the same
  way.

- **Tangent circles keep their digits far from the origin.**
  `circles_tangent_to_three` and `circles_of_radius_tangent_to_two`
  solved rows holding each target's squared distance from the origin, so
  a small circle ten metres out came back off by a percent of its radius.
  Both now solve about the targets' own middle.

- **A nearly circular ellipse is a circle, not NaN.** An ellipse whose
  minor radius exceeded its major by less than the confusion tolerance
  was accepted and reported a NaN focal distance and eccentricity. Such a
  minor radius is stored as the major.

- **Two circles' bisector bends round the smaller one.** `bisector` of two
  circles of unequal radii put the hyperbola's equidistant branch round
  the second circle whichever was smaller, the wrong branch when the
  first was. Its frame now points toward the smaller circle.

- **A cone narrowing up its axis measures its distance.** A negative half
  angle (a cone narrowing in `+z`, which tapered sweeps and chamfers
  build) is accepted, but `Cone::distance_to` measured such a cone's own
  points as off it. Its nappes now lean out from the axis whichever way
  the cone widens.

- **Cone and torus inversion past the apex and the axis.**
  `cone_parameters` returned a point's own angle and height, the chart
  point only on the cone's own nappe and only for a point already on it;
  `torus_parameters` did not round-trip on a spindle torus's folded half.
  Each now returns the nearest point's parameters, half a turn round where
  that point lies past the apex or the axis.

- **Trimmed planar curves meet as their bases do.** `intersect_curves_2d`
  sampled any trimmed curve, and missed a line tangent to an arc, a line
  crossing one twice just below the tangent, and a circle touching an arc.
  It answers on the bases in closed form and clips to the trims' windows,
  as the space-curve path already did.

- **A curve's tangencies and lying stretches come back once.** Where no
  closed form applies, `intersect_curve_surface` reported a tangency up to
  eighteen times, and a curve lying in the surface (a cone's ruling, a
  circle on a sphere) as hundreds of piercings. Piercings between which
  the curve never leaves the surface are now one contact: a tangency once,
  or a stretch in `lying`.

- **A sphere on a torus's axis meets it in exact parallels.** The pair
  went to the marcher, and a sphere seated against the tube came back as
  two fitted pieces of its contact, a hundred times past their stated
  tolerance. It has a closed form: each place the sphere's meridian meets
  the tube's is a parallel, one tangential where the sphere is seated.

- **A fitted section states the error its images reach.** A section
  fitted beside a cone's apex reported a tolerance two hundred times
  tighter than its pcurve on the cone stood from the curve: the error was
  read at the trace's samples, and the chart turns fast between them
  there. Each pcurve is now lifted through its surface at dense stations,
  and the stated tolerance covers the worst.

- **A cone through its apex is sectioned on both nappes.** The closed
  forms for a cone against a coaxial cylinder, a coaxial cone and a level
  plane answered on the chart's own nappe only: a cone face running past
  its apex lost its far-nappe section, two parallel-sided cones meeting
  past an apex read as apart, and a level plane past the apex went to the
  marcher. Each crossing is now reported, and the cone's pcurve places a
  far-nappe parallel half a turn round; a face stopping short of its apex
  drops it by its window.

- **Arcs of one circle overlap across its seam.** The overlap of two
  trims of one circle was clipped within a single turn, so a stretch
  crossing the seam was lost, and a pair walked opposite ways reported
  parameters outside the second trim. The correspondence is now carried
  round whole turns, and every stretch both trims share comes back inside
  both.

- **A spindle torus is sectioned on its folded half.** A torus whose tube
  swallows its axis has, in every meridian, two tube circles reaching past
  the axis, and the closed forms against a plane square to its axis, a
  coaxial cylinder and a coaxial torus read only one: the parallels on
  the folded half were dropped. Each now reads both, and the pcurve on
  the torus takes such a parallel half a turn round.

- **Endlessly nested glTF and VRML are refused.** Both readers parse by
  recursion, and a document nesting arrays, objects or nodes a hundred
  thousand deep ran the stack out and aborted the process. Nesting past
  256 levels is refused by name.

- **A PLY face naming a missing vertex is refused.** `read_ply` indexed
  the vertices with the file's face indices unchecked and panicked on one
  past the end; it refuses the file by name, as `read_obj` does.

- **Flat caps of a fine mesh far from the origin stay whole.** A finely
  meshed part stored in single precision a long way out broke its flat
  faces into dozens: their small triangles lean past the coplanar angle by
  the rounding of their corners alone. A triangle's allowed lean now
  widens by its corners' rounding (twice the mesh's quantum over its
  smallest altitude).

- **A mesh's void stays a void wherever its ray lands.** The converter
  decided nesting with one ray, and a ray through an edge two of the outer
  piece's triangles share counted the crossing twice, or not at all, so a
  void came back as a separate solid. A ray grazing a triangle's boundary
  is set aside and another direction tried.

- **Unusable converter options are refused.** `solid_from_mesh` read a
  NaN, zero or negative `coplanar_distance`, a negative `quantum`, and a
  NaN or zero `crease` as something else, silently. Each is refused by
  name, as the documentation said.

- **Whole spheres and tori from single-precision meshes check valid.** The
  converter widened a whole face's tolerance to its fit's deviation but
  left its seam, poles and vertices tighter, which `check` refuses. The
  widening now reaches everything the face bounds.

- **A point in a spherical void classifies out of the material.** The
  bound of a sphere face used inside out covered only its seam: its chart
  outline put both seam walks on one column, so `classify_in_solid_exact`
  took points in a void bounded by it as inside the material. The first
  seam walk of an outline takes the image its own record names.

- **A pointed cone converts up to its apex.** `to_nurbs` and
  `general_transformed_shape` widened a cone face's window by a margin
  that ran past the apex into the other nappe, and the converted cone came
  back 15 to 20% too big. The margin stops at the apex.

- **A hole wound like its outline revolves hollow.** `make_revolution` read
  the hand of the outline alone, so a hole wound the same way swept walls
  facing into the material and added the hole's volume. Each wire's turn
  is measured, as the prism does: the outline must turn one way and every
  hole the other.

- **Faces welded within microns are not kept as a membrane.** Two boxes
  stacked a few microns apart, closer than the boolean welds vertices,
  fused into one shell holding both facing walls on the same edges: a
  wall of no thickness that `check` passed. Such a pair, bounded by the
  same edges and lying on each other facing apart, is dropped as the
  exact contact's opposed pair is, and a shell left with nothing else
  goes.

- **A solid poking microns through another is not taken as nested.** Where
  the general boolean refused, the nesting fallback read the boundaries at
  a mesh's samples, a tenth of a millimetre apart, and took a drum two
  microns wider than its box as lying within it. Each inner face now
  climbs its chart from the samples nearest each outer face, and a point
  found outside refuses the nesting.

- **Fillets whose balls overlap are refused.** `fillet_edges` accepted two
  edges of a narrow face (or a short wall's two rims) whose bands overlap
  across it, and built a solid that was neither rounding. A chain whose
  contact on a shared face stands inside another edge's band is refused;
  edges meeting at a corner still meet there.

- **A fillet the faces cannot hold is refused.** `fillet_edge` never
  checked that the ball fits: a radius of 100 on a 2 mm box deleted the
  solid and reported success. It now checks each face reaches the ball's
  contact, as `fillet_edges` did, and the check itself no longer skips a
  radius so large its side test stepped off the face.

- **An offset through the part is refused.** `offset_shape` and
  `make_thick_solid` returned inside-out solids when the offset or the wall
  ran a part's faces through each other; they now refuse, as `offset_faces`
  and `move_faces` did.

- **Located solids offset where they stand.** `offset_shape`,
  `offset_faces`, `move_faces` and `make_thick_solid` rebuilt a solid
  moved by a location in its nodes' own frame: wrong geometry that passed
  `check`, or refusals. Such a solid is baked into world coordinates
  first, as an instanced one already was.

- **Over-repeated end knots are refused.** `KnotVector::new` accepted a
  domain end repeated `degree + 2` times, where evaluation gave NaN; it is
  refused as a malformed vector.

- **Unclamped splines split into Bezier segments and elevate unchanged.**
  `to_bezier_segments` and `elevate_degree` assumed clamped ends and moved
  an unclamped or periodic curve by millimetres. The domain's ends are now
  raised to full multiplicity like every interior knot, and each segment
  takes the control points its span reads.

- **Unclamped and periodic splines reverse and rescale unchanged.**
  `KnotVector::reversed` and `reparameterized` clamped every knot into the
  domain, changing a spline whose outer knots lie beyond it by up to
  millimetres. Only the knots at the domain's ends are snapped now.

- **Sections winding round a periodic surface fit.** The closed sections of
  two crossing cylinders at an angle (an oblique pipe tee) came back
  millimetres off both surfaces with a stated tolerance in the hundreds:
  each loop's image on the branch pipe ends a period from where it
  begins, so the loop fitted as an open trace and stalled. Such a loop is
  now fitted closed in space with its images' join smooth across the seam
  (`fit::fit_points_joint_winding`) wherever the plain fit misses.

- **Spheres come back whole through STEP and IGES.** A sphere face bounded
  by its seam walked both ways, as an exchange file states it, read back
  with no region to mesh: whole spheres came back empty, and a sphere cut
  on its side or fused with another lost most of its face. Where a
  boundary meets a pole at one vertex from two sides of the chart, the
  readers now close it with an edge across the pole row
  (`closed_at_poles`). Both writers also dropped a solid's voids; they now
  write them (`BREP_WITH_VOIDS`, and the void shells of IGES entity 186),
  and the IGES reader honours the shells' orientation flags.

- **glTF read under a mirroring node stays right side out.** A node whose
  transform reflects (a negative scale, or a matrix with a negative
  determinant) left the triangles wound inside out and turned every normal
  inward; the winding is turned back and the normals follow the map.

- **Booleans with mirrored curved solids.** A mirrored cone gave wrong
  results that passed `check`, and mirrored cylinders, spheres and tori
  were refused. Baking a reflected placement restates each analytic surface
  on its frame's right-handed twin, and a seam that runs against its
  chart's direction gets its pcurves inside the chart window.

- **Primitives in a mirrored frame are sound.** A box, cylinder, cone,
  torus or wedge laid out in a left-handed frame came out inside out, and a
  sphere failed to build. A round primitive is built on the frame's
  right-handed twin about the same axis, and a box or a wedge on the frame
  with `x` and `y` exchanged, which spans the same corner; the box's faces
  keep the names the caller's frame gives them.

- **A mirrored cylinder or cone converts to a spline turning its way.** The
  section circles were built on a right-handed frame, so a mirrored
  cylinder's or cone's spline turned the other way and its normals
  pointed in; they are built on the surface's own frame.

- **An offset surface's closed form is the offset.** `OffsetSurface::analytic`
  rebuilt an offset plane on a right-handed frame, turning a mirrored
  plane's normal over, and an offset cone on its basis's frame, naming
  points a little along the axis from the offset's; the plane keeps its
  frame's hand, and the cone's frame rises with the offset, so each names
  the offset's points at the offset's parameters.

- **Mirrored revolutions and offsets are the images of the originals.** A
  surface of revolution carried by a reflection turned the other way round
  its axis, and offset curves and surfaces landed on the other side of
  their basis; each now names, at every parameter, the image of the point
  it named. `Transform2::preserves_handedness` no longer reads a negative
  scale, a half turn in the plane, as a reflection.

- **A negative scale reverses directions.** `Transform::apply_direction`
  (and `Transform2`'s) left a direction unchanged under any uniform scale;
  a negative one is a point mirror with a scale, and reverses it, so
  circles, planes and frames scaled by a negative factor keep the sense
  their points do.

- **`fix_shape` turns a located face in place (#101).** A face turned right
  way out was rebuilt on a new node without its location, so a void whose
  faces are placed by one (an assembly part read from STEP) moved out of its
  solid and still checked broken. The turned face keeps its location.

- **A fuzzy boolean takes its fuzz as the faces' coincidence.** Two faces
  lying within the fuzz of one another over what they share are one
  surface to `fuse_fuzzy` and `cut_fuzzy`, as they are within the
  tolerances the faces state: a prism of a converted part's top face along
  its rounds, left as facets each leaning a few microns, fuses and cuts at
  a fuzz wider than the lean.

- **A pad beside a wall leaning a few microns off it is sound.** Where
  two exact edges (lines, circles) come within a measuring width of one
  another, the boolean now holds them to their own tolerances, and takes
  them for one curve only over a stretch ending where one of them ends:
  a pad's side meeting a wall along its top edge and parting from it by
  a lean was read as running along the wall's edge for a tenth of a
  millimetre, and the pad's bottom edge was left short of the wall.

- **The kernel builds on its declared minimum Rust, 1.90, again.**
  `ogeom-bool` used an `if let` guard on a match arm, which Rust 1.90
  refuses; 0.6.0 built only on newer compilers.

### Added

- **A history says which faces came through a boolean unchanged.**
  `History::copy` records a modification that is an exact copy (a new node
  on the same surface within the same boundary) and `History::copy_of`
  reads it back. A boolean reports so every face the tool never touched,
  and the sew every face it only moves onto shared edges; composition keeps
  the mark only where every step copied the shape or left it alone. An
  application holding names, meshes or bounds for a face carries them
  across without comparing geometry.

### Changed

- **A boolean on a face with hundreds of holes is several times faster.**
  Splitting a face picked its probe points with a scanline per vertex
  height crossing every segment, and placed each hole by asking every
  pair of cycles about every other; on a plate with four hundred bores a
  small pocket took 580 ms, most of it on faces the pocket never touches.
  The scanlines now stop once their choice is settled and meet only the
  segments at their height, each hole's containers are found once, and
  the material test skips polylines the ray cannot cross. The probes and
  pieces are the same; the pocket takes 170 ms.

## [0.6.0] - 2026-09-30

A minor release: `cargo semver-checks` against 0.5.1 finds one breaking
change, a field added to `MeshSolidOptions`, and every crate packages and
verifies in a publish dry run. A mesh now converts in two steps if wanted,
as it is and then refined (`refine_solid`), and recognition is sounder
throughout: noisy flat faces, rounds tangent to their faces, thin rims,
crests and single-precision input. The boolean holds faces to the
tolerances they state, so a pad sketched from a converted part's face
fuses, cuts and commons with it (#99, #100).

A known limit: a pad whose sides run along a solid's walls at a tiny angle
(a converted solid's facets, tilted by the mesh's rounding, under a pad
sketched from them) leaves a sliver the boolean does not yet close; a
refined solid's walls are not facets, and pads along them fuse.

### Changed

- **`MeshSolidOptions` has a new field, `keep_vertices`.** Code building
  the options with a struct literal names it, or takes the rest from
  `MeshSolidOptions::default()`.

### Added

- **A converted solid refines in a second step.** `refine_solid` finds
  the cylinders, cones, spheres and tori among a solid's facets and
  rebuilds them, as `solid_from_mesh` with recognition would: a mesh can
  be converted as it is first (`recognize: false`) and refined on top.
  `MeshSolidOptions::keep_vertices` keeps every mesh vertex a vertex of
  the solid, so the first step keeps the samples along curved stretches
  that merged edges would drop. A recognized region lying on a larger
  neighbour's surface of its own kind now joins it, where two fits of one
  cone parted over their few vertices' slop.

- **STL precision reaches the mesh converter.** `ogeom::io::stl::read_with_quantum`
  returns the mesh with its encoding's rounding: single precision at the
  largest coordinate for binary files, the finest printed digit for ASCII.
  `MeshSolidOptions::quantum` takes it (or `single_precision_quantum` for a
  mesh held in `f32`), and the default coplanar distance is then never
  below twice it. An ASCII file written to six significant digits a
  hundred millimetres from the origin kept its fillets only as facets; it
  now comes back with cylinders.

### Fixed

- **A face met along one line by two faces of the other solid is split
  there once.** Where a section runs along one face's own edge while
  reading inside both faces by a hair, it now counts as that edge's split
  of the other face, as one reading outside does, and a second face
  meeting the other along the same edge (a chamfer beside the face it
  bevels) gives way to it. A prism of a chamfered slot plate's top face
  pushed through its bottom had its slot walls cut twice along the bottom
  rim, and the boolean left an open shell.

- **A pad sketched from a converted part's face fuses, cuts and commons
  with the part.** A face states how far it may stand off its surface (a
  face fitted to a mesh, by its fit's deviation), and the boolean now
  holds two faces to that: they lie on one surface where every point of
  the one over the other stands within what the two state, though their
  surfaces differ (a sketch keeping its points in single precision draws
  an arc through them on another centre and radius, a few microns off the
  wall it was taken from). The edges carried between the two, the
  pieces' partners on the shared surface and the pieces one stands in for
  are all read at that tolerance, and an edge running past the other
  surface's window is carried over the stretch it has on it. A pad from
  such a sketch along a converted part's walls was refused, a piece on the
  other's boundary with no partner face. A box and a sheared copy sharing
  a plane with it, refused as unresolved same-domain contact, now fuse,
  cut and common to their exact volumes.

- **A section leaving a face's boundary all but along it splits the face.**
  The boolean orders the strands round each junction of a face's
  arrangement by where they leave a small circle about it, not by their
  first steps: a section starting at a tangent can set out a hair outside
  the boundary before turning in, and read by its first step it crossed
  the boundary and the face was left with no piece at all.

- **A round's crest meeting a flat face is a circle.** A mesh's rim where
  a torus meets a flat face along its tube's crest, a few slops off the
  crest after rounding, was taken for a free curve: the circle's radius
  turns infinitely fast with its height there. A ring within the reach of
  the crest is now the crest circle, and the torus is built as a band
  between two circles rather than falling to facets.

- **A thin disc's rounded rim comes back a torus from a mesh.** When a
  curved region shed the flat face's facets it had taken in and the refit
  of what was left found nothing, the whole region went back to facets;
  its rim's normals, one-sided once those facets were gone, misled the
  fit. The surface already held now stands where it still holds the rest.

- **A round between two flat faces comes back tangent to them.** The mesh
  converter puts a curved region that meets two non-parallel planes across
  smooth edges on the cylinder tangent to both, its axis along their line
  of meeting and its radius the median of the vertices' own, wherever that
  holds every vertex within the coplanar distance. A narrow round's few
  rows of vertices fitted a cylinder standing a micron off its faces.

- **Noisy flat faces come back whole from a mesh.** The mesh converter's
  default coplanar distance also follows how far the mesh's vertices stand
  off the flat faces they lie on (measured across nearly coplanar edges),
  and a flat face is grown against a plane refitted to what it has
  gathered rather than its first triangle's, whose lean carried a large
  face's far side past the distance. A drilled block exported with noise
  came back with its top split in two; it now comes back its seven faces.

- **Rounds met by flat faces along a tangent come back exact.** A flat
  face's facets drawn across the tangent circle it shares with a round have
  their corners on the round; they are now peeled from the round's region
  (which is fitted again without them) instead of making the whole region
  fail. Two regions left on one surface and sharing an edge are one. A
  plate with rounded corners and a rounded bottom rim comes back on its
  eighteen exact surfaces, four tori among them, where it took a thousand
  faces; a thin disc with a rounded rim on its four.

- **A converted dome on a cylinder is sound.** A sphere centred on the axis
  of the cylinders and cones beside it kept its own angular origin, so the
  cap's seam and the cylinder's started a quarter turn apart on the circle
  they share; the solid failed the checker, and on a coarser mesh the
  sphere was not built at all. The sphere now takes the axis's frame.

- **A curved region no longer takes a flat face's facets.** The mesh is
  first cut into flat patches, and a curved region may take a patch much
  larger than its own facets only whole, with every corner on its surface.
  A round met along a tangent by a long flat wall took the wall's facets
  whose near corners lay on it; a finely meshed plate with rounded rims
  came back in a thousand faces and now in five hundred. A plane beside a
  recognized face takes its boundary from it, so the face must also keep
  its flat neighbours on their own triangles, or it is faceted.

- **A recognized region stays one connected piece.** Dropping the triangles
  of vertices a fit refused could leave a region in pieces joined by
  nothing, each then built as a face on the first fit's surface. The
  largest piece is now the claim, and the others stay free to seed regions
  of their own or fall to the planes. A perforated cup's conversion lost
  three spheres its mesh does not hold, and came to within 0.02% of its
  volume from 0.17%.

- **The mesh converter checks each body it builds, and gives up only what
  fails.**
  - A recognized face pointing against the triangles it replaces, or
    covering a part of its surface they do not (sampled on the face, each
    point within the rise of the triangles near it), is faceted and the
    part rebuilt.
  - A body that does not hold its mesh volume gives up its own recognized
    regions; other bodies keep theirs.
  - The volume is measured at the facets' mean offset rather than a
    millionth of the part: a converted part of four thousand faces took
    55 s and takes 5 s. A thin disc with a rounded rim, meshed coarsely,
    made the conversion fail outright; it now comes back a sound solid.

- **A mesh of bodies glued along an edge converts to closed solids.** An
  edge four triangles use (two blocks touching along it, as exporters
  write glued parts) paired nothing, and both bodies came back open
  shells. The bodies are now what the edges used exactly twice join, and at
  a crowded edge each body's own two triangles pair up.

- **`solid_from_mesh` no longer builds faces over the wrong part of their
  surface.** A recognized region whose fit refused a few vertices lost the
  triangles bringing them, which could leave it in patches joined by
  nothing; its face was built from one patch's boundary and the solid
  passed the checker with the wrong volume (0.85% on a perforated cup).
  - Each patch of a region is now a face of its own.
  - A face whose bounds reach well past its triangles is faceted and the
    part rebuilt.
  - The finished solid must hold the mesh's volume within what its
    recognized surfaces may add over their facets; past that recognition
    is withdrawn and the mesh comes back faceted, flagged by the new
    `MeshSolidReport::recognition_withdrawn`.

- **`solid_from_mesh` unfolds slivers folded under their neighbours.** An
  exporter that moves a vertex across a thin triangle leaves it facing
  against all three of its neighbours while its winding agrees with
  theirs; it became a face pointing into the material, and the solid
  failed the checker. The diagonal it shares with a neighbour is now
  swapped, choosing the best-shaped pair that faces with the rest.

- **A pad through a converted part (#100, in part).** A pad cut from a
  converted part's own face and pushed through its walls and chamfered
  bottom now fuses, cuts and intersects, matching the same booleans on the
  exact part. Six causes, each where the converted part's single-precision
  slop met a rule written for exact geometry:
  - An edge lying on the other solid's surface within the two faces'
    stated tolerances takes its closed-form image at that looseness.
  - One cylinder on two frames, where the fitted image refuses, has the
    owner's pcurve carried across the charts exactly.
  - An exact section meeting a line or conic boundary edge on a curved
    face counts a near miss within that edge's tolerance as a crossing.
  - A contact measured along an edge ends at the edge's vertex, not a
    measuring width past it; the overshoot left one face a sliver its
    neighbour did not have.
  - A contact strand is folded onto the face's branch with the leaning
    ray: its middle sits level with the vertices it paved, where the level
    ray's parity is wrong.
  - A prism's rails take the tolerance of the vertex they sweep.

- **A chamfer along a rounded rim removed too much; a fillet refused it.**
  A chamfer on an arc of a rim (a rounded corner's quarter) cut a wedge
  turned the whole way round the corner's axis, taking material on the
  far side too: the solid passed the checker with the wrong volume. It
  now turns through the arc alone, as the fillet's does. The probe that
  tells a rim that stops from an arc split off a whole rim read the arc
  the wrong way round when the cap faced against the circle's normal, so
  a fillet along the bottom rim of a rounded block was refused; it now
  reads the arc's own middle. Both blends along a chain of lines and
  quarter arcs now match the closed-form volume.

- **`solid_from_mesh` recognizes one-row chamfers and exporter noise.**
  - A chamfer meshed as a single row of triangles has its corners on two
    circles, which a sphere and a torus fit as well as a cone. The cone
    is now taken, and a cone fit reads its lean from radius against
    height.
  - Coaxial surfaces now share the axis of the best determined one: a
    cylinder's before a cone's, and a larger region before a smaller. A
    fit that does not hold on the shared axis is fitted again about it.
    Neighbouring rims then meet on one circle.
  - An axis within 1e-3 of square to a plane of the solid is made square
    to the largest such plane, where the region still fits. Rims on that
    plane are then exact circles, not fitted splines, and a feature
    sketched on it extrudes into exact cylinders.
  - When no `coplanar_distance` is given and the curved fits press
    against the default, the distance widens (up to four times). A mesh
    whose exporter moved its vertices a little further no longer leaves
    thin planar slivers at the foot of each fillet.

- **A pad within a converted solid fuses (#99).** A pad sketched on a
  converted part's own outline and sunk into it runs along the part's
  walls and touches its recognized ones at every vertex, within the
  tolerance the conversion states; the boolean read that contact as
  ambiguous or as edge contact and refused. Where it refuses, `fuse` and
  `common` now ask whether one solid lies within the other at the
  looseness the solids state, and answer with the outer one or the inner
  one. A piece whose probe no ray can read is also asked at its next probe
  before the boolean gives up.

- **A traced section round a whole loop is no longer taken for a wandering
  trace.** Two surfaces meeting all the way round (a near-coaxial cone and
  cylinder), whose faces hold only an arc of the loop, gave a section
  longer than four times the faces' extent, and the guard against traces
  lapping beside a chart's pole refused it. The guard now allows a turn
  round the section's own extent as well.

## [0.5.1] - 2026-09-29

A patch release: `cargo semver-checks` against 0.5.0 finds no change
needing a bump in any of the fifteen crates, and every crate packages and
verifies in a publish dry run. It adds direct editing of a solid's faces
(`offset_faces`, `move_faces`), and a mesh converted with recognition
keeps its recognized faces on the triangles they replace.

### Added

- **Direct editing of a solid's faces (#98).** `offset_faces` offsets the
  named faces along their outward normals and `move_faces` moves them by a
  translation or a rotation; the faces around them follow, each staying on
  its own surface with its edges re-derived where the moved faces now meet
  it. A plane moves parallel, a cylinder, cone, sphere or torus changes
  radius; several faces move together; the history maps every face to the
  face it became. An edit that runs a face past the faces across from it
  or into the rest of the solid is refused by name.

- **`check_self_intersection_near`** tests only the face pairs holding one
  of the faces named: after an edit, the pairs it left alone cannot have
  begun to cross.

### Fixed

- **`offset_shape` of a bored block and of a ball.** A seam reached
  reversed (a bore a boolean cut) was rebuilt the wrong way round, and a
  sphere's poles and two-pole band had no rebuild; both now offset.

- **A surface recognized in a mesh stays on the triangles it replaces
  (#97).** A region was accepted when its vertices lay on the fitted
  surface, so long facets across a flat stretch between two rounded
  corners could be taken for a torus through the corners, whose face then
  stood far outside the part and turned the classifier's rays: the solid
  checked inside out. Every triangle of a region must now stay on the
  surface as closely as its own turn allows, and an edge between two
  recognized faces is not traced where either surface turns through a
  right angle along the mesh chord it replaces. `tight_bounds` no longer
  counts a face's surface beyond the chords its chart triangulation runs
  on, which put a dome a few hundredths below the edge that bounds it.

## [0.5.0] - 2026-09-29

A minor release, since one public enum gains a variant: `cargo
semver-checks` against 0.4.1 finds `PipeLaw::Fixed` added to an exhaustive
enum, and nothing else in the fifteen crates. It adds helical sweeps of a
profile off the axis, pipes and skinned lofts that close to a point, and a
pipe law that carries its section by translation. Unevenly scaled spheres
and revolved solids now mesh, measure and cut as the ellipsoids they are,
through their axis and on their seam's plane included; a profile turned
past a quarter turn revolves outward; and `check` finds a solid inside out
as a whole.

### Changed

- **`PipeLaw` has a new variant, `Fixed` (#90).** A `match` over
  `PipeLaw` needs an arm for it.

### Added

- **Helical sweeps of a profile off the axis (#89).**
  `make_helical_sweep` takes a planar face in any plane clear of the axis,
  not only one through it: a level profile square to the axis climbs as a
  ramp or a coil of a flat section. The caps are the profile's plane and
  that plane turned with it. A taper is taken for a profile through the
  axis or square to it; a profile the axis runs through, or one the screw
  carries partly forward through its plane and partly back, is refused.

- **A multisection pipe closes to a point (#92).** `make_pipe_sections`
  takes a vertex as its first or last section where the spine starts or
  ends: the section shrinks onto the spine and the skin closes on the
  point, uncapped there. Down a straight spine it is the exact pyramid or
  cone over the section.

- **`PipeLaw::Fixed` (#90):** the section keeps the frame it has at the
  spine's start, carried along the spine by translation and never turned.
  Each section edge sweeps the exact B-spline `edge(u) + spine(v)` (the
  control points summed, the weights multiplied), so nothing is fitted. A
  spine that crosses the section's plane one way and back the other is
  refused, as is one running in the plane.

### Fixed

- **A revolved solid scaled on every axis meshes (#96).** Converting a
  face with a pole under a general transform looked for the pole's row on
  the new surface at the pole's old place, so a scaling that moved the
  pole put its degenerate edge on the wrong end of the chart and the face
  enclosed nothing. The row is now found at the pole's mapped place: a
  hemisphere scaled to radii 8, 5, 3 meshes and measures as half the
  ellipsoid. A pcurve refitted by the conversion now holds its target in
  space, through the chart's stretch and between its samples, so a quarter
  turn scaled that way is integrated exactly instead of falling back to a
  mesh 4% light.

- **A profile turned past a quarter turn about the axis revolves
  outward.** A profile placed by a rotation had its missing pcurves read
  from its curves as placed on its surface as stated, so past a quarter
  turn `make_revolution` built its solid inside out, and at exactly a
  quarter turn refused it as lying on the axis. The curves are now read in
  the surface's own frame.

- **Cutting a solid of revolution or an ellipsoid through its axis.** A
  plane through a revolution's axis now meets it in the profile turned
  onto the plane, exactly, where marching stalled at the poles. A plane
  holding whole columns or rows of a spline surface meets it in those iso
  lines, so an ellipsoid from a scaled sphere halves on the plane of its
  seam. A point is classified against a revolution whose profile crosses
  the axis at either of its two places in the chart, so the half of the
  profile a face does not use no longer reads interior points as outside.

- **`check` finds a solid inside out as a whole.** Faces whose flags all
  agree with each other were taken as facing out, so a shell turned
  inside out everywhere read as valid. Where the flags agree, one face of
  each shell is now probed against the material; a shell found inside out
  is reported face by face, and `fix_shape` turns it back.

- **A fillet chain holding closed edges rounds each (#94).** A piece of a
  rim the solid holds in parts is rounded the whole turn round, so the
  rim's other pieces (a second half passed in, or a whole rim a boolean
  split in two on an earlier step) were gone when their turn came and the
  chain refused as though its members interfered. A member lying on a
  circle the chain has already rounded now counts as done. Eleven fillet
  stress cases refused before now round.

- **Booleans of whole revolutions whose seams coincide off the profile
  plane (#95).** A piece of a face was classified at a probe placed by the
  widest scanline interval across it; a scanline in the thin gap just past
  an arc's chords, and inside the arc, crosses no chord of it and reports
  an interval as wide as the piece, with its midpoint in a hole. Room is now
  the lesser of the interval's width and its scanline's gap. Four fillet
  stress cases refused before now round.

- **A skin narrowing to a point along a bend was built inside out.** Its
  wall was turned away from the centroid of every row, which a bent skin
  puts outside itself; it is turned away from its own section instead.

- **A skinned loft through several sections closes to a point (#91).** A
  family fit refining past one control point per datum left spans with no
  data and gave up on its best earlier round, which missed the tolerance;
  it now takes the interpolating spline there. A loft's sections are read
  off chords kept within a micron of their curves instead of 64 per edge,
  so a loft of circles is a solid of revolution to its tolerance.

- **Unevenly scaled spheres and revolved solids are what they look like
  (#93).** A revolution whose profile runs on past its pole (a whole circle
  turned) was converted to a spline over the whole profile, and its pole
  put at the chart's end instead of where the surface reaches it; each
  face is now converted over the rows it uses, a pole is placed where the
  surface reaches it, and under an affine map whether the conversion
  turned the normal is measured before the map. A partial sweep's two
  sides share their edge on the axis through the conversion. A boolean
  with a sphere spelt as a spline no longer tears a closed section at the
  seam, so a scaled sphere cuts as the ellipsoid it is.

## [0.4.1] - 2026-09-28

A patch release with no API change needing a bump; `cargo semver-checks`
finds nothing against 0.4.0. A boolean's history now reaches every face of
its result, and a fillet's or chamfer's history its blend or bevel face.

### Fixed

- **A boolean's history reaches every face of its result (#88).** The
  pieces were recorded, but sewing rebuilds them onto shared edges and
  its record was not carried on, so only faces the sew left alone were
  reached. A fillet or chamfer also records its blend or bevel face as
  generated from the edge it replaces.

## [0.4.0] - 2026-09-28

A minor release, since two public types change shape; `cargo
semver-checks` against 0.3.4 finds those two and nothing else. It follows
an adversarial review of the kernel: hostile files are
refused instead of aborting the process, results of several shells nest
by where they are, faces built inside out are found by `check` and
turned by `fix_shape`, and most of the time went out of booleans,
meshing, mass and fitting (the stress harness runs in 20 s where it took
143 s, every outcome the same or better). It adds revolutions up to a
face, arc joins for thick solids, pipe frame laws, corner styles and
multisection pipes, flat spirals, sections of faces and open shells, and
seams for imported periodic faces, with which booleans cut parts they
refused before.

### Changed

- **`FixReport` has a `faces_turned` field**, so a literal building one
  needs it.
- **`step::parse::Arg::Typed` holds its keyword and arguments boxed
  together,** `Typed(Box<(String, Vec<Arg>)>)` instead of
  `Typed(String, Vec<Arg>)`, which makes every argument 32 bytes instead
  of 56.

### Performance

- **`check` walks the shape once** for every type it examines, instead of
  once per type.

- **A boolean places each edge's curve once,** shared by the faces on
  either side of it, instead of copying and transforming it per face.

- **An exact mass integral holds at most a million samples.** Each run's
  comparison measures are summed as its samples are drawn; a run past the
  bound is drawn a second time into the result instead of held, the same
  samples in the same order.

- **Fewer surface evaluations per Newton step.** A march's correction
  takes the walked point and its gradient from the evaluation its system
  already made (four evaluations a step instead of six, the same values),
  and extrema read each surface's point and derivatives from one jet.

- **A walk's tangent is found without allocating:** the null vector's
  minors are expanded in place, with the same arithmetic, and each accepted
  state is stored once instead of copied.

- **A node's identity is a slot lookup,** not a hash: the table is a
  vector indexed by node, cheaper to copy for undo, and `shape_of` answers
  in node order.

- **Sewing bins its few loose edges apart,** so one edge with a wide
  stated tolerance no longer coarsens the twin search for every other.

- **A planar section skips faces wholly on one side of its plane** by
  their boxes, before intersecting their surfaces.

- **Meshing measures each chart edge's sag once** across the refinement's
  rounds and the triangles sharing it, and a fillet chain finds its next
  edges by projection only among edges whose boxes reach them.

- **A face's arrangement has no quadratic steps left.** Strand ends snap
  through a grid, dangling chains peel through a queue, each dart's angle
  is computed once before sorting, and hole nesting reads each cycle's
  nodes, area and box once and tests boxes before polygons. Every choice
  is the one the old scans made.

- **A march's seeding tests each cell only against cells it could meet,**
  binned by box, in their own order, instead of every cell of the other
  surface.

- **A STEP argument is 32 bytes instead of 56:** the rare typed value's
  keyword and list are boxed rather than sizing every argument of a file.

- **A boolean splits and classifies its faces in parallel,** each face's
  pieces and junctions joined in face order, and intersects a contact
  edge only with target edges whose boxes it can reach.

- **Checks and bounds repeat less.** `shape_bounds` bounds each shared
  edge and vertex once per call; the face-orientation agreement places
  each edge's curve once and keeps one golden-section point per round;
  `check_self_intersection` gathers each face once, skips pairs whose
  boxes stand apart, and keys adjacency by node rather than by a hash
  that could collide.

- **`cells` runs one boolean, not three:** the three cells are three
  choices of the same classified pieces. Pipes and evolved solids fuse
  their pieces in a balanced order, neighbours with neighbours, instead of
  each onto everything before it.

- **The intersectors' small Newton solves allocate nothing.** Curve
  crossings, curve and surface extrema, curve-surface polishing and the
  marcher's contact solve run on `newton_system_fixed`, the same damped
  iteration on stack arrays; a curve crossing reads each point from its
  derivative table.

- **Hidden lines test each sample against the triangles and faces in
  front of it only.** The mesh drawing bins triangles by where they
  project; the exact drawing skips faces whose projected box misses the
  sample or that lie wholly behind it. Mesh-to-solid finds the triangle
  across a sliver from an edge map instead of a scan.

- **The STEP reader borrows its arguments instead of copying them,** and
  the writer finds a coloured product's entities by node instead of
  scanning every written entity per product.

- **Projection onto a plane, cylinder or sphere starts from its closed
  form** where the foot lies in the window, instead of a 33 by 33 grid.
  `SolidMesh` meshes a solid once to classify many points; the fillet's
  run-out probes use it instead of meshing per probe.

- **Smaller savings.** A surface's derivative table is filled only to the
  total order asked for, each control row summed across `v` once; the
  thick-solid edge scan maps edges to faces in one walk; `ancestors_of`
  searches each candidate for the target's type only; `fix_shape` skips
  measuring an edge whose ends already stand further apart than the reach.
  A marched fillet's chord and the offset guard's finer retry follow the
  part's own size, not fixed lengths.

- **Spline basis derivatives run on fixed arrays below degree eight,** bit
  for bit the general recurrence without its nested scratch vectors: a
  third off the thread-groove tests. Curve crossings skip segment pairs
  whose boxes stand further apart than the reach.

- **Exact volumes of spline faces cost a third of what they did.** The
  chart integral reads a polynomial patch along each inner integral's own
  `v` from the net summed across `v` once, and a polynomial patch's first
  and second derivatives skip the rational quotient. `Surface::point_d1_at`
  gives the point with its first derivatives from one evaluation.

- **A boolean intersects its face pairs in parallel,** each pair's
  findings joined in pair order so the result is the same at any thread
  count. A parallel stage inside another now runs on the worker it lands
  on instead of spawning threads per outer item, and the machine's thread
  count is read once.

- **Spline fits solve their normal equations instead of inverting them.**
  An open curve's normal matrix is banded, and a banded Cholesky costs its
  size times the band squared where the inverse cost its size cubed. A
  system so near singular that the data does not decide a span's control
  points now reports itself singular, as an exactly singular one did,
  instead of solving to control points anywhere.

- **Curve crossings seed only where the curves come nearest.** Two curves
  running near each other (a boolean's section beside the edge it was cut
  along) put hundreds of segment pairs in reach, and each was polished by
  Newton to the same few crossings. Seeds are now the pairs nearest along
  either curve, and stretches the curves share are not polished at all.
  The stress harness runs in 16 s instead of 143 s with every outcome
  unchanged; one corpus drill went from 18 s to 1.5 s.

- **Walks, classification and evaluation do less of the same work.**
  `explore` no longer descends below the type it is looking for; a
  placed child is built in one step; a solid's prepared boundary answers a
  point near a face from what it prepared, not by rebuilding the face's
  trim; interior mesh points are tested against a banded ring index;
  de Boor evaluation and mass simplices use stack storage; the identity
  transform returns its input; the walker measures each tangent once; the
  boolean builds each face's chart outline and trim samples once.

### Added

- **`fix_shape` turns faces that face into their solid's material,** and
  reports how many (`faces_turned`); `inside_out_faces` names them.

- **`seam_periodic_faces`** gives imported periodic faces the seam or pole
  edge their chart needs: a band between rings that each go round once
  (the whole circle re-anchored and the other ring's crossing edge split
  where needed), a sphere's cap bounded by one ring, and a sphere's face
  reaching a pole with no edge along the pole's row. Each faces the way
  its loop runs round the chart. `reanchor_periodic_rings` runs it after
  the rings. `make_band_of_rings` builds a band from rings of any number
  of edges.

- **`native::compacted`:** a fresh model holding only what given roots
  reach, for reclaiming what a long session's discarded and failed
  operations left in the model.
- **`Document::set_undo_limit`:** at most 256 checkpoints are kept by
  default, the oldest dropped first.

- **`solve::newton_system_fixed`:** the damped Newton iteration for a
  fixed number of unknowns, allocation-free.

- **`SolidMesh`:** a solid's boundary meshed once, to classify many points
  against.

- **`Surface::point_d1_at`:** the point and first derivatives from one
  evaluation, for callers that need both at every sample.

- **`fit_surface_grid_sections`:** a grid fit with rows placed across by
  chord length, for sections lofted at spacings of their own.

- **`check` finds faces turned inside out.** A flipped face leaves every
  edge used twice, so nothing topological saw it and a mesh volume repaired
  it silently. Where the faces' flags disagree along their shared edges,
  each face is probed a step off both sides; a face whose outside is
  material is reported broken.

- **Revolutions up to a face (#73).** `make_revolution_until` turns each
  point of a profile about its axis until its circle first meets a limit
  face's surface, refused by name where a circle never meets it or meets
  it only past half a turn.
- **Arc joins for thick solids (#75).** `make_thick_solid_with` takes a
  `Join`: the walls extended to meet, or rounded about each edge convex on
  the side they grow toward, a cylinder of the wall's thickness about the
  edge and a ball at a corner.
- **Pipe frame laws (#76).** `make_pipe_shell_law` and
  `make_pipe_shell_with` take a `PipeLaw`: rotation-minimizing, Frenet, an
  auxiliary spine the section's normal points at, or a fixed binormal.
- **Extended and round pipe corners (#77).** `make_pipe_shell_with` takes
  `PipeCorners`: the mitre, legs run on past the corner and fused, or the
  section turned about the corner. Built for a face down an open spine of
  straight legs.
- **Multisection pipes (#78).** `make_pipe_sections` sweeps several planar
  sections down one spine, each standing where the spine crosses its plane
  and blended between them in the spine's moving frame.
- **Flat spirals (#79).** `make_helical_sweep` with a pitch of zero and a
  taper per turn turns the profile in its plane square to the axis while
  moving it out, refused where one turn would meet the last.
- **Pipe shells round straight corners are exact.** A face swept down an
  open spine of straight legs with corners is each leg's prism, run on
  past its corners, trimmed at the mitre planes and fused: planes and
  drums, measured exactly, whatever the profile's seam (#82).
- **Sections of a face or an open shell by a plane (#84).**
  `section_face` gives the edges where a face, shell or compound of faces
  crosses a plane, each trimmed to the face it lies on: exact where the
  surface and the plane have a closed form, fitted otherwise, a face's own
  boundary where it lies in the plane. `section` takes such a sheet
  against a half space bounded by a plane.

### Fixed

- **`canonical_simplify` rebuilds faces whose edges are used reversed.**
  An edge it left unchanged kept its occurrence's direction and had it
  applied again, so a bore spelt as a spline running the other way round
  failed with "edge 0 ends where edge 1 does not begin".

- **Mirror-placed bodies measure exactly.** The exact mass path refused
  every reflecting placement, so a mirrored part was measured from its
  mesh (a small drum 2.5% out at the default chord); a reflection is now
  taken with the chart's normal turned back. The bound-filter audit
  reports an unfiltered fill that refuses instead of failing the boolean.

- **A face's chart outline closes across the openings its edges' ends
  leave**, so a point level with one no longer reads inside the face. The
  boolean's paving kept a section piece off a narrow torus wedge that way;
  the bound-filter audit (`OGEOM_BOOL_AUDIT_BOUNDS`) now passes on the
  curved-corner blends.

- **Booleans on parts with seamless periodic faces.** Every face of a
  part is arranged, so one such face refused any boolean on it; three of
  the six stress drills into `nist_ftc_06` now cut. A spherical pocket
  imported facing out of the sphere comes out of healing facing in.

- **Spline walls whose pcurves stray past the chart's border by rounding
  are measured exactly.** The exact volume path and its orientation check
  bring such points back into the surface's domain instead of falling back
  to the mesh, which read a binormal-law pipe 1.5% high at the default
  deflection (#87).

- **`make_pipe_sections` matches each section's start and sense to the
  section before** in the spine's frame, so sections drawn from different
  seams or running opposite ways blend without a twist (#86).

- **An auxiliary pipe guide ending a hair short of the last station** (a
  sketch's end in single precision) is carried along its end tangent to the
  station's plane when within the sweep's tolerance, instead of refused
  (#85).

- **`check` names each part too tight for its bound once.** A vertex
  reached through both of a face's edges was reported twice.

- **STEP names past ASCII round-trip.** The writer spells them in `\X2\`
  escapes and doubles a backslash; the reader decodes `\X2\`, `\X4\`,
  `\X\`, `\S\` and `\\`, and reads raw bytes as UTF-8 where they are.
  A number that is not finite is refused rather than written as `NaN.0`.
- **Smaller fixes.** `section` shares one vertex where consecutive pieces
  meet, so the section is a wire's worth of edges. A pair of faces marched
  whole because one curve had no chart image no longer drops the pair's
  later curves or repeats its earlier ones. A strand end inside two
  junctions names the nearer. `check_tolerances` holds a face to its edges
  across the wire between them. Area and volume fall back to the mesh
  alike where the closed forms cannot evaluate, and both report a
  cancelled watch. Exact hidden lines draw a closed silhouette to the
  deflection instead of as an octagon. A face edge with only a pcurve is
  meshed at the chord the faces agreed on.

- **A loft through unevenly spaced sections sagged between them.** The
  skin's parameters across the sections ignored their spacing, so a square
  twisting through steps of 0.625 and 0.3125 measured a fifth light or
  would not fit. Sections are now placed across the skin by chord length
  (`fit_surface_grid_sections`). A section of a face whose ends two
  neighbouring faces place a hair apart shares one vertex that reaches
  both.
- **A frame-law pipe asked every edge of its spine for the whole station
  count,** and a law whose frame jumps doubled its stations twelve times.
  Stations are now shared out by length, and a law that outruns 2048
  sections is refused.

- **Trims scaled with their basis only sometimes.** A trimmed surface
  whose basis domain started at zero, and a trim of a trim of a line, kept
  their old parameters under a scaling. Trims are now carried affinely from
  the basis's old domain to its new one.
- **A tangency contact on a full drum read as outside its face.** The
  contact's trim test saw only one seam column; it now reads the same
  welded outline every other trim test does.

- **Faces built inside out, found by the new check.** A fillet that adds
  material (a concave seat) on a marched blend faced its blend band into
  the fill; a draft about an oblique neutral faced its wall inward whenever
  its hinge chained against the old surface's direction; a skinned tube
  round a bend in any plane but XY faced its wall into the bend;
  `to_nurbs` of a mirrored part turned its bore walls; canonical
  recognition of a spline facing its axis kept a flag the cylinder no
  longer agreed with. Each now measures the side it builds on.
- **Newton on an unevaluable point read as a root.** Where a surface or
  curve could not be evaluated mid-solve, the residual fell back to zero
  in the walker and to a made-up point elsewhere; it is now infinite, so
  the damped step backs off.

- **A part standing in another's notch was swallowed as a void.** Results
  of several shells were nested by bounding box, so a disjoint part inside
  another's box became a void of it, and an island in a cavity a second
  void. Shells are now nested by classifying one against the other, by
  depth: a shell inside a void is a solid of its own. `make_volume` shares
  the same nesting.
- **Two planes a microradian apart were called one plane.** Parallel was
  decided on one minus the cosine, a million times looser than the angle it
  was compared with. Planes are now parallel when they part by no more than
  the confusion distance across the region their windows cover.
- **Surface extrema never polished a farthest point, and reported stalls
  as extrema.** Seeds are the local minima and maxima of the near and far
  fields over each surface's lattice, thinned apart, and every polished
  point in curve and surface extrema is checked for stationarity by angle.
- **A line touching a sphere or a cylinder was missed about half the
  time.** A tangent discriminant that rounding leaves a few ulps below zero
  is now the tangency.
- **Smaller fixes.** Periodic curves and surfaces refuse a parameter that
  is not finite instead of answering a point that is not a number. Two
  boolean contact paves cap an imported edge's tolerance as the others do.
  `fuse` decides half spaces at the caller's tolerances. Sewing carries
  pcurves across in a fixed order, so repeated runs build identical tables.

- **Malformed and hostile files are refused, never aborted on.** The STEP
  reader no longer panics on entities with too few arguments, overflows its
  stack on entities that refer in a circle or on deeply nested argument
  lists, or allocates by the largest entity id, a knot multiplicity or a
  degree in the file. IGES reads its fixed columns as bytes (a label past
  ASCII no longer panics), bounds Hollerith lengths and every count by the
  entity's own parameters, and refuses transform and curve cycles. The
  native, BREP and glTF readers reserve no more than the file holds; 3MF
  archive offsets are bounded by the file and a deflated entry is refused
  the moment it inflates past its promised size. A placement's power costs
  its bit length, not its size.

- **Sweeps sampled a placed section where it was stored, not where it
  stood.** A section wire carrying a placement was read at its unplaced
  position when a pipe or loft sampled it.
- **A boolean between a helical sweep and a box whose face holds the helix
  axis did not close (#81).** The sweep's walls were fitted a quarter turn
  at a time with their borders on the profile's own plane and the planes
  square to it, so such a face met the walls only along those borders,
  where no section starts. The borders now stand an eighth of a turn off
  those planes, and the box keeps exactly half the coil either side. A
  junction between fitted sections now closes in the exact integral within
  both sections' tolerances or the vertex's, whichever is wider.
- **A helical sweep trimmed by a boolean measured about 3% out (#80).**
  The fitted sections the boolean leaves on its walls meet within their
  own stated tolerances, and the exact integral refused them for missing
  a millionth of the chart, falling back to a coarse mesh. Chart loops now
  close, and pcurves lie on their edges, within those tolerances, so the
  common and the cut measure their closed forms at the default deflection.
- **A groove cut into a bore drilled by a cylinder primitive failed
  (#74).** The primitive's wall is stored over exactly its own height, so a
  section crossing its rim stopped a hair short of it and the cut refused to
  close. The boolean now stretches a cylinder's window a little past an end
  the other face reaches beyond.
- **A half space clear of a solid sectioned it anyway (#83).** The box
  standing in for a planar half space reached twice the other shape's
  diagonal from the plane, so a plane further off than that put the box's
  far face inside the shape: a section returned edges and a common trimmed
  material. The box now reaches past the shape's far side.
- **Re-anchoring periodic rings turned rebuilt faces inside out.**
  `reanchor_periodic_rings` reversed a reversed face once when rebuilding
  it and again when rebuilding its shell, so the face came back facing
  into the solid. `check` passed it and a mesh volume did not notice, but
  a boolean's exact volume lost the material behind it: a drill through
  an imported plate measured 1.5% light. Found by the new robustness
  harness.

## [0.3.4] - 2026-09-26

A patch release with no API change needing a bump; `cargo semver-checks`
finds nothing against 0.3.3. It adds typed DXF reading, helical sweeps,
exact pipe shells along lines and arcs, the medial axis of any planar
face and the middle path of a pipe-like solid, dividing and remodelling
shapes, NURBS conversion within a tolerance, tight bounds, VRML reading
and more of STEP and IGES. Fillets pinch at tangent poles, round corners
where a curved face meets and close in more orders at apexes and curved
corners, and a fillet or chamfer too wide for its face is refused.

### Added

- **Typed DXF reading (#72).** `read_dxf_entities` reads lines, arcs,
  circles, ellipses, splines and polylines with their bulges and closed
  flags, the drawing's `$INSUNITS`, and which layers and linetypes are
  dashed. `read_dxf` now closes a closed polyline and no longer reads a
  `POLYLINE` header's elevation point as its first vertex.
- **The middle path of a pipe-like solid.** `middle_path` walks a solid
  from one named face to another and returns the curve through its
  cross-section centroids: a line for a straight tube, a fitted spline
  otherwise. Each point is the centroid of a section cut square to the
  path, holes subtracted. The result carries the measured deviation, held
  under the tolerance asked for.
- **Fillets pinched at a tangent pole.** A crease whose two faces turn
  tangent at one of its ends, such as the seam of two equal drums crossing
  where the drums touch, blends: the band narrows to a point at the pole.
  Both the fused drums' notch and their common part's ridge round, one
  seam edge at a time or all in one call.
- **Corners where a curved face meets round.** `round_vertex` rounds a
  convex vertex whose three faces include a drum, cone, ball or spline:
  one ball touches all three and its patch closes the corner. `fillet_edges`
  on the three edges of such a corner closes it with that ball.
- **Exact fillets where a wall meets a drum along a ruling.** A straight
  edge between two planes or drums parallel to it (the vertical edge where
  an extruded profile's line meets its arc) blends with an exact drum band,
  and no longer fails for want of a planar pair.

- **Pipe shells round a skew corner against a curved leg.** Where a
  spine turns a curved leg's end out of its plane, the sweep is built in
  pieces: each side runs on straight past the corner, is trimmed by the
  mitre plane and the pieces are fused. This was refused.

- **The medial axis of any planar face.** `medial_graph` returns the
  axis of a face with holes, reflex corners and curved edges as a graph of
  branches, each the exact bisector of its two boundary elements: lines,
  parabolas, hyperbolas and ellipses for segments, corners and circular
  arcs, and a fitted curve held to the tolerance for any other edge. Each
  branch point carries its clearance, and every branch is checked to keep
  no other boundary element nearer than its own two.

- **More of IGES reads.** A conic arc whose axes turn, parametric
  spline surfaces (114) as exact bicubic B-splines, offset curves whose
  distance varies, and trims given only in a B-spline surface's parameters.
  Model-space notes, leaders, labels, symbols and dimensions read as PMI
  callouts, a dimension whose text states a number also as its value.
  Subfigure instances place their shared definition, and levels and
  groups read as layers.
  Constructive solids read too: the primitives, solids of revolution and
  extrusion, boolean trees, assemblies and instances.
- **Small faces and small solids go.** `fix_small_faces` collapses a face
  that fits in a ball of the given size to a point, and a thin strip to one
  of its long sides, rebuilding the faces around it. `remove_small_solids`
  drops the solids enclosing less than a given volume.
- **More STEP geometry reads.** Curves: Bezier, uniform and quasi-uniform
  splines, hyperbolas, parabolas, trimmed curves, offset curves,
  polylines, composite curves and curve replicas. Surfaces: surfaces of
  revolution, offset surfaces, rectangular trimmed surfaces, curve-bounded
  surfaces, degenerate tori, surface replicas, the Bezier and uniform
  spline forms, and a rectangular grid of patches joined into one spline.
  The writer puts hyperbolas, parabolas, offset curves, surfaces of
  revolution, extrusions and offset surfaces out exactly instead of as
  splines.
- **An offset of an analytic surface is that analytic surface.**
  `OffsetSurface::analytic` returns the plane, drum, cone, ball or torus
  the offset is.
- **A trimmed surface converts to a B-spline** as its basis over the
  trim's window. `BSplineCurve::segment` and `BSplineSurface::segment`
  cut a spline to a piece of itself, parameters kept.
- **Dividing shapes.** `divide_face` cuts a face along a line of its own
  parameters, with the surface's exact iso-curve as the new edge.
  `divide_by_continuity`, `divide_by_angle` and `divide_by_area` cut edges
  and faces at weak knots, past an angle or above an area. `to_bezier`
  leaves every curve and surface a single Bezier span.
- **Every face divides.** A face with no boundary of its own (a whole
  ball or torus) is first bounded by its chart's sides, seams and poles
  included; a face on a surface with no closed-form iso-curve (an offset
  of a spline) is cut along one fitted at its own parameters.
- **IGES trims on analytic surfaces.** A trim given only in the
  parameters of a plane, cylinder, cone, sphere or torus is lifted through
  the format's parameterization, and those surfaces take the reference
  direction the file names for where their angle starts. So is one on a
  tabulated cylinder, and on a surface of revolution of a line or spline.
- **Remodelling.** `restate_geometry` rebuilds a solid with its surfaces
  and curves restated by the caller, every trim re-derived.
  `restrict_degree` refits every spline above a degree at that degree,
  within a tolerance. `swept_to_elementary` names an extruded or revolved
  line or circle as the plane, drum, cone, ball or torus it is.
- **NURBS conversion within a tolerance.** `to_nurbs_within` converts a
  solid exactly where a closed form exists and fits the rest (offset
  surfaces, helices, offset curves) within the tolerance given, where
  `to_nurbs` refuses. `SurfaceGeometry::fitted_bspline` fits any bounded
  surface.
- **A curve continues to a named point.** `BSplineCurve::extended_to`
  joins on a piece that carries the curve's end derivatives to the order
  asked and ends at the point. `BSplineCurve::extended` now meets the
  length asked exactly rather than to first order.
- **Global minima.** `global_minimum` finds the least value of a
  function over a box by branch and bound, certified to a tolerance under
  slope and curvature bounds it estimates; `swarm_minimum` is a particle
  swarm; `minimize_local` is a Nelder-Mead descent in several variables.
- **VRML reads.** `read_vrml` reads VRML 2.0 and 1.0 scenes into placed
  meshes: face sets and the box, ball, drum and cone primitives, with
  `DEF`/`USE`, transforms, switches and material colours honoured.
- **Helical sweeps.** `make_helical_sweep` turns a profile lying in a
  plane through an axis along a helix about it (a thread or a spring):
  every point runs its own helix, cylindrical or tapered, either hand.
- **Pipe shells along lines and arcs are exact.** A spine of lines and
  arcs meeting tangent sweeps as extrusions and revolutions of the
  section, fused, so its walls are planes, drums, cones and tori and the
  solid measures in closed form.
- **A wire swept round a skew corner.** A closed planar wire profile
  sweeps round a corner that turns a curved leg out of its plane as the
  walls of the face it bounds, where it was refused.
- **Skinned lofts keep corners and circles** (#57). Two sections loft
  ruled and exact; more with matching corners loft as one strip per
  edge, corners kept; coaxial circles loft as a solid of revolution. Two
  squares give the prism and two circles the drum.
- **Tight bounds** (#65). `tight_bounds` gives the smallest axis-aligned
  box holding a shape, each side where the shape actually reaches: a
  revolved tube of outer radius 6 bounds at 6, not at its carrier's 8.29.

### Fixed

- **A washer's medial axis did not settle (#71).** A face bounded by two
  full circles has no branch point, and its one branch, the ellipse
  between the circles, was never emitted. It now comes back as one closed
  branch naming a single vertex at both ends.
- **A fillet on a prism's base edge was refused as too wide (#70).** The
  check stepped into the reversed base face on the wrong side of the
  edge. Each face's inner side is now read from the face itself.
- **A rim blend beside a sphere's pole could not be built.** Where the
  ball's rail on the sphere passes over the pole while the rim does not,
  the sphere's leg was asked to be a band between two loops that both
  wind round the pole. It is now the band from the rail to the pole with
  the rim cut from it.
- **A section through a point where two surfaces touch was fitted loosely.**
  The tracer could report one sample's position with its neighbour's
  parameters there, and the joint fit stalled near 1e-5 mm at that
  sample. Such samples are dropped before fitting, and a torus against a
  cylinder now fits to about 3e-7 mm.
- **A straight band after a rim arc's band could not close at a curved
  corner.** Its section through the rim's torus runs tangent to the
  torus's end meridian into the corner, stopped short of the vertex, and
  split the meridian into a sliver. A tangent touch at an edge's end is
  now that end's vertex, and the section's chart image is bent onto it.
- **A flush fillet beside a rounded apex could fail to close.** A straight
  band end tangent to the corner sphere's rim was read as running along
  the rim for a short stretch, which split the rim at a vertex the band
  never had. A line and a circle now share no stretch. The four flush
  fillets close after an oblique apex's corner, and the five after an
  irregular pentagonal one.
- **A holed profile round a closed spine filled its hole.** Where the
  hole's ring wound the other way from the one the sweep assumed, its
  tunnel came out as material and the solid measured the hole's volume
  too many. Each shell is now turned by the volume it encloses.
- **A partial revolve could build its end caps inside out.** Where the
  profile's ring wound against its plane's normal, the walls followed the
  ring and the caps the plane, and the caps faced into the solid; a cap
  in a plane through the axis adds nothing to the volume, so only its
  orientation was wrong. The caps now face the way the prism's do and the
  walls are turned to meet them. Surfaces of revolution of a spline
  curve trimmed by their chart now measure exactly too.
- **Spline faces with any trim measure exactly.** Mass properties
  integrate a spline face, or a spline curve extruded or revolved, round
  its own trims on the exact surface, the panels broken at its knots; a
  converted box drilled through keeps its volume to rounding.
- **A helical sweep measured 0.65% light at the default deflection**
  (#69). Its fitted walls were measured from their mesh. A spline face
  trimmed by its own chart's borders is now integrated on the exact
  surface, knot span by knot span, and the thread measures to Pappus.
- **A chamfer wider than its face cut through it** (#67): a 12 mm bevel
  on a 10 mm face sliced a slab off the box. A setback running past a
  face's far side is now refused, naming the edge and the distance. A
  fillet is held to the same on any edge, curved faces and concave edges
  included: its ball's contacts, set back the radius times the tangent of
  half the faces' turn, must land on the solid's boundary.
- **A refined solid of pads measured from its mesh** (#64). The face
  `unify_same_domain` merged lay on one of the merged planes, whose window
  did not cover the union, and the plane refused to be read past it. A
  plane is now read anywhere; its window only says where its face was
  built.
- **A wall drafted a degree and a half could not be measured** (#66): its
  trim reached past the plane's window. The same change reads it.
- **A pad placed again at a scale measured and fused wrong** (#63): a
  cube padded from a square and copied 1.5 times the size measured 144
  where it is 216. An edge's trim stored under the pad's own placement
  was not found once the whole solid was placed again, and the trim of
  the edge's unplaced twin was read instead.
- **A partial revolve of a profile touching its axis failed** (#62). The
  edge on the axis is its own image; a wire now joins one vertex under two
  placements that land on one point, so the start and end faces close. A
  profile face built without trims gets exact ones before a prism or
  revolution keeps it as an end, and is read at its own middle rather than
  at its plane's origin.
- **A pipe shell failed on a profile ring walked by reversed edges**
  (#58). A reversed edge's ends were swapped twice.
- **A profile square to a curved spine's start was refused as leaning**
  (#59). The test now reads the spine's exact start tangent, to an
  angle's tolerance.
- **A circle piped along a line was 2.6% short of its cylinder** (#60),
  and a square along an arc 0.6% short of Pappus. Both are now exact.
- **A thread groove failed to cut from a bored block at some lengths.**
  Where the sweep had turned a flank just inside the bore, the flank met
  the bore's wall along a short run at a grazing angle, shallower than the
  sampled cells could see, and that section was missed. Sections are now
  also found from where they cross a spline's border, so 5.75 and 6.75
  turns cut.
- **Rebuilding a torus broke its trims.** `to_nurbs` and `baked_shape`
  took every seam for a `u` column, so a torus's `v` seam was given
  pcurves on the wrong edges of the chart and the face could not be
  meshed. A seam on the `v` rows now gets its trims there.
- **A pipe shell failed along an arc that runs on into a straight line.**
  One skin fitted across the join, where the curvature steps, and could not
  reach the tolerance. Each spine edge now skins its own run, the runs
  sharing the section at the join, and a run whose sections share a plane
  is built on that plane.
- **IGES curve-on-surface entities read their model-space curve from the
  wrong field.** The reader took the parameter-space curve as the
  model-space one; it now reads the fields in the order the format gives
  them.
- **A boolean on a prism split its near cap where its far cap was cut.**
  The far cap of an extrusion is its near cap moved, sharing its edges
  under a displacement, and the boolean shared edge splits by edge alone.
  Edges now count as one only at the same placement.
- **`analyse_blend` failed on a corner patch whose corner is its ball's
  pole.** A pcurve ending on the pole stood a rounding past the sphere's
  chart there; the station is now read at the chart's edge.

## [0.3.3] - 2026-09-25

A patch release with no API change; `cargo semver-checks` finds nothing
against 0.3.2. `make_half_space` accepts any face, so a cut or common can
stop on a drum, ball, cone, ring or spline patch. A thread groove cuts
from a plain or bored block. Holes drilled through large converted parts,
beside rounded corners, near a ball's rim and all but parallel to a bore
cut valid, and booleans on parts of tens of thousands of faces spend far
less time on bookkeeping.

### Added

- **Half spaces bounded by curved faces.** `cut`, `common` and `section`
  accept a half space from `make_half_space` whatever its face's surface.
  A drum, cone, ball or ring divides all of space and is exact on its
  surface; a drum and a cone are read as unbounded along their axes. An
  open surface (a spline patch, an extrusion, a revolution, trimmed or
  offset surfaces) divides space across its own extent: the other
  argument must lie across it, and the surface must not fold back along
  the direction into its material, or the operation says which. An
  extrusion can stop exactly on a curved face.

### Fixed

- **A thread groove would not cut from a block.** A section swept along a
  helix skins one wall per profile edge, each fitted at its own pace
  along the sweep, and a wall adopting its neighbour's rail took that
  rail's image on itself to be straight; along a helix it drifted a fifth
  of a millimetre off. The image is now read off the adopting wall. And a
  ray starting a few microns from a wall winding several turns crossed it
  in the gap between the wall and the flat cells seeding its crossings,
  so a point in the groove read as outside; a curve near a cell by its
  surface's bow from it is now seeded there. The groove cuts from a block
  for any number of turns.
- **A thread groove would not cut from a bored block.** Where the groove
  runs in and out of the bore, four things went wrong. A section winding
  several turns round the bore was fitted as one curve that wandered
  hundreds of millimetres between the points it was checked at; fits are
  now checked between their samples too, and a trace no single fit
  follows is fitted in pieces. Such a section crossed the bore's seam
  once a turn and was read as running along it. A section ending along a
  rail, or two meeting a rail that grazes the bore, ended beside the
  rail's split point by more than the weld reached in the face's chart;
  they now meet where the junction there says. And a ray crossing a wall
  near the end of its patch was solved from a seed that stepped off the
  patch, so a point in the groove's end read as outside. The groove cuts
  from a bored block for one to seven turns.
- **`make_half_space` could pick the wrong side of a closed surface.** The
  side was read from the normal at one sample of the face, which on a rod
  or a ball can face the given point across the surface. It is read where
  the surface comes nearest the point.

- **A hole drilled through a large converted part failed, and slowly.**
  One loose edge on a face set how far every end on that face was welded,
  so short edges and short sections collapsed whole. On a small sphere the
  weld was also taken in radians as if they were millimetres. The weld
  now reaches only as far as each piece's own doubt, scaled into the
  face's chart. Where a circle or ellipse crosses another in a different
  plane, the crossing is solved in closed form rather than by sampling.
  On a converted part of 25 000 faces a drill cuts in 6 s where it
  failed after 18 s.
- **A hole drilled beside a rounded corner often failed.** Where a hole
  crosses a fillet's tangent line, the section on the flat and the section
  on the fillet meet tangentially, and their crossing with each other or
  with the line landed up to a couple of microns off. The flat was left
  unsplit or a wire had two vertices where it needed one. Such a crossing
  now stops where the sections cross the shared edge, and a face welds its
  ends as far apart as the stops on its edges were found. A seeded scan of
  a hundred holes beside the corner of a rounded box went from a third
  failing to none.
- **A converted part's volume could not be measured.** Its faces meet on
  curves threaded through the mesh a few hundredths off, and the drawn
  mesh did not weld shut there, so the volume was refused as having no
  closed boundary. Closure is now asked of the topology, and where the
  drawn mesh stays open the volume is summed face by face from each
  face's own mesh. Cracks in a drawn mesh narrower than its chord are
  sealed, and triangles laid twice facing opposite ways are cancelled.
- **Booleans on parts of tens of thousands of faces spent seconds on
  bookkeeping.** Junctions were merged by comparing every pair, deleted
  faces were found by scanning every kept one, sewing projected edge
  middles before checking the ends matched, and a point classified
  against a solid asked every face. Each now looks only at what is near
  it. A line crossing a circle or ellipse is solved in closed form rather
  than by sampling. A drill through a converted part of 25 000 faces
  cuts in 2 s.
- **Sewing could open a face's wire.** Two edges within their own
  tolerance of each other were merged, but when their ends were two
  vertices farther apart than either vertex's tolerance, the faces
  around the dropped edge kept ending on its vertex, and their wires
  broke. Merged edges now end on one vertex, widened to reach both.
- **Holes through a converted part failed or took the wrong volume.**
  Two causes, together half of a scan of random holes through a real
  printer's part. A plane all but parallel to a hole's axis meets it in
  an ellipse kilometres long, whose short stretch inside the hole was
  missed by sampling, so the face was never cut. And a bore the converter
  opens along a slit, its chart starting at the slit rather than at the
  cylinder's zero, had points just short of the slit read as outside it:
  sections were cut short, a ray missed a wall and read material as air,
  and holes were attached a period away from the face they cut. Points
  are now placed on the side of the seam the face's trim is on, and a
  section is kept wherever it may pass through the surfaces' extents.
- **A hole's rim passing close to a loose edge was dropped.** Splitting a
  face, a loop was taken for part of the face's boundary, not a hole in
  it, whenever any point of it came within the face's weld of the
  boundary, and the weld on a face with loose edges is a few hundredths.
  A loop is part of the boundary only when it meets it at a node.
- **Circles and ellipses in all but one plane missed their crossing.**
  The closed form for conics in different planes looked for crossings
  only where one crosses the other's plane. Planes a hair apart, a hole's
  rim fitted to one facet group and a section through another, let the
  curves pass within the gap elsewhere, and the crossing was missed. Such
  pairs are measured by sampling.
- **A hole's rim cut a few microns from its own start stayed open.** A
  rim is cut where it crosses the hole's seam, and when that falls a few
  microns from where the rim's curve starts, the sliver between collapses
  in the chart of a face with loose edges. It did not collapse in space,
  so the rest of the rim ended on two vertices. A strand a face's weld
  collapses now makes its ends one junction.
- **Holes all but parallel to a bore or wall of a converted part failed,
  and slowly.** Two drums whose axes lean a ten-thousandth apart were
  marched, and their fitted crossings cost seconds each to intersect; a
  plane leaning a hundred-thousandth off a drum's axis met it in an
  ellipse tens of metres long, too coarse a ruler for the drum's height.
  Both are now solved in the drum's cross-sections along its height and
  kept as a line, or a cubic through the stations, with their departure
  stated as the section's tolerance. A hole that took 17 s and failed
  cuts in under a second, and a scan of a hundred random holes through a
  converted printer's part fails in one, from four.
- **A mesh with a flat sliver at a T-junction converted open.** A
  triangle with three distinct corners and no area was dropped, and the
  long side it sealed and its two short sides were each left used once.
  The triangle across the long side is now split at the sliver's middle
  corner, and the solid closes.
- **A hole near a ball's rim failed.** A drum passing through a ball
  close to its rim, all but grazing its far side, was marched, and the
  trace wandered along the long thin loops into a curve kilometres long.
  Where every line along the drum meets the ball twice, each loop is a
  function of the angle round the drum; the loops are now sampled exactly
  so, fitted with their images on both surfaces, and checked midway.

## [0.3.2] - 2026-09-24

A patch release with no API change; `cargo semver-checks` finds nothing
against 0.3.1. `project_edge_onto_plane` is new. Meshes convert far more
completely: coarse fillets, curved faces meeting along any curve, balls
and rings with holes, and a real printer's part all come back exact and
valid. Mass properties integrate most faces exactly whatever their trim,
and booleans cut cleanly along a closed surface's seam and through a
countersink whose bore is its narrow end.

### Added

- **`project_edge_onto_plane`** projects an edge orthogonally onto a plane
  as an exact curve in the plane's coordinates, for a sketch to constrain
  against. A line comes back a line or a point; a circle or ellipse a
  circle, an ellipse or a segment, its arc range counter-clockwise; a
  B-spline the spline on its projected control points. Any other curve
  is fitted and the fit's error returned.

### Fixed

- **The rim of a bore cut with a prism would not fillet or chamfer.** A
  circle face pushed down through a block, against its own axis, swept a
  general extrusion rather than the cylinder it is, and the blend refused
  its rim. It now sweeps a cylinder whichever way along its axis it
  travels, and the rim rounds with one torus, or bevels with one cone,
  exactly.

- **A converted printer's part came out invalid, a chamfer and fillet
  corners faceted.** Edge tolerances were recorded equal to the gaps they
  were measured from, and the checker, measuring again, found them a
  rounding over; they are recorded a millionth wider. A band round its
  axis cut open along a slit (a few triangles off its surface, left by
  the exporter) is built as a patch running round to the slit instead of
  falling to facets. And where two faces meet all but tangentially, a
  fillet running into a corner the mesh leaves free-form, the boundary is
  a curve threaded through the chain's own vertices, its tolerance how
  far it strays, instead of faceting the fillet whole.

- **Balls cut by two planes, and bores crossing at a coarse mesh, came
  back wrong or faceted.** A ball's axis could put both poles inside the
  cut-away part, where no loop goes round them, and the face covered the
  wrong side; both poles now stand on the ball's own surface. A flat face
  with two edges was taken for a sliver strip whenever it had two, and is
  now only when both are straight. An edge with an exact curve gets its
  image on a curved face interpolated at the edge's own parameters, close
  enough for its volume to be integrated exactly. Where a coarse mesh
  leaves skewed triangles on a surface next to a hole, or a few slivers
  together round a pole, they are taken into the surface when their
  middles sag no more than its curvature explains.

- **A countersink whose bore was exactly its narrow end refused to cut.**
  The cone's rim lay in the bore's wall, and the rim edge and the section
  along it were split at different points: the section at its middle, the
  rim where another section crossed it. Once every section is paved, a
  piece running along an edge is cut wherever the edge is split, and the
  edge wherever the piece ends.

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

[Unreleased]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.8.0...HEAD
[0.8.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.6.1...v0.7.0
[0.6.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.5.1...v0.6.0
[0.5.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.4...v0.4.0
[0.3.4]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.3...v0.3.4
[0.3.3]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.2...v0.3.3
[0.3.2]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/gilbertorconde/ogeom-rs/releases/tag/v0.1.0
