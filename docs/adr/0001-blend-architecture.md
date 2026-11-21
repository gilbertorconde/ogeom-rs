# ADR 0001: One blend model for fillets and chamfers

- **Status:** proposed (#902). The behavioural questions are settled by survey (§13).
- **Scope:** `ogeom-fillet`: edge rounds, chamfers, asymmetric and variable
  sections, vertex blends. Face–face blends (`blend_faces`) and sheet rounds
  follow later and are out of this record's first steps.

## 1. Context

`ogeom-fillet` is 40 files and about 51,600 lines, with about 650 refusal
sites. It has had 187 commits since 2026-10-01, 77 of them `fix(fillet)`.

- **Wrong solids:** in the last two release cycles it has built wrong solids
  that `check` passed: #713/#839, #846, #868/#871, #843.
- **Order and radius dependence:** results depend on the order the edges are
  asked in, and on the radius (#838, #809, #878).
- **Speed:** it is slow (#870).
- **Refusals:** it refuses simple parts. The stepped bracket in #902 is one:
  chamfering every edge gives a concave bevel hanging past both side faces,
  unbevelled side profiles and inconsistent vertex facets, and filleting
  every edge doesn't build.

What the code does today. Each fact is checked against `main`; files are
under `crates/ogeom-fillet/src/`.

- **Edge by edge, against a changing solid.**
  - `fillet_edges_own` (`fillet.rs`) builds the bands of edges that share a
    vertex one after another, each against what the last one left. Its
    comment says so: the order the edges were asked in decides everything
    else.
  - Chamfers apply one boolean per wedge in request order
    (`chamfer_edges_with_own`, `chamfer.rs`).
  - The request is sorted into a canonical order only when a signed or cusp
    corner exists (`fillet_edges_meeting`).
- **Ends and corners as a catalogue.** Each configuration has its own
  construction, its own refusals and its own cell cut and re-sew:
  - corners: setback, mixed and signed (`setback.rs`), cusp (`cusp.rs`), one
    edge left sharp (`sharp.rs`), ball-and-block (`corner.rs`),
    chamfer-facet (`chamfer_corner.rs`);
  - junctions where creases turn (`junction.rs`);
  - a ball on a far edge (`rolled.rs`), a staircase (`stepped.rs`), a contact
    changing kind along the edge (`carried.rs`), a strip taken off whole
    (`consumed*.rs`), a section shrinking to a pole (`pinched.rs`);
  - cut-off-flush ends (`cut_off_bands`, `fillet.rs`).

  A configuration outside the catalogue is refused, or built wrong where a
  case is misapplied.
- **Many rebuilds per request.**
  - A request costs one corner boolean, at least two booleans per "round" of
    vertex-disjoint edges, about ten per junction corner, one cut and one
    full re-sew per special corner, and two booleans per cut-off.
  - Retry ladders multiply this by up to six: exact rims, then fitted rims;
    crossed junctions, then flush, then none.
- **Validation after the fact.** Corner patches are fitted by iterative net
  refinement, up to 135 controls, with hold passes and 202-control nets
  (`fill_n.rs`). Each is then measured against a region and refused if it
  strays. Soundness rests on those checks being complete, and they have not
  been.
- **Planar chamfers share the curved machinery.**
  - Planar chamfers use the shared wedge-and-boolean scaffolding.
  - Where a bevel meets a round, they go through a cell cut, an N-sided G1
    patch and a re-sew (`blend.rs`, `CornerEnd::Bevel`).
  - A chamfer of a polyhedral part is an exactly solvable plane problem, yet
    it can be slow, flaky or refused.

## 2. Decision

Replace the case catalogue with **one blend model**, built in fixed phases
from the **original solid and the whole request at once**:

```
request ─► P0 normalise & classify ─► P1 sections & contact curves (per edge, on the original solid)
        ─► P2 setbacks & ends (per face corner, by intersection) ─► P3 vertex blends (one method)
        ─► P4 one assembly (trimmed faces + strips + vertex patches, sewn by known adjacency)
        ─► P5 independent guards (check, region, overhang)
```

The model rests on one topological fact, the **face-setback structure**,
shared by every blend kind:

- Every original face shrinks to the region inside its **setback curves**:
  the contact curves of the blends on its edges, or the edge itself where
  the edge is not blended.
- Every blended edge becomes a **strip** between its contact curves on its
  two faces.
- Every vertex where blended edges meet becomes a **vertex blend** bounded
  by the strips' end sections and the faces' setback corners.

This is the polyhedral bevel and the setback vertex blend of Várady and
Rockwood [5, 6], made the backbone of the whole operation. The topology of
the result is known from the request alone, before any geometry is
computed. That pseudo-topology is the same idea as the SMLib fillet
executive's "fillet brep" built ahead of geometry [14].

**Exact polyhedral chamfers are the first phase delivered.** They prove the
backbone at no geometric risk (§5).

## 3. Semantics: what a request means

A request is a set of `(edge, section)` pairs. A section is one of:

- a circular round of radius `r`;
- a conic round with distances `d1`, `d2` and a shape parameter (#877);
- a chamfer with distances `d1`, `d2`, or a distance and an angle;
- any of these varying along the edge under a law.

The result is defined as a function of the set, not of a sequence, as
follows.

1. **Contact curves.**
   - A blend's contact curve on face `F` is where its section touches `F`:
     - for a round, the rolling ball's contact, so the spine is where the
       two faces' offsets meet [2, 4];
     - for a chamfer, the line or curve at the given distance from the edge,
       measured on `F`.
   - On a plane beside a straight edge it is a line parallel to the edge.
     On curved faces it is the marched contact the crate already computes.
2. **Setback corners.**
   - Where two edges of `F` meet at a vertex `v`, the face's setback corner
     `S(F, v)` is where their contact curves on `F` meet.
   - An unblended edge contributes the edge itself, at distance 0.
3. **Strips** run between `S(F1, v)` and `S(F1, w)` on one side and
   `S(F2, v)` and `S(F2, w)` on the other, for an edge `vw` between faces
   `F1` and `F2`. Their ends are where the section sweep meets the vertex
   blend.
4. **Vertex blends** fill the polygon of setback corners around `v`, by the
   rules settled in §13:
   - **Chamfers:** a planar facet wherever the setback corners are
     coplanar, always so where three bevels meet. Otherwise one patch
     through them, not a triangulation.
   - **Rounds:** a sphere where an equal-radius convex or concave trihedral
     corner allows one. Otherwise one n-sided G1 patch (§4.4).
   - **No setback** beyond the bands' own ends unless the request asks for
     one, with its distances.
5. **Admissibility.** A request fits when:
   - every face's setback curves leave it a simple region, so no setback
     curves cross and no shrunk face inverts;
   - every strip's contact curves stay on their faces;
   - every vertex polygon is non-degenerate.

   These predicates are the domain. `largest_fillet_radius` and
   `largest_chamfer_distance` become their exact boundary, not a separate
   search. Outside the domain the request is refused by name, naming the
   face or vertex that collapses.
6. **Chamfers and rounds in one request** compose as two features: the
   chamfers first, as one set, then the rounds, as one set, on that result
   (§13, question 3). Each kind is order-free within itself, and the order
   of the kinds is fixed, so the request's order never matters.

Rounds have a global reference semantics: rounding every convex edge is the
morphological opening by a ball of radius `r`, and filleting every concave
edge is the closing [3]. This is the **oracle** for tests on pure-convex and
pure-concave regions. It is not the construction. It has no meaning for
mixed vertices or partial requests, and building it is far costlier than the
local model.

## 4. Architecture

### 4.1 P0: normalise and classify

- **Chains.** Group tangent-continuous edges into chains, so a chain is one
  strip with one section law.
- **Vertex classes.** Classify every vertex as `N × M` (N edges, M blended),
  with the count of convex, concave and tangent edges. This is the
  classification commercial systems publish for vertex blends [7, 14]. It
  selects the vertex method in P3 and nothing else; no geometry is decided
  by class.
- **Stable keys.** Give every element a key derived from the geometry and
  topology of the input, not from request order. All iteration below runs in
  key order, which makes the result bit-identical under any permutation of
  the request.

### 4.2 P1: sections and contact curves, per chain, on the original solid

- Each chain's section is solved against the **original** faces. Nothing
  built in P1 changes the solid another chain sees, so order independence
  holds by construction, not by sorting.
- **Exact where a closed form exists.** Plane–plane is a cylinder (round) or
  a plane (chamfer); plane–cylinder and plane–cone coaxial are a torus; and
  the other closed forms in `support.rs` (`planar_seat`, `revolved_seat`),
  `fillet.rs` (`section_at`) and `asym.rs` (`seated_conic`, `conic_pieces`).
- **Marched otherwise.** `march.rs` (`march_blend`, `march_blend_law`,
  `seat_section`) and `pair_marched.rs` (`pair_seat`) already take two
  surfaces, a section law and a guide, and return the spine and contacts in
  each support's own parameters. They need no topology.
- **Strips are extended.** Each strip is computed past its edge's ends far
  enough to reach its vertex blends, as the SMLib sheets are made "way long"
  [14]. Its true ends come from P2, not from the section solver.

### 4.3 P2: setbacks and ends, by intersection

- On each face, intersect neighbouring contact curves to get the setback
  corners `S(F, v)`. On planes this is a 3×3 linear solve. On curved faces
  it is a curve–curve intersection in the face's parameter space.
- A blend ends where its contact curve meets:
  - another blend's contact curve on the same face (a setback corner);
  - the face's boundary (an unblended neighbour edge);
  - a face it was not seated on (roll-over: the contact leaves its face).

  This is one intersection query per end, not a case. The kind of end only
  labels the result.
- **Roll-over** is first refused by name, as today, then handled in a later
  phase. The strip continues on the new face with the contact re-seated
  there, the cliff-edge and roll-over geometry of [14, 15].
- **Admissibility predicates (§3.5)** are evaluated here, on curves and
  corners, before any surface is built. A request that fails them is refused
  with the colliding pair named. Nothing is built and then measured.

### 4.4 P3: vertex blends, one method with stated exact cases

Every vertex blend is a hole bounded by the strips' end sections and the
faces' setback curves near `v`. Its boundary is closed and its corner data
consistent by construction, because each corner is a setback corner shared
by exactly the two curves that meet there.

1. **Exact cases, chosen by geometry, not by catalogue:**
   - **Planar facet:** chamfers whose setback corners are coplanar, always
     so where three bevels meet.
   - **Sphere:** equal-radius rounds at a trihedral vertex where one ball
     touches all three faces. That is today's ball-and-block corner, kept as
     a special value of the general method.
   - **Torus or plane continuation:** a tangent `2 × 2` vertex, where one
     strip continues the other.
2. **General case:** one n-sided G1 patch interpolating the strips' end
   sections and the faces along the setback curves.
   - The method is a ribbon-based transfinite patch [8, 9], or Gregory's
     n-sided patch [10, 11] for n ≤ 6.
   - It is evaluated, or fitted once at a fixed resolution with the residual
     reported. It is not refined until a check passes.
   - Its hole is the strips' natural ends: each strip ends on its surface
     through its two faces' setback corners (skew where those stand at
     different stations), so the hole's corners are the setback corners.
   - The patch is Kato's blend of the sides' ribbons over the triangle.
     Each ribbon leaves its side along a cross tangent in the strip's
     tangent plane, its ends fixed by the corners and its slopes there
     making the two ribbons' twists agree in the face both strips touch,
     and bends as the strip does, so that near its side it keeps to the
     strip's own continuation. It takes each side's tangent plane exactly,
     however the strips' curvatures at a corner disagree (see the amendment
     below).
   - The face carries it as a clamped bicubic B-spline taking the exact
     patch at a fixed net of points and its cross-boundary derivatives
     along the sides, closest at the corners. Its gap and tangency to every
     strip, and its region, are measured on that B-spline, never assumed
     from the exact patch.

   Setbacks are the knob Várady and Rockwood give for patch quality: larger
   setbacks give rounder patches and fewer bad corner angles [5, 6]. They
   are applied only where the request asks for them, with its distances
   (§13, question 2).
3. **Chamfers at vertices of more than three bevels:**
   - The facet is planar where the setback corners are coplanar.
   - Otherwise the corner is one patch through them, never a
     triangulation (§13, question 1). It is G0, not the general case's G1
     patch (see the amendment below): the bilinear patch at four corners,
     the corners' Wachspress blend at five or more.

### 4.5 P4: one assembly

- Every original face is trimmed once, along its setback curves.
- The strips and vertex blends are built as faces whose boundary edges are
  the setback curves and section curves computed in P2.
- The shell is assembled by **known adjacency** from the P0 pseudo-topology:
  each new edge is created once and shared by its two faces. No tolerance
  sewing searches for matches, and no boolean runs on the local path.
- The one boolean that remains is for **non-local** blends [12]: a strip
  that, extended, crosses a face that is not its neighbour (a fillet
  reaching a boss). The strip set is closed into a blend solid and one
  boolean is applied (the "sheets then one boolean" design, §9 B). The
  first steps refuse it by name, as non-local; it is a later phase.

### 4.6 P5: independent guards

These are kept, and kept independent of the construction:

- `check`, and orientation consistency of every new face.
- **Region:**
  - every convex strip and vertex blend lies inside the original solid;
  - every concave one lies inside the solid's closing at the largest radius
    asked, within the fill tolerance;
  - every vertex blend lies inside the convex hull of its boundary
    ribbons' region.

  These bounds come from the semantics in §3, not from per-case cells.
- **No overhang:** no new face looks down past its strips' contact normals
  by more than the fill wobble.
- `analyse_blend` on every strip and vertex blend: gap, tangency within
  2e-4 rad, curvature where G2 is asked.

A guard that fails is an engine bug and refuses the request. It is never the
way the engine searches for a working construction.

The region guard proves faces, not samples:
`ogeom_algo::certify_face_within` bounds every face over rectangles of its
parameters by the convex hull of its exact B-spline form there, splitting
what it cannot clear, against the region as half-spaces under maxima and
minima. Its one input from the mesher is the trim, read from the face's
chart rings drawn within a chord.

## 5. Polyhedral fast path: exact chamfers (prototype result)

A throwaway prototype (Python, not committed) built the §3 structure for
planar solids. Per face and vertex it solves for `S(F, v)` (three planes),
then builds:

- shrunk faces;
- edge quads, which are planar because the two setback lines are parallel to
  the edge;
- vertex polygons.

It checked the result: closed (volume independent of the origin), no two
triangles crossing, and exact volumes. On a stepped bracket like #902's
(profile 60 deep, plate 6 thick, step to 16, block to 30, 80 long; 16
vertices, 24 edges, all 3-valent):

| request | result | time (Python) |
|---|---|---|
| every edge, d = 1, 2, 3 | valid, closed, 0 crossings; every edge bevelled, side profiles included; a triangular facet at every vertex, the mixed ones included | ≈ 4 ms |
| every edge, d = 4.9 | 13 crossings: past the domain, since the 6-thick plate admits d ≤ 3, and detected | ≈ 4 ms |
| the concave step edge alone, d = 2 | ends flush with both side faces; volume 69,760 = 69,600 + d²/2 × 80, exactly | ≈ 3 ms |
| asymmetric distances 1 / 2.5 | valid, closed | ≈ 4 ms |
| five random halves of the edges | all valid and closed | ≈ 3.5 ms each |

- **Order independence.** The construction reads a function of
  `(edge, face)`, never a sequence, so every permutation of the request
  gives the same solid.
- **No booleans and no fitting.** Every coordinate is a 3×3 solve.
- **The domain** is the face-polygon offset staying simple plus the facets
  staying non-degenerate. Both are exact predicates, and the largest distance
  is computable in closed form per face.

Planar-faced requests whose sections are all chamfers take this path
outright. That covers the #902 bracket's chamfer-all-edges case and every
chamfer of a prismatic part.

**Rounds on polyhedral parts** use the same backbone with exact geometry:

- **Strips:** cylinders, between tangent lines at `r · cot(α/2)` from the
  edge for an interior dihedral angle `α`.
- **Vertices:** spheres where one ball fits; otherwise the §4.4 patch.

Their oracle is the closed-form section area along each strip, plus the
opening and closing reference [3] on pure-convex or pure-concave parts.

## 6. Soundness by construction

| Failure seen | How the design rules it out |
|---|---|
| Order and radius dependence (#838, #809, #878) | P1 reads only the original solid. Every phase iterates in stable-key order. No phase depends on another chain's output except through P2's symmetric intersections. |
| Blend past a face's side over air (#868, #871) | A contact curve that leaves its face is an end found in P2. The strip is trimmed there, or roll-over is refused, before any surface exists. There is no "run off and cut later". |
| Corner patch outside the solid (#713, #839, #846) | The vertex hole's boundary and region are fixed before the patch. Exact cases need no fit. The general patch is evaluated or fitted once against known corner data, and P5 measures it independently. |
| Feature lost inside a corner cell (#843) | No cells are cut: faces are trimmed along curves, so nothing inside a box is discarded. |
| Concave bevel past both side faces, unbevelled profile edges, inconsistent facets (#902) | §5: the face-setback structure gives every edge its quad and every vertex its facet by one rule. |
| Fill nets refined to 202 controls, then checked (#870) | Exact cases need no fitting. The general patch has one fixed resolution with its residual reported. |

## 7. Complexity targets

- **P0–P2:** `O(E + V)` local work, plus one march per curved chain. The
  march's cost is the solver's, unchanged.
- **P3:** constant work per exact vertex. One patch per general vertex, at a
  fixed resolution.
- **P4:** one trim per face and one assembly: `O(F + E + V)`.

Targets, release build:

| Case | Target |
|---|---|
| Chamfer every edge of a 50-edge prismatic part | under 10 ms |
| Round every edge of a 50-edge prismatic part | under 1 s, as #902 asks; analytic strips and sphere or single-patch corners make well under 100 ms reachable |
| Curved chains | bounded by the march; no retry ladders and no whole-request reruns |

## 8. How #877 and variable radii fit

- **The section is the only per-edge law.**
  - #877's asymmetric distances and conic profiles (`asym.rs`,
    `asym_law.rs`, merged through #898) become section laws in P1. A conic
    section's contact curves sit at `d1` and `d2` on the two faces.
  - Asymmetric chamfers are a chamfer section with two distances. The
    prototype already builds them on the polyhedral path.
  - Variable radius is a section whose parameters follow `RadiusLaw`
    (`law.rs`) along the chain.
- **P2–P5 never look inside a section.** Setback corners are where contact
  curves meet, whatever made the curves. The vertex patch interpolates
  whatever ribbons the strips present.
- **One open semantic question:** what a chamfer meeting a round at one
  vertex should look like (§13).
- **Coordination:** the #729 agent's #877 work keeps landing on the current
  path. Its section math (`seated_conic`, `conic_pieces`, the law) moves
  into P1 unchanged. Its end and corner handling (`after_symmetric`, corner
  routing in `blend.rs`) is what retires.

## 9. Alternatives weighed

- **A. Keep extending the catalogue.**
  - Rejected. Each case adds constructions, refusals and checks. The wrong
    solids came from cases misapplied or checks incomplete, and order
    dependence is inherent to sequential application.
- **B. Blend sheets, then one boolean or merge** (the SMLib executive's
  global merge [14]).
  - Strong on generality: non-local blends fall out.
  - Booleans on blend sheets meet near-tangent faces at every contact curve,
    which is the hardest input a boolean has. The result's topology is
    discovered, not known.
  - Kept as the P4 fallback for non-local blends only.
- **C. Morphological offsets** (opening and closing [3]).
  - Exact rolling-ball semantics for whole regions.
  - Undefined for mixed vertices, partial requests, chamfers and asymmetric
    sections. Needs exact solid offsetting, which is costlier than the
    blends themselves.
  - Kept as the test oracle.
- **D. Implicit or potential blends** [16].
  - Smooth vertex blends for free, but they don't produce the trimmed
    parametric faces a b-rep needs. Rejected.
- **E. Direct topology surgery per edge** (Euler operators, edge by edge).
  - Same fragility as A: each step changes what the next sees.
  - The adopted design is surgery done once, from the precomputed
    pseudo-topology.

## 10. What is reused and what retires

Classes from a survey of `main`.

**Reused unchanged or nearly:**

| File or function | Role |
|---|---|
| `march.rs` (`march_blend`, `march_blend_law`, `seat_section`) | Marched sections; topology-free |
| `pair_marched.rs` (`pair_seat`) | Guide where the surfaces or their offsets cross |
| `support.rs` (`planar_seat`, `revolved_seat`), `fillet.rs` (`section_at`), `ruled.rs`, `asym.rs` (`seated_conic`, `conic_pieces`), `law.rs` | Closed-form sections and section laws |
| `marched.rs` (`conic_section`, `fit_open_band`, `fit_closed_band`) | Band fitting. Not the wedge or run-out builders. |
| `fill_n.rs` (`make_filling_n_region`) | Fitter for general vertex patches, fed consistent corner data |
| `fill_patch.rs`, `bridge.rs` | Hermite blend surfaces |
| `between.rs`, `support.rs` (`ball_fits_by`, `bands_clear`) | Admissibility predicates, read on the original solid |
| `fits.rs` | Domain boundary, reworked onto the §3.5 predicates |
| `analyse.rs` | P5 measurement |
| `setback.rs` (`certify_within`, `least_setbacks`) | Closed-form corner maths |
| `corner_curved.rs`, `turned_pole.rs`, `chamfer_marched.rs` | Section and chart utilities |

**Retired in steps:**

| File | Case |
|---|---|
| `setback.rs` (corner constructions), `cusp.rs`, `sharp.rs`, `corner.rs` (the corner booleans, not its frame helpers), `chamfer_corner.rs` | Corner catalogue |
| `junction.rs`, `carried.rs`, `rolled.rs`, `stepped.rs`, `consumed*.rs`, `pinched.rs` | End and junction catalogue |
| `fillet.rs` (`cut_off_bands`, the round/retry orchestration), `chamfer.rs` (sequential wedges), `blend.rs` | Sequential orchestration |

What these do correctly today becomes corpus cases for the new path.

## 11. Migration plan

Every step keeps the current public API. A step routes a class of request
to the new engine only when that class's corpus passes there. Everything
else stays on the current path.

- **M1: polyhedral chamfers.**
  - New engine for requests whose faces are all planar and whose sections
    are all chamfers, symmetric or asymmetric.
  - Oracles: exact volume (divergence over the built faces, against an
    independent plane-arrangement volume), closedness and no crossing, and
    every permutation of the request bit-identical.
  - Ships with the #902 bracket chamfer case.
  - The old chamfer path stays for curved faces.
- **M2: request model and P0–P2 for planar faces with rounds.** Exact
  cylinder strips. Sphere and planar exact vertices. The general vertex
  patch of M4, brought forward (see the amendment below).
  - Corpus: box, L, T, U, stepped block, the #902 bracket fillet-all-edges,
    mixed rounds and chamfers.
- **M3: curved faces.** P1 on `march.rs` and `pair_marched.rs`. P2
  intersections in face parameter space. The plate with a boss, tangent and
  not (the #868 and #878 families).
- **M4: the general vertex patch** (ribbon-based transfinite or Gregory),
  replacing iterative fill refinement for vertex holes. Its three-sided
  case for rounds lands before the M2 switch; more sides and other
  sections follow here.
- **M5: roll-over and non-local blends.** Roll-over re-seats the contact on
  the new face. Non-local blends go through the one-boolean fallback.
- **M6: retire old paths,** file by file, once no routed class reaches them.
  The real-world check (Voron 2.4, DoomCube) runs at every routing change.

**"Keeps every currently correct result"** means that every case the current
tests accept passes on the new path within its oracle. Bit-identity to the
old path is not promised where the old construction was a fitted patch,
since the new one is exact or fitted differently. It is promised for
closed-form strips (cylinders, tori), which both paths build exactly. A case
the old path refused may now build. A case the old path built and the new
one refuses is a migration blocker.

### Amendment: the general patch before the M2 switch (#902)

M2 first plugged the fixed-resolution `make_filling_n_region` patch into
unequal-radius and mixed convex and concave corners. On the corpus it missed
the bands' fit target (2e-4) at every such corner, with either hole: skew
strip ends, or square ends with conic trims on the faces. The patch stood up
to 0.29 outside the box, or met its strips up to 9.5e-3 rad off tangent;
refining its net did not mend it.

The cause is twist incompatibility. At a hole corner two strips of
different radius or sign meet, both tangent to the same face. A G1 patch's
shape operator there must be symmetric, which asks
`k_A (t_A . c_A)(c_A . t_B) = k_B (t_B . c_B)(c_B . t_A)` of the two
cylinders' curvatures `k` and circumferential directions `c` along the hole's
sides `t`. Different radii or signs break it, and the normal part of the
mismatch cannot be mended by any choice of cross tangents. A least-squares
polynomial fill cannot take both supports' tangent planes there; a
Gregory-type patch exists to.

So the three-sided case of M4's patch moves ahead of the switch, with the
natural hole at every corner (§4.4). Mixed and unequal corners stay on the
old path until it meets 2e-4 at a fixed resolution, and the switch waits for
the whole corpus, as before. At corners of radii 1, 1.5 and 2 and of 1, 3
and 2 (#713's) and at an L-prism's mixed corner, the B-spline at a net of 33
meets its strips within 6.7e-5, 7.7e-5 and 5.5e-5 rad and 4.3e-6, 6.4e-6
and 2.2e-6 in position, and is proved inside its region. Straight ribbons
stood the 1, 3 and 2 corner 3.0e-4 past a face: off a strip tangent to a
face, a straight line crosses it, so each ribbon bends as its strip does. A Gregory patch converts exactly to rational Bézier form,
but with base points at its corners (weight zero there), which a
tensor-product B-spline kernel cannot evaluate. The interpolated form is
taken instead, its error measured on it.

### Amendment: chamfer corners are G0 (#902, M4)

§4.4 had the off-plane chamfer corner built as the general case with
straight sections, on the round side's transfinite patch. That patch meets
each side's support tangent along the side, which a round's hole allows:
at each hole corner the two strips are tangent to the face they both touch,
so their tangent planes agree there. A chamfer's hole does not. At each
corner two bevels meet on a face at an angle, and their planes share no
tangent plane at that point, so no surface is tangent to both. On a
rectangular pyramid's apex chamfered on its four slant edges the patch, its
supports the bevels, misses the bands' fit target in angle, as it must.

A chamfer is G0 everywhere, between its bevels and its faces alike, so its
corner patch is G0 too: it takes each bevel's end segment exactly and need
not meet the bevels tangent.

- **Four corners:** the bilinear patch through them. Its four sides are the
  segments between them, it lies in their hull, and it halves their
  tetrahedron, which gives an exact volume oracle (#952). SMLib builds the
  same corner [14].
- **Five or more:** the corners blended by the Wachspress coordinates of
  the regular polygon of as many corners, `S(x) = sum_k l_k(x) P_k`. On each
  side of the polygon the coordinates are its two corners' alone and linear
  between them, so the patch takes each straight side exactly; on a square
  they are bilinear, so four corners give the bilinear patch again; where
  the corners are coplanar the patch is their plane. The face carries it as
  the round side's five-or-more patch is carried, a clamped bicubic
  B-spline on each side's sector, its pole at the polygon's centre,
  neighbours sharing their line exactly; its miss from the exact blend is
  measured and held to the bands' fit target.
- **P5:** the patch is proved within the material every face at its vertex
  bounds.
- **A curved bevel there:** its end is where its surface meets the plane
  through its two corners and the hole's middle (a cone's conic), and the
  patch takes curved sides by adding each side's departure from its chord,
  weighted by its two corners' coordinates. On side `i` the sides beside it
  stand at an end of theirs, where a curve meets its chord, so every side is
  taken exactly; on a square this is the bilinearly blended Coons patch.

### Amendment: what M6 refuses by name (#902, M6)

When M6 deletes the old paths, the engine builds a request or refuses it
by name. The rule: refuse only what the mainstream applications refuse as
well (the survey below, of their public documentation); build the rest.

**Pieces of one edge asked with different chamfers.** A boolean or an
import can leave one edge of the form as several edges that continue each
other: collinear, or pieces of one curve, meeting at a vertex of no other
edge. Asked alike, the pieces chamfer as the one edge (the seam layer on
planes, #915). Asked with different chamfers:

| Application | Different sizes on edges that continue each other |
|---|---|
| SOLIDWORKS | A tangent chain takes one value ("tangent propagation"). Different sizes in one feature only for symmetric offset-face chamfers ("Multi Distance Chamfer"); the multiple-radius fillet note says sizes cannot differ on faces or loops with common edges. [help.solidworks.com, Chamfer PropertyManager; Constant Size Fillets] |
| Fusion | "Tangent Chain" is on by default; cleared, each piece takes its own selection set with its own size and the feature builds. The junction's shape is not documented. [Fusion Chamfer reference: help.autodesk.com/cloudhelp/ENU/Fusion-Model/files/SLD-REF-CHAMFER.htm] |
| Onshape | One set of distances per feature; per-edge overrides change only the direction. [Onshape Chamfer help] |
| Creo | Separate chamfer sets are allowed, with documented "transitions" where sets meet (Stop Case, Stop at Reference, Blend, Continue, Corner Plane, Intersect, Patch); Blend and Continue are shown "between two collinear chamfer pieces". Which one applies between two sets of different sizes is not stated. [Creo 9.0 help, About Chamfer Transitions: support.ptc.com/help/creo/creo_pma/r9.0/usascii/part_modeling/part_modeling/About_Chamfer_Transitions.html] |
| NX | Not documented publicly (the help sits behind a login). |
| Solid Edge | Refused within one feature: "multiple edges in one operation if they have the same setback value… you must construct separate chamfer features for each chamfer size". [Solid Edge ST6 help] |

Practice treats a chain of edges that continue each other as one edge with
one value by default. Onshape and Solid Edge refuse different sizes inside
one feature, but Fusion (separate selection sets) and Creo (separate sets,
with transitions between collinear pieces) build them, and every
application builds two chamfer features on the two pieces. No application
documents a taper between two sizes there; a size that varies along an
edge is what a variable chamfer is for. So (revised in M6: the step
replaced #1003's refusal):

- **Rule:** pieces of one edge asked with different chamfers are built, not
  refused: each piece is an edge of its own, chamfered as asked, and the
  bevels step at the vertex between two pieces, each ending on the plane
  through it square to the edge, as two features on the pieces would build.
  The step's face is the deeper bevel's end section less the shallower's;
  where two distances swapped between the pieces make the sections cross,
  it is a triangle on each face's side of the crossing, each facing the
  side whose bevel is deeper. A piece asked alone still stands for the
  whole edge (one value per tangent chain, the default above); once pieces
  are asked apart, an unasked piece stays sharp. On planes the rings build
  the step; pieces of a curved edge asked apart go to the one boolean
  (amendment "M5b"), each piece's tool ending on the plane square to the
  edge at the joint: rims on any coaxial surface of revolution (a cylinder,
  cone, sphere or torus) step on the plane through the axis, creases and
  straight edges on curved solids on their strips' closed tools. Only a
  pair no tool closes is refused by name.
- **Degenerate planar topology** (an edge bounded by other than two faces,
  as where two blocks touch along an edge; an edge closing on itself; two
  edges between one pair of vertices; a vertex whose faces close round it
  more than once, as where two blocks share only a corner in one shell)
  is refused by name. The mainstream
  applications do not model such a body, and a chamfer needs every edge
  between exactly two faces.
- **A compound holding anything but solids** (a shell, a face, a nested
  compound, or nothing) is refused by name: a chamfer cuts solids. A
  compound of solids is chamfered solid by solid (#960).
- **A rim's round as wide as its cap, or wider** (a drum's rim rounded at
  its own radius or more: the ball would touch the cap at or past the axis)
  is refused by name. SOLIDWORKS diagnoses a fillet whose radius is too
  large for the face it rounds onto [help.solidworks.com, Fillets
  diagnostics, "fillet radius too large"]; Onshape refuses a fillet larger
  than the face it lies on, reporting "radius of fillet on face too large"
  [forum.onshape.com]. No application documents a round that consumes its
  cap past the axis.
- **Rounds that no ball fits** (a ball set back past its face's far side
  along an edge that shares a face or a vertex with another of the
  request, two bands overlapping across a face they share, a plate rounded
  past half its thickness) are refused by name, as every mainstream
  application refuses a fillet whose faces cannot hold it. A request no
  engine route builds is checked before it is refused: a request of two
  edges or more as a whole (its corners read, each ball seated, no two
  bands overlapping, a concave band running past its face's side kept for
  the cut-off at that side), a request of one edge as its entry reads it
  (`request_fits` through `fillet_edges`, the ball on its two faces through
  `fillet_edge`). One failing is refused in that check's own words, or,
  where the ball rests on a narrow face's far edge, in the words the
  engine's reading of that run gives for why it does not hold; one passing
  is refused as a round the engine does not build. With M6 there is no
  other construction to hand it to.
- **Chamfers past what the part holds** (bevels whose setbacks on one face
  cross, taking it past one of its own sides: a plate 6 thick bevelled 4.9
  on every edge) are refused by name where no route builds them, the
  setbacks read as the engine reads them. A face only taken to a line or a
  point, or a bevel rolling over a narrow face that the notch builds, is not
  refused. SOLIDWORKS's chamfer diagnostics hold that chamfer faces cannot
  eliminate adjacent faces [help.solidworks.com, FeatureWorks diagnostics,
  "Chamfer"]; Fusion reports "Chamfer could not be created at the requested
  size" where a distance runs past the faces it bevels
  [forums.autodesk.com]. Before, the old facet path refused the request
  for an unrelated reason (a facet closing a concave corner).

### Amendment: M5a, the one-boolean fallback (#902, M5a)

§4.5 keeps one boolean for what the one assembly cannot place. M5a builds
it for planar chamfers first.

- **The blend solid, as convex pieces.** Each bevel takes the material
  past its plane between its edge and its two contact lines, as far as its
  setback corners: the convex hull of the edge's ends and the four corners
  (a step's corner on its own line where a face steps). Each vertex of two
  bevels or more takes the corner between the vertex and its facet: the
  hull of the vertex and its faces' corners there. These are the regions
  the assembly's rings bound, so where both build a request they agree.
- **Applied by booleans.** Concave pieces are fused into the solid, then
  convex ones cut from it, one boolean per piece, as the current path
  applies its wedges: pieces that share faces (two bevels' ends at a
  vertex) are not fused into one tool first, since the boolean refuses
  same-domain contact between its arguments but melts a tool's planes
  into the solid's own. At a vertex mixing convex and concave bevels the
  corner piece is cut where the facet faces the vertex (the vertex stands
  outside the result) and fused where it faces away.
- **Held to the rings' domain.** A bevel whose contact runs past its
  face's far side, a setback corner standing off its face, or a bevel the
  others leave nothing of is not built: these are past what the part
  holds, refused by name as before.
- **History** is traced through the booleans: each tool names the faces
  that carry a blend's face (a bevel's, a facet's) with the edge or vertex
  they come from, and what the booleans make of those faces is generated
  by it.
- **Shared with rounds.** `fallback::apply` takes closed tools (a solid, cut
  or fused, and its named faces) from any section law: the round side
  closes a round's strip, its two flaps and its end caps into its tool
  (M5c), and the ordering and history are this one.
- **Rounds (M5c).** A planar request of rounds alone that the non-local
  test labels goes there too. Each round's tool is the region between its
  edge and its strip: the strip is the engine's own round of the same edge
  on the wedge its two faces and its end faces bound alone (the hull of the
  edge and its setback corners, pushed twice as far out, which no other
  face reaches), closed by its two flaps (its faces' planes between the
  edge and its contact lines) and its two caps (its end faces' planes). At
  a concave edge that wedge is the room in front of the edge, whose edge is
  convex, and its round's spandrel is the concave round's fill, fused.
  Every contact must lie on its own face end to end for these tools; one
  running where its face does not reach (a wall cut away under a boss)
  goes to the general route below (amendment "M5c, non-local rounds keep
  their surfaces"). Rounds joined at vertices, or a round ending where other than three
  faces meet, close into one tool where the solid is convex round every
  vertex they end at: the convex polytope the faces round those vertices
  bound (within a box well past the setbacks), less that polytope rounded
  by the engine, one boolean between the two (which share every face but
  the rounds'). Their strips and corners are the engine's own, as where
  nothing else reaches them. Such rounds at a vertex not convex all round
  are not taken yet: the sheet builds them, as before.
- **Curved strips closed (M5c).** `trimmed::close_tool` closes any strip
  the trimmed assembly reads (`Strip`: contacts, surface, section) into
  its tool: the strip, each face's piece between the edge and its contact
  on the face's own surface, and a cap at each end on a plane, either the
  strip's own end (the one other face there, a plane square to the edge)
  or a joint where another strip runs on (`EndCap::Plane`, square to the
  edge); each cap's sides straight on the faces. The chamfer side's
  curved classes close their strips with it too. A rim's round reaching a
  face it does not end at (a face standing under the round within its
  contact ring, or one not proved outside the box it can take up) is
  closed so, arc by arc, read on the solid with its seams dissolved and
  applied on the solid as given, each tool's strip named for the edges
  asked. A straight round along a cylinder's ruling is labelled so too:
  every face but its two and those at its ends proved outside the box it
  can take up (between its ends, within the spine's cylinder through the
  edge, and at a convex edge on the edge's side of the spine and inside
  both faces' tangent planes there); its caps' sides on a wall square to
  the cap are the arcs about the wall's axis. A strip along one closed
  edge (no ends to cap) declines for now.
- **The general route (M5c).** Where no strip's tool closes a non-local
  request of rounds, its far feature comes off a copy of the solid
  (`ogeom_bool::remove_faces`): the faces the rounds reach without ending at
  them, grown across their outer rings, never into a face the rounds end at
  or a face the feature stands in through a hole of that face. The engine
  rounds the copy as the local request it is there; the copy less its
  rounding is cut from the solid as given and the rounding less the copy
  fused, each strip's faces named for its edge. The rounds keep their own
  surfaces and the feature they reach is changed, as with the strips' own
  tools. This takes a corner not convex all round, strips meeting at an
  angle on a curved solid, a strip along one closed edge, and one along a
  cone, and a corner of mixed convexity (concave and convex rounds at one
  vertex), whose signed corner's patch the boolean of the copy less its
  rounding closes along the patch's pole.
- **What takes it.** The non-local test (P2/P5) routes there, and so does a
  local request the assembly cannot place: chamfers set back differently
  in one line where a third face meets the step (the step's face is then
  that face, less and plus the two bevels' end sections, which the
  booleans resolve); a consumed face's hole not bevelled all round; a
  bevelled concave edge at a consumed face's hole or reflex corner; a
  straight bevel between planes on a solid with curved faces reaching a
  face it does not end at (each bevel's hull, the reach proved by
  `certify_face_within`); a bevel of a solid with voids reaching a void
  (the outer shell chamfered alone less each void's body, a boolean a
  void, exact since what the bevel takes there is already empty); a
  mitred request (`ChamferCorner::Point`) reaching a face it does not end
  at, each bevel its edge's whole prism, whose union at a convex vertex is
  the mitre to the bevels' shared point. This is a boolean on a local request, which §4.5
  rules out for the one assembly; it is kept to the cases the rings cannot
  place, each listed here as it is routed, so that M6 deletes the current
  path with every planar chamfer built through the engine.

### Amendment: M5b, a chamfer rolling over its face (#902, M5b)

A chamfer whose contact runs past its face's far side (a leg longer than
the face is wide) was refused by name ("runs past its face"). The rule is
to refuse only what the mainstream applications refuse (amendment "what
M6 refuses by name"); their public documentation:

| Source | A chamfer meeting an edge it was not seated on |
|---|---|
| T-FLEX CAD 17 (Parasolid) [P27] | Overflow processing applies to chamfers as to roundings, but for the cliff-edge options: "Notch option is applicable both to Roundings and Chamfers", "Process Smooth Overflow option is applicable both to Roundings and Chamfers", "Retain Cliff Edges option is not applicable to chamfers". Notch keeps the blend and trims or extends the element it overflows. |
| SOLIDWORKS [P32] | Its "cliff chamfer" ("the chamfer face replaces one of the support faces") is a FeatureWorks recognition diagnostic for imported bodies, not a modelling refusal. |
| Solid Edge, NX, Onshape, Fusion, Creo | Not found. |

So a chamfer rolls over as a notch: its bevel carries on past the face's
far side and trims the faces it meets there, the face it ran past gone
where the plane passes all of it. On planes the rings build it: a bevel
rolling over (its setback past its face's width) is marked so in P0, and
its plane may cut any ring of its cluster (`Cluster::rolling`), the faces
beyond included; a request rolling over is not labelled non-local (the
notch trims what it reaches, by definition). On a convex plate every edge
chamfered past the plate's thickness is the polytope of the plate's, the
bevels' and the facets' planes. A rim between a planar cap and a coaxial
cylinder, cone or sphere the curved path's rings do not place (rolling
over its wall, reaching a face it does not end at, or pieces of one rim
asked apart, each stepping at the joint on the plane through the axis)
goes to the one boolean: its tool is the bevel's meridian section turned
about the axis over the rim's own arc, or once round the whole turn where
the arcs of one circle asked alike close it (so the bevel is one face and
no two tools meet on a plane through the axis). A rim between any two
surfaces of revolution about its axis that no other route builds (a
torus's rim, its arcs asked apart) takes the same tool, its section the
faces' traces in the half-plane (a line on a plane, cylinder or cone, an
arc on a sphere or torus) and the bevel's line, turned by
`make_revolution`. Curved bevels that meet each other and roll over or
reach a face are built on a convex solid (every face a plane or a
cylinder, cone or sphere holding the material inside, every edge
convex): there the notch takes everything past each bevel's surface, so
the solid is taken in common with each bevel's material side, a
half-space or the inside of a cone or cylinder about the rim's axis, as
long as no vertex has three asked edges (those close with a facet).
Elsewhere such bevels take their local tools, below. Any other curved bevel the rings do not
place (a crease rolling over or reaching a face, or pieces of one crease
asked apart) goes to the one boolean too where it stands alone: its strip
closed into its tool (`close_tool`, shared with rounds, a chamfer's caps
cut along the segment between its contacts), capped on its own end faces
or, at a joint between pieces of one curve, on the plane square to it. A
bevel the others leave nothing of
stays refused by name where the one-boolean fallback builds (amendment
"M5a"), as before.

**What a rolling chamfer takes.** Fusion defines a chamfer as material
removed from an exterior edge, or added to an interior one [P34]; T-FLEX's
notch keeps the blend and trims the element it overflows [P27]; SOLIDWORKS
documents its overflow types for fillets (the blend trimmed against the
adjacent surface, or the adjacent edge kept) [P30], and none for chamfers.
So a chamfer rolling over takes the material between its edge and its
bevel, and nothing else: its **local tool** is bounded by the bevel, the
two seat faces' surfaces carried on past the faces, and at each end of the
edge either the neighbouring bevel's tool, where two bevels meet at the
vertex, or the plane square to the edge there. The tools at a vertex are
taken together and the union applied in the one boolean (cut at a convex
edge, fused at a concave one). A wall rising at an edge's end is cut only
where the tool reaches it (the wall square to the edge there is the cap's
own plane, so not at all), and a far wall only where a tool reaches it.

- **Planes.** #1028's rings build the notch where they place it: a rolling
  bevel's plane cuts the rings of its cluster, which is the same material
  as its tool where the rings build (a plate's edge rolling over, a plate's
  edge ending at a rising wall, a reflex corner mitred without a notch).
  Where the rings decline, the union of the tools is applied. A bevel
  rolling down a front whose foot meets a floor at a concave edge is
  refused instead (amendment "a bevel stops at a concave floor").
- **Curved solids.** Where the rim and crease tools decline because the
  asked edges meet, and the solid is not convex, the union of the tools:
  a rim's tool is its section turned over its own arc, any other strip's
  is closed by `close_tool` with its caps square to the edge at its ends.
- **Corners.** Mainstream applications close every such corner: Creo's
  chamfer transitions name the corner closings ("Corner Plane", "Patch")
  where chamfer sets meet, and the planar engine closes a vertex of three
  bevels or more with the facet through their setback corners by default.
  So the rolling tools close them too:
  - At a convex corner of two bevels the tools capped square at the vertex
    cover it between them.
  - At a reflex corner each tool runs on past the vertex until it meets the
    other's. The piece they meet in is beyond both caps, on the edge's side
    of both bevels and on each seat face's side that holds the strip's other
    contact. Between plane bevels on planar faces it is the polytope of
    those planes; otherwise the box round the vertex is taken in common
    with each side kept (a half-space, or the inside of a cylinder, cone or
    sphere) and less each side dropped, the box grown until the piece
    stands clear of it. A kept cone stops short of its apex, and the piece
    must too. Concave edges' pieces are fused like their tools. The piece
    is bounded by the bevels carried on and the face both edges share, not
    by each edge's other face: its lump reaching the vertex is kept, so each
    bevel runs on until it meets the other, flush, as the planar case
    closes (on a tab whose arc meets a straight edge at a reflex foot, no
    sliver of the square stays outside the tab's cylinder).
  - At a vertex of three rolling bevels or more the facet closes the corner:
    the plane through the setback corners (where the curves bounding each
    face at the vertex meet, an asked edge's contact or the edge itself),
    or, where four or more corners leave one plane, a fan of triangles about
    their centroid (G0, amendment "chamfer corners are G0"). Its piece is
    what lies inside every bevel on the vertex's side of the facet: between
    plane bevels the tetrahedron of the three setback corners and the point
    the bevels share. The planar rings build the same corner where they
    place it.
  - **Refused by name:** two rolling bevels meeting at a reflex corner, one
    convex and one concave, and a corner of three or more mixing convex and
    concave bevels (a cut and a fill meet there, which no one tool closes);
    a corner piece that does not stand clear of a box of 64 times the
    strips' reach. A bevel the others leave nothing of stays refused, as in
    amendment "M5a".

- [P32] SOLIDWORKS 2018 Help, Chamfer - Cliff Chamfer:
  https://help.solidworks.com/2018/english/SolidWorks/fworks/IDS_FW_DIAGMSG_CHAMFER_CLIFF.htm
- [P34] Autodesk Fusion, Chamfer reference:
  https://help.autodesk.com/cloudhelp/ENU/Fusion-Model/files/SLD-REF-CHAMFER.htm

### Amendment: M5c, non-local rounds keep their surfaces (#902)

Three behaviours of the non-local rounds were chosen in M5c and are
settled here by survey of the mainstream systems' public documentation.

**The question.** A round whose ball, rolled on its two faces, meets
another feature of the part: its fill reaching over a pocket's mouth, its
spandrel holding a hole, or its contact running where its face does not
reach (a wall cut away under a boss, the ball touching the wall's plane
inside the boss).

**What the systems do.**

| System | Behaviour |
|---|---|
| Creo [P33] | "When you create a convex round and attach it as a solid or as a same quilt, the radius of the round is maintained and the adjacent geometry is trimmed." For a concave round, "If you attach the round as a solid, the neighboring surfaces are extended." (Extended surface rounds: a radius larger than the adjacent surface.) |
| SOLIDWORKS [P30] | Overflow type **Keep surface**: the fillet stays smooth and the adjacent edge is changed; **Default** picks it or keep edge from the convexity met. |
| Parasolid family [P27, P28, P29] | After smooth overflow and retain cliff edge, **notch**: the round is kept and the element it overflows is trimmed or extended to meet it. |
| Onshape [P7] | **Allow edge overflow**, on by default: the fillet may change the edges it meets to stay one smooth surface. |
| Fusion | Its fillet documentation describes no overflow option; not found. |

**Rule.** A non-local round keeps its own surface, the one its ball rolled
on its two faces sweeps, and the feature it reaches is trimmed or
extended to meet it. Hence:

- a convex round's spandrel is taken whole, a hole within it going with
  it, and the part beyond its contact trimmed by its surface (the boss
  over a cut-away wall keeps the spandrel cut under it);
- a concave round's fill is added whole, roofing a pocket's mouth where it
  reaches over it (the neighbouring surfaces extended), the pocket below
  kept.
- a contact crossing a hole in its face (a cross drill coming up into a
  rim's band) runs on across it: the hole's faces are the far feature,
  taken off the copy and put back, so the hole cuts the round's surface.
  A contact leaving its face across the face's outer boundary is the
  roll-over of §13, question 8, and is handed back as before.

Retaining the cliff edge instead (the ball resting on the boss's
underside edge) needs the ball to pass the vertex where that edge meets
the wall's cut, a vertex seat this milestone refuses (§13, question 8;
amendment "overflow handling" holds an edge only where no vertex seat is
needed);
the notch, which every surveyed system with an overflow option builds, is
built instead. The engine builds all three by the general route of
amendment "M5a" (the far feature off a copy, the rounds built there, the
differences applied by booleans).

- [P33] Creo Parametric 10.0 Help, About Extended Surface Rounds:
  https://support.ptc.com/help/creo/creo_pma/r10.0/usascii/part_modeling/part_modeling/About_Extended_Surface_Rounds.html

### Amendment: marched seats through the engine (#902)

A round whose crease no closed-form section reads (a spline edge, the
seam of two cylinders crossing, a converted solid's every edge, the
rims of a pinched pair of drums) has its strip from the march: the
ball's arcs marched along the crease and fitted into the band, closed
into its wedge by its legs on the two faces and its caps. That
construction exists (`marched.rs` and the constructions it calls) and is
not rewritten; the engine takes it as the strip builder for these seats.

- **Scope.** A request one edge of which at least is a marched seat,
  where every ball fits beside its edge (no strip taken off whole, no ball
  on a far edge). Each edge's wedge is built on the solid as given, with
  the request's other edges as its mates, and set aside unapplied.
- **Junctions.** Where two of the edges' creases turn at a shared vertex
  (a rib's spline crease into a boss's rim), the bands meet as the
  current path's junctions meet them: each fitted to the plane the two
  meet in (carried on past its cap across a re-entrant turn, trimmed back
  across a convex one), a flush pair joined across it into one wedge, a
  crossed re-entrant pair's corner built from the two. A marched edge
  meeting one of another kind anywhere but at such a junction, or with a
  third at a vertex, is still built after it on the current path.
- **Application.** The wedges are the one-boolean fallback's tools
  (amendment "M5a"): cut where the crease is convex, fused where it is
  concave, one boolean each, every face of a wedge named for its edge
  (or the edge its builder credits it to). A result that does not check
  valid is handed back.
- **A strip taken off whole.** One marched edge whose ball does not fit
  beside it, or whose band the boolean refuses, takes the strip beside it
  (a bevel) off whole as the current path reads that round: the corner of
  the faces either side of the strip, marched, its wedge set aside and
  applied as the engine's tool; a ring beside a rim likewise. Where the
  ball does not fit and no such round is read, the current path's refusal
  is the engine's, by name.
- **A ring beside a rim.** One rim (a revolved seat) beside a ring a
  round takes off whole (a chamfer or a smaller round on the rim, as
  `consumed_rim` reads one): the round of the cap and wall either side of
  the ring, its wedge set aside and applied by the seated booleans, which
  trim the faces a rim's band sits on.
- **Still the current path's:** a band the boolean refuses where no strip
  is taken off whole (its refusal follows the current path's retries).
- **A leg on a fitted band.** A face beside the crease may be an earlier
  marched blend, such as a bevel round a bore's rim. Its surface closes
  round `v` alone and is not periodic. The leg on it is built the way every
  leg is built: its chart is transposed so its loops wind in `u`, and the
  face is then rebuilt on the patch's own chart. The loop's first station
  stands on that patch's seam, so the leg's chart stays within the
  patch's domain, where the boolean finds the solid's face it lies on.
  Where the ball fits beside such a bevel, the round rolls between the
  bevel and the face beyond it, and the bevel stays narrowed. This is what
  mainstream applications build. The bevel comes off whole only where the
  ball does not fit beside it (`consumed_marched.rs`).

### Amendment: the marched seats' tools are the engine's own (#902)

The amendment above took the marched seats' wedges from the current
path's edge-by-edge dispatch, which built each wedge and applied it unless
a caller was collecting. M6 deletes that dispatch, so the engine builds
these tools itself. What is accepted, refused and built does not change;
only who builds it does.

- **Seats.** The engine reads each edge's seat: a straight crease between
  planes or a ruled pair, a rim on a planar cap and its coaxial wall, or
  the march's. The fit and overlap guards (`request_fits`) are the
  engine's, and blame and the largest-fit searches ask them there.
- **Strip builders.** The march's band (`marched.rs`, `pinched.rs`) and a
  rim arc's exact torus band (`rim_arc.rs`) hand their wedge to a sink the
  engine passes them; nothing is intercepted. Junctions are fitted and
  joined by the engine (`engine/marched_junctions.rs`) with the junction
  geometry (`junction.rs`).
- **A strip or ring taken off whole** is built as its rule says: the plain
  crease's round. The strip comes off a copy of the solid
  (`remove_faces`, its two neighbours grown back to meet along the crease
  the strip cut), the engine rounds that crease at the asked radius, and
  the copy rounded is the result where every point of the strip lies on or
  past it on the round's side (the strip wholly in what the round takes).
  Elsewhere the round is not this one, and the one-edge refusal stands by
  name. One construction serves a bevel beside a marched seat and a ring
  beside a rim, and needs no seated booleans. A bevel round a bore's rim on
  a drum comes off so: the drum and the bore grow back to meet along a
  closed spline, cut where their seams reach it, each regrown face holding
  both columns of its seam, and only the edges the removal made are
  rounded (the bore's far rim is the two faces' too).
- **Checked.** The old edge constructors (those that built and applied
  their own wedge, or read a strip off whole for the one-edge entries) are
  removed, so the marched seats have no old code to go back to.

- **A cone rim fitted.** Where no route closes a rim on a planar cap and
  a coaxial cone asked alone, its exact torus band is taken by the march
  instead, fitted through the ball's arcs, and the request goes through
  these seats as above. This is the current path's retry
  (`with_fitted_rims`): the boolean resolves some meetings of the fitted
  band, within its fit of the torus, that it does not of the exact one
  (a drum drafted out, on a block whose leaning crease an earlier round
  has rounded).

### Amendment: mixed rolling corners (#902)

Chamfers rolling over their faces (amendment "M5b") that meet at one
vertex, some convex and some concave, were refused by name: their corner's
facet is not one tool's.

**What the systems do.**

| System | Behaviour |
|---|---|
| Creo Elements/Direct Modeling [P40] | "You cannot create two chamfers that meet at a vertex, where one is concave and one is convex." |
| COMSOL (Parasolid kernel) [P41] | A "Y-shaped chamfer at a vertex where three or more edges meet and there are two chamfer surfaces of different convexity"; "in some cases, using this option is necessary for the operation to succeed". |
| ACIS [P42] | At a vertex of mixed convexity a convex vertex chamfer is assumed, and a vertex of three edges of mixed convexity "may not chamfer successfully". |
| Cimatron [P43] | Vertex selection is ignored where the result is added material on a concave and convex chamfer at one vertex. |

**Rule.**

- **Three or more** at a vertex, of mixed convexity, are built (the
  Parasolid-based system builds them). The corner is the rings' own: the
  facet through the setback corners (ADR 0001, "chamfer corners are G0").
  Through the union of local tools, the vertex's side of the facet takes
  the sense most of the vertex's bevels take (chamfering the complement
  gives the complement, so the rule is the one that commutes with it).
  The piece is that side inside every bevel, less the solid where it is
  fused (within it where it is cut), its lump at the vertex. Each bevel of
  the other sense spares that side; each of the piece's sense spares the
  far side. On wide faces this is the rings' corner exactly.
- **As many convex as concave** at one vertex is refused by name: the
  bevels do not settle which side of the facet the material keeps.
- **Two** meeting at a reflex corner, one convex and one concave, stay
  refused by name, as Creo Elements/Direct refuses them and no surveyed
  system documents building them.

- [P40] Creo Elements/Direct Modeling, Modify3D Chamfer:
  https://support.ptc.com/help/creo/ced_modeling/r20.2.0.0/en/ced_modeling/OSDM_Main/Modify3D_Chamfer.html
- [P41] COMSOL Design Module User's Guide, Chamfer:
  https://doc.comsol.com/6.3/doc/com.comsol.help.design/design_ug_function.5.40.html
- [P42] ACIS, solid:chamfer-vertex:
  http://www-isl.ece.arizona.edu/ACIS-docs/HTM/DATA/BLND/BLND/02SC/0033.HTM
- [P43] Cimatron, Chamfer options and results:
  https://help.cimatron.com/en/2026/Chamfer_Options_Results.htm

### Amendment: narrow faces and extended flush ends

The narrow-face cut-off requests use the existing **keep-surface/notch**
semantics, not a new refusal based on the width of a neighbouring face.
The public evidence is:

| Source | A round reaching the end or side of a neighbouring face |
|---|---|
| Creo [P33] | Extended surface rounds maintain a convex round's radius and trim adjacent geometry. Concave solid rounds extend neighbouring surfaces. Stop, Blend and Intersect transitions are supported. |
| Onshape [P7] | Allow edge overflow is on by default and permits changing the neighbouring edges to keep a smooth, continuous surface. |
| SOLIDWORKS [P30] | The existing survey records Keep surface as maintaining the round and changing the adjacent edge. The current public page serves a JavaScript shell rather than its example text. |
| Fusion [P10] | Documents material removed at exterior edges and added at interior edges. No cut-off or overflow end rule is stated. |
| NX [P29] | The existing survey records maintaining the blend and moving sharp edges. The public API URL could not be retrieved for this amendment. |
| T-FLEX [P27] | Notch leaves the rounding surface undeformed and trims or extends the affected element's surfaces to it. |

The papers supply the geometric construction, not an application dialog's
end default: the offset intersection gives the constant-radius strip
[1, 2], and trimming resolves its interaction with other faces [12].
In particular, an old face's outline is not a new limit on the round's
radius once that face is extended.

- **Flush ends.** A concave strip extends its neighbouring face to its end
  section, in that face's plane even outside the original outline. Its
  contacts must still end on its two support faces. The same rule applies
  where two cylindrical strips meet in a crease: their crossing must meet
  the end section in the far face's plane, but that face may extend to the
  crossing. A convex strip retains amendment "M5f"'s local-material proof
  where its end leaves the original face.
- **Support first.** Amendment "M5d"'s far-edge round takes precedence over
  the ordinary sheet and its cut-off tools: where the ball rests on a
  narrow face's far edge, the changing support determines the round, not
  a full cylinder extended across that face.
- **Independent overflowing strips.** Where a contact leaves a support
  face over part of the edge and the sheet cannot place its rings, the
  engine closes its own cylindrical strip into a finite tool. Its ends
  are the sections square to the guiding edge at its original vertices,
  so it does not round material past the selected edge. The convex tool
  is cut from the original solid; the concave tool is fused, extending
  the neighbouring faces to meet it. This is the finite-tool convention
  of the non-local route, not a projection of the contact back onto a
  narrower face and not a reduction of the requested radius.
- **Domain and guards.** The original ball-fit and band-clearance
  predicates apply before tool construction. A face taken off whole is
  not this class. Bands sharing a vertex still use the sheet's common
  vertex blend. The migration retains the existing named overlap and
  mixed-sharp-corner refusals without enlarging their domain. The tools
  and result are closed, and the result passes `check`. An exposed flap
  of a concave tool can meet its cylinder at a zero-angle edge: the
  sheet's opposite-normal fold test is not a test for this solid's
  validity.

Independent checks are the roof's closed-form spandrel moment, the
straight dihedral's exact section area times its length, and a voxel
volume bound at two concave strips' joint. The latter uses only the
original three planes and the two tangent circles. Plane distance and
axis distance are 1-Lipschitz, so each voxel's half diagonal bounds the
uncertainty, independently of the blend and volume implementations.

### Amendment: round section laws, conic and elliptical (#902)

The round side's section laws beyond the rolling ball: a round meeting
each face its own distance, its section an ellipse or a conic of a given
`rho` (`fillet_edge_distances`, `Blend::RoundDistances`).

- **The law.** `Section::Conic` names each face with its distance, as a
  chamfer does; P1's setbacks are those distances. The strip is the
  section swept along the edge (the planar section's own conic spec, the
  current path's exact section), built by the rounds' sheet.
- **Its vertices.** Alone at a vertex of three planes it ends flush in the
  third face, as a ball's round does. With rounds at every edge of a vertex
  it closes on the transfinite patch: its section is no ball's, so no
  sphere, ridge or mitre is read from it. Where its corners there stand at
  different stations its end is skew, the strip (an extrusion along its
  edge) cut by the plane through both, exact as the ball's. Beside an edge left sharp, or
  reaching a face it does not end at, it is not built yet and goes on as
  before.
- **Routing.** `fillet_edge_distances` and `blend_edges` (rounds of both
  kinds, no bevels) take the engine first.
- **On a rim** (a circle between a planar cap and a coaxial cylinder or
  cone, `fillet_edge_distances`): the rings stand at the two distances, the
  cap's along it and the wall's along its generator, and the strip is the
  meridian conic (its control point the rim's own point, its weight the
  profile's for the legs' angle) turned the whole way about the axis, a
  surface of revolution, built by the rims' sheet as the ball's torus is.
- The patch reads the normal and curvature of any parametric strip surface
  from its chart, an extruded conic among them.
- **A varying radius** (`fillet_edges_variable` on one straight edge
  between planes): `Section::Varying` carries the law and the vertex its
  fraction starts at. Each face's distance is the ball's radius over the
  half-angle's tangent along the law, the section the circle (the
  elliptical section of equal distances), the strip the current path's
  exact varying band. Each end sets back its own distance (P2's corners
  reset to the law's there, on end faces square to the edge), and the strip
  ends flush; its chart's sense is read where it stands. Its contacts are
  the section's contact curves, bending with a law of points, each shared
  by the strip and its face.
- **A varying radius over a chain** of straight edges between planes,
  open or closed, the law laid over the chain by arc length from its
  first edge's start (or, asked from the middle of an open chain, from
  where the chain begins): each edge takes its own stretch of the law,
  held level past its ends, and the bands meet at the chain's turns as
  constant rounds meet there (amendment "M5e"): carried on through an
  inside corner, running out past a convex turn. At any other vertex a
  varying round is not built yet.

### Amendment: varying rounds on curved seats (#902)

A varying radius (`fillet_edges_variable`) over a chain the planar
engine does not read: a drum's rim, an open arc, a chain of arcs of one
circle, a drum's ruling, the seam of two drums, a marched crease, the
top loop of a pad whose sides are planes and whose corners are drums.

- **Survey.** Mainstream applications take a variable radius over any
  chain of edges, curved and closed ones included, and lay its points by
  arc length. NX stores a variable-radius point as a share of the whole
  chain's length since NX 11.0.1
  (https://support.sw.siemens.com/en-US/okba/PL8015439/The-variable-radius-points-cannot-be-specified-by-percent-Arc-Length-of-each-edge/index.html),
  and T-FLEX places points as percentages of the chain's total length
  (https://tflex.com/help/eng/T-FLEX%20CAD/17/14_edgeblend_types_roundvarrad.htm).
  SOLIDWORKS places control points by percentage of an edge's length
  (https://help.solidworks.com/2024/English/SolidWorks/Sldworks/r_Variable_Size_Fillet.htm).
  Inventor notes that a closed loop has no start or end point
  (https://help.autodesk.com/cloudhelp/2014/ENU/Inventor/files/GUID-A72E7D7B-80EB-4D4F-B6F9-C623104FE3E9.htm).
  Onshape propagates a fillet along tangent edges by default
  (https://cad.onshape.com/help/Content/PartStudio/fillet.htm).
  This is the law the current path already lays, and it does not change:
  - by arc length from the chain's first edge's start, held level past an
    open chain's ends;
  - periodic round a closed one;
  - each edge taking its own window of it;
  - the ball one radius from either side of a junction.
- **The band.** Each band is the march under the law (`marched.rs`), or
  the exact envelope on a straight edge between planes. An open chain of
  arcs of one curve is one band over them all, named for the first arc;
  the others go with it. Each band is closed into its wedge and built on
  the solid as given, with the chain's other edges as its mates, as
  amendment "marched seats" builds a constant one. Bands joined across
  tangent junctions are one wedge. A loop closing round a sharp corner
  opens on a straight edge between two of its junctions, the two halves
  meeting flush. Bands meeting at a convex or re-entrant turn are fitted
  to the plane they meet in, as the constant bands of the meeting's
  radius are.
- **Application.** The wedges are the one-boolean fallback's tools, as a
  marched seat's are: cut where the crease is convex, fused where it is
  concave, every face named for its edge. A result that does not check
  valid is handed back.
- **Still the current path's:**
  - a chain whose turn the bands meet in rounds, one built against what
    the other left (a drum's ruling turning into its rim);
  - a law's ball running off a narrow face (the cut-off family);
  - a law refused by name.

### Amendment: round distances on curved seats (#902)

A round set back its own distance on each face, with an elliptical or
conic section (`fillet_edge_distances`), over a seat the planar engine
and the rims' revolved conic do not read: a bore's rim on a drum, a rim
arc that stops, an open spline crease, a slanted edge whose dihedral
changes along it, a round running up under a roof or out along its wall
onto a ramp.

- **Survey.** Mainstream applications build asymmetric and conic rounds
  on any edge, curved ones included. The SOLIDWORKS tutorial's
  asymmetric fillet is on a cylindrical edge, with Distance 1 and
  Distance 2 and an optional Conic Rho profile
  (https://help.solidworks.com/2025/English/swconntutorial/t_tut_fillets_create_asymmetric.htm).
  NX's edge blend takes a conic section with two offsets and a rho
  (https://docs.plm.automation.siemens.com/data_services/resources/nx/12/nx_api/custom/en_US/nxopen_python_ref/NXOpen.GeometricUtilities.ConicCrossSection.html).
  Onshape's conic fillet is not symmetrical either
  (https://cad.onshape.com/help/Content/PartStudio/fillet.htm).
  The section is the current path's, unchanged: in the plane square to
  the crease, each contact its distance along its own face, the conic
  between them.
- **The band.** Each band is built as the current path builds it (the
  seated conic, marched where no closed form reads it, its ends as the
  constant round's), on the solid as given, closed into its wedge and
  collected rather than applied.
- **Application.** The wedges are the one-boolean fallback's tools, as
  in amendment "varying rounds on curved seats"; a result that does not
  check valid is handed back.
- **Tangent chains.** A `blend_edges` request of asymmetric rounds alone,
  meeting nothing but each other along tangent chains (the top loop of a
  pad, a pocket's floor), closes each chain as the current path does:
  its bands capped in their shared end sections, all set aside on the
  solid as given and applied together.
- **Still the current path's:**
  - asymmetric rounds in `blend_edges` beside rounds or bevels, or at
    their corners;
  - the largest-pair search;
  - refusals by name.

### Amendment: engine-owned variable-law and conic tools

The variable-radius amendments above and M5f's end-plane rule retain
their public meaning. Equal endpoint radii are a level law, not a reason
to enter a different constructor. The arc-length parameter, periodic
closed laws, endpoint holding and each edge's law window are unchanged.

- **Construction.** The engine closes the planar section's exact sweep
  into its own tool. It reuses the radius profile and varying-distance
  rational net as section mathematics, not the old solid constructors.
  Curved strips use the low-level march or exact revolved section;
  asymmetric-only tangent chains join those tools at shared end sections.
  This supersedes the collected-constructor implementation described in
  the two preceding amendments, including its narrow-face and curved-turn
  exceptions. Mixed asymmetric requests remain outside this migration.
- **End resolution.** Intersecting curved strips retain their original
  support surfaces. The surviving crease interval restricts the original
  law by its original arc-length fraction; the law is not restarted on the
  shorter crease. Guides resolve in a fixed order, straight before curved
  and then by source identity, independent of request enumeration. Planar
  end faces intersect the section's control hull to bound its extension;
  an overhanging end plane clips the tool, as M5f requires.
- **Compatibility.** No new public law or refusal rule is introduced.
  The earlier mainstream survey already supports variable circular and
  conic sections and edge overflow; Onshape's public fillet documentation
  also explicitly permits variable values for every cross-section type
  (https://cad.onshape.com/help/Content/PartStudio/fillet.htm). Existing
  named computation limits remain limits, not claims that mainstream CAD
  rejects the same geometry.
- **Evidence.** The variable and asymmetric suites retain their geometric
  and analytic oracles. An instrumented subprocess additionally exercises
  all three law APIs and curved conic entry points, requiring no old
  constructor event; result validity alone cannot distinguish the builders.

### Amendment: engine-owned mixed asymmetric requests

A request mixing asymmetric rounds with rounds and bevels is built by the
engine, not by a sequence of asymmetric API constructors. It follows the
section, tangent-chain and corner rules above, plus what a short survey of
user documentation not cited elsewhere in this record shows:

- **All edges blended means the vertex is blended.** T-FLEX CAD rounds the
  vertex when every edge meeting there is rounded in one operation, and
  recommends one operation over several
  (https://www.tflex.com/help/eng/T-FLEX%20CAD/16/blending_edges.htm). A
  convex corner of three planes whose three edges are all in the request
  closes on one setback patch whatever mix of circular, conic and bevel
  sections meets there, and the request's order carries no meaning.
- **Each edge has its own setback; the default is the least that fits.**
  Inventor sets a setback per edge at a vertex and computes the minimum
  allowable one on request
  (https://help.autodesk.com/cloudhelp/2014/ENU/Inventor/files/GUID-55B9D5C0-3859-4146-8866-B4E4AE89469E.htm,
  https://help.autodesk.com/cloudhelp/2014/ENU/Inventor/files/GUID-0C3675D8-7CE5-4664-AAE9-360EB324F363.htm).
  Each band's setback is the least the corner's widths admit, read from
  each band's width on each face, so an asymmetric band sets back by its
  own two distances.
- **Two of three blended edges of one convexity meet by intersection.**
  T-FLEX CAD leaves such a vertex unblended, with the two blend surfaces
  meeting along a new edge
  (https://www.tflex.com/help/eng/T-FLEX%20CAD/17/14_edgeblend_theory_2of3.htm).
  The engine does not build that intersection for an asymmetric band yet:
  an asymmetric round meeting another blended edge at a vertex where not
  every edge is blended stays refused by name, as before.
- **A bevel meets the corner in position only.** T-FLEX CAD cuts a chamfer's
  vertex with a plane offset from it; no source asks for tangency there.
  The patch is tangent to each round's band and meets a bevel's end in
  position.

- **Tools and provenance.** Circular, conic and bevel tools are built from
  the original seats in source-edge order. The convex planar corners that
  need setback patches are read before any boolean. Their independent tools
  remove the bands' union; section planes bound the cell replaced by the
  existing constrained G1 patch. Tool faces retain their source edges, also
  when an earlier step of the request rebuilt those edges, and the patch
  retains its source vertex through cutting and sewing.
- **Disjoint features.** Away from those corners, bevels precede circular
  rounds and asymmetric groups, with each edge and named face traced through
  the earlier engine result. Compatible tangent conics use separate capped
  tools with seated booleans, not calls back into the asymmetric constructor.
- **Periodic charts.** A partial revolved tool anchors its circle and wall
  charts at the arc's start. Its wall chart reaches past both contact levels,
  so the trim does not coincide with the chart's artificial boundary. This
  preserves the analytic pad and pocket loop volumes without asking the
  boolean to resolve a periodic seam as a material boundary.
- **Evidence.** Existing asymmetric tests are unchanged and keep their
  volume, tangency, validity, face-swapping and order oracles. New tests
  cover a bevel beside an asymmetric round, source-edge history at a mixed
  corner, a round away from that corner (which takes exactly its own
  section's volume) and the refusal of a distance of nothing. A subprocess guard reruns
  accepted corners, tangent loops, disjoint rims and named refusals with zero
  old constructor events.

### Amendment: rounds on a solid bounded by planes holding seams (#902)

A solid bounded by planes that a boolean built from boxes (an L-bracket
fused from a base and an upright) holds seams: edges between two faces
of one plane, which mainstream CAD merges after the boolean. The
chamfers read such a solid with its seams dissolved (M1, seam layer);
the rounds did not, so the planar engine read a corner of three faces
there as five, matched none of its vertex classes, and handed the
request to the current path.

- **The rule.** Rounds read the solid with its seams dissolved, as
  chamfers do: an edge asked that is a seam names nothing, the rounds
  are built on the dissolved solid, and the history is told of the solid
  as given. A solid with curved faces was read so already.
- **What it opens.** The mixed corners of §13 questions 6 and 7 on a
  fused bracket close as the engine closes them on a prism, its signed
  corner where an edge is left sharp. The current path's construction
  differs: it ends each band on its section plane and turns a shrunk
  patch round the sharp edge. The engine ends the bands on each other's
  rails, and the edge left sharp stays sharp up to the vertex. Both
  close one corner, G1 to the bands. The tests now hold any correct
  corner:
  - bit-equal orders;
  - a valid solid within the bracket and the concave round's fill;
  - every patch tangent to its neighbours but along the edge left
    sharp;
  - the patches' outline on the faces and bands;
  - the solid outside the corner's cell equal to the bands built one
    after another.
- **No refusal for a cell.** The current path refused a mixed corner
  whose rebuild box held a void, or reached through the upright's far
  face, because of how it builds the corner, not because of the
  geometry. The engine refuses a voided solid's rounds by name only where
  a round reaches a void. These corners are built (§13 question 7: no
  class of mixed convexity is refused), and the void comes through
  whole.
### Amendment: a later round reaching an earlier band at its end (#902)

Rounds asked one call after another are separate operations, each built
on what the one before it left (amendment "M5e", one request is one
operation). A later round's edge can end where an earlier round's band now
stands: the band ran on past the vertex, or ends flush in the later edge's
face. Where the band crosses the later edge's faces at an angle, the
engine already carries the later band out past that end until its section
has left the solid (amendment "M5f", "a crease carried out at both
ends"). Where the band is **tangent** to one of the later edge's faces
(it rounds the crease that face shares with the earlier edge), #1069 left
the end to the current path. That path stops the later band flush at its
own end section and keeps the corner's vertex for the corner tool
(`round_vertex`), so one edge at a time followed by the corner tool comes
out as the corner rounded in one call.

- **Survey.**
  - SOLIDWORKS' default corner solve extends the fillet faces to their
    intersection at the corner
    (https://blogs.solidworks.com/tech/2021/01/solidworks-advanced-part-modeling-partial-fillet-option-through-partial-edge.html).
    Stopping a fillet short of a face is an explicit option (Partial
    Edge, with a distance or reference offset end), not the default.
  - Onshape's default lets a fillet change the edges it meets to stay
    smooth (Allow edge overflow); keeping them is an option
    (https://cad.onshape.com/help/Content/PartStudio/fillet.htm).
  - SMLib's 3x1 corner trims the end face back to the fillet surface
    (amendment "corners where a round's rail misses").
  - No source found stops a round at its own end section so that a later
    corner tool can close the vertex.
- **The rule.** A later round reaching an earlier band at its end runs on
  past that end, as at an oblique end face, and the earlier band's surface
  trims it. The material it takes is its spandrel carried on, in the solid
  as it stands. Two convex rounds of one radius meeting so close in the
  mitre the two cylinders cut, the solid one call with both edges builds
  (§4.3). The bands are tangent to the face they share where their
  contacts cross, and the booleans meet that point tangentially; the
  intersector's tangent-contact fixes (#1054, #1080) make it build. Where
  a boolean still refuses, the request is handed back, not refused.
- **Domain.** (Widened by amendment "a later round reaching a band seated
  on a curved face".) Built where the earlier band is a cylinder seated on planes
  (every face it is tangent to is a plane). A band seated on a curved face
  (a drum's ruling rounded first on a cube the drum shaves) stays with the
  current path for now. Run on there, the rim's arc rounded after it no
  longer seats a ball at its start on the current path.
- **The corner tool afterwards.** `round_vertex` still closes such a
  corner. The current path's own corner tool already carries a band that
  ran out past the vertex back to the corner, so "one at a time, then the
  corner tool, is the corner in one call" still holds on every sheared
  block and placement in `crates/ogeom/tests/blends.rs`. Its port (the
  vertex agent, #902 tracker) inherits that contract unchanged.

### Amendment: a crease's seat past an earlier blend (#902)

Extends amendment "a later round reaching an earlier band at its end"
to the faces an earlier round cut back. A straight crease's ball must sit
on both its faces along the crease, and the way to each contact must lie
on the solid's boundary: a groove beside the crease breaks that, and such
a band is cut off, not carried on. Near an end that meets an earlier
blend, the crease's faces are cut back by that blend, and the way to the
contact stands in air there although no groove is near.

- **A blend of planes.** A cylinder or a sphere tangent along its edges
  to at least two faces, each a plane or, two levels down, itself a blend
  of planes: a band on a crease between planes, and a corner tool's sphere
  and ridge beside the bands it meets. A face tangent to fewer than two
  (a drum a band runs along) is no ball's seat and is not one.
- **Cut back, not grooved.** A point of the way in air is cut back by an
  earlier blend where, straight in from the face's plane, the first
  surface met is the own cylinder or sphere of a blend of planes bordering
  that face, solved in closed form, the point met on the solid's boundary
  and the way to it in air. Over a groove the first surface met is its
  floor or wall, neither such a blend, so a groove still cuts the band
  off.
- **Past an end face.** For a crease carried out, a point of the way in
  air on the outer side of a planar end face at one of the crease's own
  ends is past that end, where the band runs out across the face (a
  pyramid's lateral edge meeting its base obliquely).
- **Ending at a blend.** A crease end that meets a blend of planes
  tangent there to one of its faces needs no run-out past it: the band
  carried on takes what lies under that blend's surface, the two bands'
  mitre.

### Amendment: a later round reaching a band seated on a curved face (#902)

Widens amendment "a later round reaching an earlier band at its end"
(#1098) and amendment "a crease's seat past an earlier blend" to an
earlier band seated on a curved face: a drum's ruling rounded first on a
cube the drum shaves, a drum's rim rounded first beside a straight edge of
the cube's top. #1098 left these to the current path and named the cost:
run on there, the rim's arc rounded after the ruling no longer seated a
ball at its start on the current path.

- **Survey.** #1098's own: SOLIDWORKS extends the fillet faces to their
  intersection at the corner by default, Onshape lets a fillet change the
  edges it meets, SMLib's corner of one rounded edge of three trims the end
  face back to the fillet surface, and no source stops a round at its own
  end section for a later corner tool. Nothing in it depends on what the
  earlier band is seated on.
- **The rule, unchanged.** A later round reaching an earlier band tangent
  to one of its faces at its end runs on past that end, and the earlier
  band trims it: it takes its spandrel carried on, in the solid as it
  stands. A rim's band runs on so too, its run bounded as amendment "a
  band's run past its end" bounds it.
- **A band, read.** An earlier band is a cylinder, a sphere or a torus
  meeting two faces or more tangentially along its edges, whatever their
  surfaces. A crease's or a rim's seat near such a band, cut back by it,
  stands in air there and is no groove (amendment "a crease's seat past an
  earlier blend").
- **Built so far.**
  - A straight edge between a plane and a cylinder along its ruling
    (`engine/ruling_out.rs`): its spandrel is a prism, carried on past
    both ends and cut from the solid as it stands, its seat read on the
    cylinder's surface and in air where the band cut the face back. It is
    tried before the curved routes, which would stop it at the band.
  - A straight crease between planes reaching a cylinder or sphere band,
    whatever that band is seated on (amendment "a later round reaching an
    earlier band at its end"'s route, its blends read so).
  - A rim reaching an earlier band it runs tangent into: the band is no end
    face. Past that end the ball runs on seated on the rim's own faces, its
    contact on the cap in air only where the band cut the cap back, until a
    contact is in the material. The run stops where the cap's contact first
    meets the band's other seat, a plane, as it would have met that face's
    edge with the cap had the band not been there (amendment "a band's run
    past its end").
  - A straight crease reaching a rim's torus band, read as a band like a
    cylinder or a sphere. The band, the crease's band and the top touch at
    one point; the boolean cuts that since #1139. Where the crease's ball
    does not seat along its faces up to the torus (a drafted wall standing
    in the way of its contact), the torus stays the face the band ends
    flush on, as a free end (amendment "M5e"), as before.
- **Order.** Rounds asked one call after another are separate operations
  and may depend on their order (amendment "M5e"; in every application).
  Convex rounds taking their spandrels carried on take the union of their
  cuts, which no order changes. One call rounds the corner with one ball
  (§4.3); one round at a time leaves no corner for the corner tool where a
  band run on takes the corner's vertex, and does not equal the one call
  there.

### Amendment: M5d, a round resting on a far edge (#902)

The current path's rolling rounds (`rolled.rs`, `stepped.rs`,
`carried.rs`), which M6 deletes, are §13 question 8's rule applied where the
seat change covers a whole seam. Their first class moves to the engine
here.

- **The case.** A concave straight seam between planar faces whose narrow
  face is narrower than the ball's setback all along. The ball touching the
  other face rests on the narrow face's far edge, a straight convex edge,
  from one end of the seam to the other: its seat on that side is the edge
  throughout, and no crossing starts or ends the run. Its centre stands a
  radius off the touched face's plane and a radius from the edge's line,
  in closed form at each station square to the seam.
- **Past the narrow face.** Clear of the narrow face's far edge, the ball
  may rest on the far edge of a strip beyond it (a staircase), or come
  down on the face past the strips (a floor). The strips are walked out
  from the seam, and the first support proved all along is taken. The
  strips before it go into the fill whole: the fill closes on the crease
  where the touched face's plane meets the plane of the strip whose far
  edge holds the ball (or the floor's), and each strip it takes must be
  gone from the result.
- **Carried across the change.** Where the narrow face widens past the
  setback along the seam, the ball rests on its far edge up to the station
  where the face is as wide as the setback, and touches the face past it
  (§13 question 8: the face's offset and the edge's pipe are tangent
  there, so the two pieces meet on one section). Each piece is its own
  tool, the two capped on the plane square to the seam at that station.
- **Reaching a face.** Where the balls come within their radius of another
  face (a post standing on the touched face), the round keeps its surface
  (amendment "M5c"): the faces reached come off a copy, the seam is rounded
  there, and the differences are applied, the post left standing in the
  fill or beside it. The current path refused these. A face joined to those
  the run holds, or to the touched face's outline, is not reached: a ball
  that does not clear it is not this run's.
- **Beside other rounds.** Asked with edges that are no such run, sharing
  no face and no vertex with the runs, those are rounded first by the
  engine and the fills fused into what they leave. A plain concave round
  (the ball on both faces) that the rings leave because its end stands in
  the open past a corner of the face it ends on is built the same way,
  its fill ended flush on the end face's plane; along the seam both its
  touches must stand on their faces within the end faces' planes.
- **Built** as one fused tool (amendment "M5a"): the balls' arcs from the
  touched face to the support, closed by the two planes (the far edge lies
  on its strip's), capped in the arcs' planes past the seam's ends and cut
  back to the end faces' planes. Where the centre runs on a line (a floor,
  or a far edge along the seam) the round is the cylinder about it, exact;
  otherwise the arcs are fitted within 1e-8 of their exact stations.
- **Ends.** Flush on the end faces' planes, as every round of the engine
  ends (§4.3). The current path capped the fill in the end arcs' planes,
  which stand proud of, or short of, an end face not square to the ball's
  path; on the slanted test part that was 5% of the fill. Its exact
  oracle now clips each swept section to the end faces' planes.
- **Proof** (§13 question 8). At stations between the caps, wherever the
  end faces' planes keep any of the arc:
  - each ball rests on the far edge alone, its centre within the edge's
    normal cone (on its side only where the ball comes down on a face) and
    its touch within the edge;
  - its other touch lies on its face;
  - every other face stands clear by the radius, with the distance's
    Lipschitz bound up to the next station, the interval halved where the
    bound is close (a planar face wholly behind a plane the ball touches is
    clear outright; the end faces are not asked, the fill being cut to
    their planes);
  - the arcs' planes cross the seam in order.
- **A narrow face taken off whole.** Where no ball touches both of an
  edge's faces because one of them is narrower than the setback, and the
  ball touching the other face and the face beyond the narrow one's far
  edge holds the narrow face in the corner it takes off (or under the fill
  it adds), the round is that crease's round: the narrow face comes off a
  copy, its neighbours closing over it, the crease is rounded there by the
  engine, and the differences are applied to the solid as the general
  route of amendment "M5c" applies them. This is that amendment's rule (the
  round keeps its surface, the feature it reaches is trimmed), the narrow
  face trimmed to nothing. It is built only where nothing of the narrow
  face is left.
- **Not yet built.** Each of these is handed to the current path until its
  own layer: a ball carried across a change of support past the narrow
  face, an end met by other than one plane, two such seams meeting. A
  ball on a vertex is amendment "a round rolling round a vertex".

### Amendment: M5e, rounds carried on through an inside corner (#902)

Two rounded creases meeting at a corner that is reflex in plan (a notch's
inside corner, a pocket's rim) are each carried on past the corner, and
the solid keeps whatever lies under either band. A point is taken only
where it lies above both bands carried on, or above one and past the
other's wall, where the corner is open. On square walls this is the
rings' flush meeting in the bisecting plane. Off square in plan, or where
the walls lean (drafted or moulded pockets), the bands meet along their
own intersection: no flush cap or patch at the vertex holds that meeting,
and capping there left the corner standing.

- **Built** (straight convex creases between planes, where the rings
  decline) by the one boolean of amendment "M5a":
  - each crease's band carried on is its spandrel prism, the wedge of its
    two faces' planes far longer than the solid, less that wedge rounded
    by the engine (an exact cylinder);
  - each crease's tool is its spandrel where, at each end shared with
    another asked crease, it also lies past that crease's wall or within
    its spandrel;
  - the solid is cut by the tools.
- **Ends.** A crease end met by no other asked crease must carry its band
  out of the solid, so the band ends flush on the faces it meets.
  Anything else is handed back.
- **Rims.** A rim (an arc between a plane square to its axis and a
  coaxial cylinder or cone, convex: a drum's, a bore's, a D-shaped block's)
  meets a straight crease or another rim at such a corner by the same rule.
  Its band carried on is its section turned about the axis over the arc,
  and past each end the end section run straight on along the arc's
  tangent there. Past its wall is the open side of the wall's cylinder or
  cone. The section's other sides stand off the wall's surface, so no tool
  face lies on the solid's own cone.
- **A convex turn.** In a request holding an inside corner, two creases
  meeting at a corner convex in plan each run out past the vertex as they
  do alone. A request of convex turns alone is left to the routes that
  close it, in one piece where they can. Where none does, each crease
  meeting the others only at vertices of three edges, across a convex
  turn on the face they share with the third edge left sharp, is rounded
  as it is alone, and the cuts take their union (`engine/alone.rs`).
  The bands then meet along their own intersection: the mitre of two
  equal straight rounds, and the general form of it where their sections
  differ (a D-shaped block's arc and chord). This is SOLIDWORKS' default
  for fillets meeting at a corner, which extends the fillet faces until
  they intersect. A rim's round alone runs on as amendment "corners where
  a round's rail misses" sets out.
- **A free end into the material.** Where a crease end met by no other
  asked crease runs on into the material, the band ends flush on the one
  plane the crease meets there, cut to the crease's side of it; past the
  crease's middle the clip does not reach, so the band carried on past
  the other end stays that end's.
  - **The same for a rim, and on a curved end face** (2026-11-16, by
    survey). SMLib's corner of one rounded edge of three, mixed in
    convexity, says "the end face will be extended out to the fillet
    surface". Where the rounded edge shares its convexity with one
    neighbour, that edge is extended to meet the rail too. Either way the
    band ends flush on the end face's surface, carried on past the face.
    The one face a crease or a rim's arc runs into at such an end is that
    end face, whatever its surface: a plane, a cylinder, a cone or a torus
    (`engine/flush.rs`), an earlier round's band among them. A rim's arc
    also ends on a band it runs tangent to; a straight crease ends only on
    a face standing across its way, and runs out past a band it merely
    touches (a block's corner rounded one edge at a time). Near the end, the band loses what lies on the far side
    of that surface from the crease. A lone straight crease with such an
    end goes through the carried-out route, carried only a few radii past
    that end. A rim goes through the run-on rims' route.
  - **A rim's band takes nothing past its wall.** The rim's tool stands
    off its faces' surfaces, so that no tool face lies on the solid's own.
    Where the solid holds material just past the wall near an end of the
    arc (a block leaning against the drum, at an inside corner), the tool
    loses that material: its standoff takes nothing the solid holds beyond
    the wall.
  - **Figures that change.** The current path capped a rim's band in the
    plane through the axis at such an end, so the keyhole's rounds lost
    the whole corner over their arcs (`chain.rs`). Rounded one edge at a
    time, the drafted drum on a leaning block takes a corner in each order
    by that order's own rule. Separate rounds depend on their order, as in
    every application.
  - **A crease past a rim's flush end.** A rim's band ending flush on a
    plane leaves a piece of that plane beside the face it ended on; the
    route merges the two, and the crease along them is one edge again, to
    its real end on the band. A crease ending at an earlier band it runs
    tangent to ends there only while its line carried on stays on the
    solid (one of its faces runs on past the band); where the line leaves
    the solid over the band, it runs out as #1098 says. Past a torus band
    from outside lies its tube and, within the tube's height, the ring's
    hole beyond the tube's centre circle: a crease carried on through the
    tube does not come out again on the other side (`engine/flush.rs`).
- **Not built.** Two bands carried on that would pinch the solid where they
  meet (a crease turning back on itself, its notch and the bands' overlap
  past the corner meeting on the corner's line) fail the solid's guard and
  are handed back.
- **Refused by name: zero thickness.** A rim's band run on past its
  corner whose own figure leaves the solid touching itself along an edge
  (an earlier band and this one meeting their walls at one depth, the gap
  between them closing to the edge, material on both sides) is refused by
  name. SOLIDWORKS does not allow such zero-thickness geometry either
  ("Zero Thickness Geometry": edges or vertices that do not join the
  adjacent geometry, non-manifold; help.solidworks.com, Error Message -
  Zero Thickness Geometry). `guards::touches_itself` reads it: two faces
  back to back at an edge that part away from it with void between them.
  Faces folded flat on one another (no gap) are still handed back.

### Amendment: M5f, a strip's end cut by its end faces' planes (#902)

A strip of the rounds' sheet ends flush in the third face of a vertex of
three (§4.3), cut by that face's plane. Two more ends are cut the same way,
by the end faces' planes, as mainstream CAD caps a blend on the faces
around its end.

- **Two faces beyond, meeting on an edge left sharp.** One round at a
  vertex of four faces (a pyramid's base corner, its ridge rounded): the
  strip carried on is cut by the plane of the face beyond each of its own
  faces, from its corner on the face beside it to where the two cuts meet
  on the sharp edge, past the vertex. The sharp edge's end comes off with
  the strip, and each face beyond ends on its cut. Each cut is the
  section's exact plane end, trimmed where it crosses the other plane,
  found once on it. Built for a convex round beside a convex sharp edge,
  each cut within its face and the meeting point within the sharp edge.
- **Past the face, under a roof.** A leaning end face whose plane runs on
  past its edge through the round's material (a rib running up under a
  roof): the strip's end there stands in that plane, which caps it, the
  face carried on over the cap. Each point of the end must lie in the face,
  or in its plane with the solid's material on both sides of it at the
  vertex, read from the vertex's local cone.
- **Local cones at reflex corners.** The vertex's local cone reads each
  face's sector the wide way round where the face's corner there turns
  past a half turn, from the turn of its ring against the ring's winding.
- **A crease carried out at both ends.** A straight convex crease between
  planes that no other route closes (a chord's crease on a D-shaped block,
  its ends on the drum, square, drafted or overhanging), whose ball sits on
  both faces along it and whose band carried on past each end leaves the
  solid: the band ends flush on the faces it meets there, as a free end
  of amendment "M5e" does. Each end's faces meet the crease's at an
  angle; a face running on tangent from them (a band an earlier round
  left) is that round's to meet. Its tool is its spandrel carried on, which must
  meet the solid in one piece, so it takes the band and nothing farther
  along its line. Asked edges meeting or sharing a face are left to the
  routes that check their bands against each other.
- **A rim asked by one of its arcs.** Where a boolean leaves a closed rim
  as several arcs and the request asks some of them at one radius, the
  rim's other arcs join at that radius for the curved routes: one value
  rounds a tangent chain (amendment "what M6 refuses by name"), and the
  current path rounded the whole turn where the rim closes. A rim that
  stops, or arcs asked at different radii, are as asked.
- **Unproved faces.** A round's face whose proof within its material
  neither clears nor finds a point past it (a patch tangent to a face only
  to its fit) is handed back, not refused. One found standing past the
  material is still refused.

### Amendment: vertex-only rounding is a local corner feature

The vertex APIs are not shorthand for an edge selection. This distinction
matters when the last old constructors are removed: rounding all incident
edges removes entire strips that `round_vertex` must leave sharp.

**Public documentation surveyed.** No implementation source was read.

| Source | Evidence and limits |
|---|---|
| Blender 3.6, [Bevel Edges](https://docs.blender.org/manual/en/3.6/modeling/meshes/editing/edge/bevel.html), Affect: Vertices | Explicit vertex-only operation: “Only the areas near vertices are beveled, the edges remain unchanged.” It distinguishes vertex and edge selection. Its tessellated profiles are not an oracle for an exact B-rep sphere or for the size of its caps. |
| Fusion, [Fillet reference](https://help.autodesk.com/cloudhelp/ENU/Fusion-Model/files/SLD-REF-FILLET.htm) | Solid fillets select edges, faces or features. Rolling Ball and Setback describe corners where fillets meet; there is no documented standalone solid-vertex selection. This is not evidence that vertex-only rounding is refused. |
| Creo, [Round Transition Types](https://support.ptc.com/help/creo/creo_pma/r9.0/usascii/part_modeling/part_modeling/Round_Transition_Types.html) | Corner Sphere applies where three round pieces overlap. Setback distances modify that transition. Stop at Reference terminates at a specified datum point or plane. None defines rounding every incident edge merely by selecting a vertex. |
| Papers [5–7], and Salvi's [public abstract](https://arxiv.org/abs/2002.11347) | The setback literature separates the edge bands from the vertex opening. The abstract describes boundary interpolation from positional information, not a G1 guarantee or a vendor's vertex-only default. The papers and existing surveys therefore do not settle the standalone feature's cap placement or its refusal domain. |

The public ACIS documentation mirror could not be consulted because its
certificate did not match its hostname. The publisher's paper page returned
403. Neither supplies evidence for an additional refusal. The earlier
vendor survey in §13 still governs edge-round transitions; it must not be
extended to this different selection mode by assumption.

**Rule and compatibility boundary.**

- `round_vertex` changes only the corner. Where one ball fits, the corner
  is its spherical polygon. Each unrounded incident edge ends on its
  plane through the ball centre square to the edge, leaving a planar cap.
  These caps are intentional G0 stops, not failed G1 joins. Existing bands
  consume their caps and meet the same spherical polygon at G1. A band
  stopped by a different feature is not extended across that feature.
- Where no single ball touches all planar hosts, retain the exact envelope
  of feasible ball centres: spherical pieces at the offset region's tips
  and cylindrical pieces along its ridges. Do not replace it with an
  all-edge request's unconstrained vertex patch. Three curved hosts seat
  their ball by their surface offsets; radial end planes close that local
  spherical corner against the original solid.
- `round_vertex_setback` explicitly requests the incident edge bands as
  well as the corner. It retains one generated filling face over its
  opening; `round_vertex_setback_sharp` retains the point patch at the
  sharp edge and one filling face over the remainder. Their caps, trims,
  requested distances, G1 checks and named G2 refusal are unchanged.
  Both use the engine's sections, face trims and assembly. The common
  filling and point-patch fitters remain geometric primitives, not the
  old cut-a-cell-and-resew corner constructors.
- Every newly exposed face, including a vertex-only planar cap, is
  generated by its source vertex. Inherited host faces remain inherited.
  Composition must retain distinct ownership when several vertices are
  rounded in turn. No new refusal follows from an undocumented vendor
  feature; existing domain and foreign-feature safety checks remain.

The independent box oracle is local: the volume removed at one square
vertex is `r³ (1 − π/6)`, not the volume removed by three full edge bands.
The sphere has centre a radius from each host plane, and its three caps
lie on the planes a radius along the three edges. These geometry and
ancestry contracts survive replacement of the construction.

### Amendment: corners where a round's rail misses (#902)

Two vertex classes the earlier amendments left open, settled by survey
(2026-11-15). The SOLIDWORKS, Fusion, Onshape and NX documents describe
the options at a corner (setbacks, rounded or smoothed corners, rolling
ball or setback corner types, overflow handling), not what a round does
where its rail misses the edge it should meet. The SMLib manual's
fillet corners do, classing each vertex by how many of its edges are
rounded and by the convexity of all its edges:

- one rounded edge of three, all of one convexity: "the end face will be
  trimmed back to the fillet surface";
- two rounded of three, of one convexity, the third of the other: "the
  rail curves that are in the same face must be extended and/or blended
  to meet", and a corner patch is created.

Sources: SMLib, *Fillets* (docs.nvidia.com/smlib/manual/smlib/fillets);
SOLIDWORKS, *Constant Size Fillets*; Onshape, *Fillet*; Autodesk Fusion,
*Create a fillet*; Siemens NX Open, `EdgeBlendBuilder`. One request
rounds its edges as one operation, so the order its edges are asked in
changes nothing (as in every application, where only separate fillet
features depend on their order).

- **Two rounds at an inside corner whose rails meet past the face.** Two
  rounds of one sense meeting where the edge left sharp at their vertex
  is of the other sense (an inside corner): their rails on the face both
  share are extended to meet. Where they meet within that face, the
  bands meet along their own intersection, as amendment "M5e" (#1062)
  builds them. Where the extended rails meet only past the face (a crease
  turning back nearly on itself, its notch narrower than the bands), the
  rails are blended to meet on that face and the corner closed by the
  ball rolling round the sharp edge (section 13, question 8: the ball
  rests on the edge):
  - each band ends in the plane through the sharp edge square to its own
    crease, where its ball touches that edge;
  - between the two, the ball rolls round the sharp edge on the
    material's side, touching it and the shared face: the corner is the
    bands' section turned about the edge through the angle between the
    two end planes behind the vertex, G1 to both bands (all three sweep
    one section square to its plane where they meet);
  - the blended rail is the circle of the radius about the vertex on the
    shared face, tangent to both rails.
  Built for rounds of one radius whose sharp edge stands square to the
  shared face, the turned section's surface exact (a horn torus about the
  edge); other such corners are handed back. The bands are not carried on
  into one another there: carried on, their overlap past the corner and
  the notch meet on the corner's line, which pinches the solid.
- **One round whose rail crosses the adjacent edge away from the vertex.**
  One rounded edge at a vertex of three edges of one convexity (a boss's
  rim where a rib's wall leaves it, turned towards the material). SMLib's
  case has the end face trimmed back to the fillet surface where the
  rails meet the edges beside them at the vertex. Here the rail on the
  cap meets the edge it shares with the end face only well back from the
  vertex, and the boss's circle carried on can run back into the rib
  before it does (where the rib holds the boss's centre), so the case is
  extended:
  - the round's surface runs on past the vertex along its own spine (a
    rim's torus about its axis) only while the end face trims it: while
    its contact on the cap stands in the solid and its contact on the
    wall's surface carried on stands outside it;
  - it stops at the first station where that trim reaches either
    contact (the cap's contact leaving the solid, or the wall's contact
    entering it), and where it stops inside the material it ends in the
    plane through the axis there;
  - over the stretch it runs on, it takes what of its section the solid
    holds, and the end face is trimmed back to it.
  It is not capped in the plane through the vertex, which leaves material
  the ball takes standing, nor run on into the material past where the
  wall's contact enters it, which would cut a groove with a cliff on the
  wall's surface carried on.
- **Not a rim that stops.** An arc whose circle's arcs on the solid run
  the whole turn (a rim a boolean's seams split) has no end: one radius
  rounds the tangent chain ("what M6 refuses by name"), so the rim is
  rounded whole, as the curved routes round it with its seams dissolved,
  and not by this rule, which would leave the arcs not asked sharp.

### Amendment: a band's run past its end (#902)

Amendment "corners where a round's rail misses" runs a rim's band on
past its vertex while its contact on the cap stands in the solid and its
contact on the wall's surface carried on stands outside it. Amendment
"M5f" and the oblique rim arc run a band on "until its section has left
the solid". On a drum cutting a cube's corners the two part: the cap's
contact circle touches the cube's side at one point, a quarter turn from
the corner, and stays on the cap past it. Read by contacts, the band runs
on along the whole side to the next corner, taking a sliver that meets
the corner's own only at that point. Read by sections, it stops there.

**Survey.** SMLib's manual (*Fillets*) trims a surface fillet "to the
first intersection with a boundary of one of the base surfaces" (minimal
trim), and at one rounded edge of three of one convexity "the end face
will be trimmed back to the fillet surface". SOLIDWORKS' "Keep surface"
and the Parasolid family's notch trim the face a round overflows to it
(amendment "M5c, non-local rounds keep their surfaces"). None runs a
round on past the first boundary its contacts meet.

**Rule.** Past an end with no face standing across its way, the band runs
on along its spine and the end face is trimmed back to it, up to the
first station where either contact meets a boundary:
- the cap's contact meets an edge of the cap, crossing it (the contact
  leaves the solid) or touching it (the section's material shrinks to a
  point there: the section has left the solid);
- the wall's contact, carried on, enters the solid.
That station is found exactly: the cap's contact circle is intersected
with the cap's edges, touching points included. At a crossing the tool
runs half a step on, holding nothing; elsewhere it ends in the plane
through the axis there. This is one rule for both readings: where the
contacts cross the boundaries, the sections leave the solid at the same
station; where the cap's contact only touches, the section leaves and the
contacts alone do not see it.

### Amendment: a corner patch held within a convex corner (#902)

A corner's transfinite patch is interpolated on a net and held to the
bands' fit target against its sides. Between the sides it can lean past a
face it is tangent to along its side: a pyramid's apex at alternating radii
stood 4.7e-3 past the face, against a fit of 2e-4. P5 does not prove such a
patch, and declines it.

Over a convex corner, whose material is the half-spaces of its faces'
planes, a patch P5 did not prove has its net held within those planes
(`ogeom_algo::net_held_within`; docs/research/vertex-patches.md, §B):
- **Kept rows.** The net's two outermost rows each way are kept as they are.
  They carry the patch's sides and its derivatives across them, so its
  tangency to the bands holds by construction.
- **Moved points.** Every other point standing past a plane moves to the
  corner's nearest point. That point is the one minimum of the squared
  distance under the planes' inequalities, found by trying the planes'
  active sets in a fixed order over the planes sorted. The held net is the
  same to the bit whatever order the corner's edges are asked in.
- **Proof.** A B-spline with positive weights lies in the hull of its net,
  so the patch lies within the corner. P5 proves it again on the stored
  surface, as the independent guard. A bound it cannot settle is not read as
  inside.
- **Infeasible.** Where a kept point stands past a plane by more than the
  fit, no such net exists. The patch is not held and is declined as before;
  scaling the cross-derivative rows is not tried.

### Amendment: rounds either side of a cusp (#902)

Where a flat side runs tangent into a drum's wall (a plate's side into a
boss), the top between them narrows to a cusp at the tangent point. The
top's rim on the wall (the boss's foot) and the top's edge with the side
continue each other through the cusp, one concave and one convex. Where an
asked rim and an asked straight edge leave the cusp on the same side, the
two bands turn through an S, the trimmed assembly's cusp patch. This
amendment covers the requests where they do not: the foot alone, or the
foot on one side of the cusp and the edge on the other.

| Source | What it says |
|---|---|
| SOLIDWORKS, Fillet PropertyManager, "Overflow Type" | Where a fillet runs into adjacent geometry, "Keep surface" trims the fillet by the adjacent face; "Keep edge" keeps the adjacent edges and splits the fillet. [help.solidworks.com, Example: Overflow Type] |
| Fusion, Tangent chain reference; fillet best practices | A tangent chain is selected as one set; fillets that pass from concave to convex edges are best made as separate fillets, concave first. [help.autodesk.com, SLD-REF-TANGENT-CHAIN; Autodesk blog, Fillet Best Practices] |
| Onshape, Fillet help and forum | Tangent propagation extends a fillet along tangent edges; fillets passing from convex to concave edges are reported as unreliable. [cad.onshape.com/help, Fillet] |

No application documents a blend between a convex and a concave band of one
tangent chain other than as separate fillets. Where a band runs past the
face it stands on, the documented choice that keeps the round whole is to
trim it by the adjacent face. So:

- **Rule.** Each band is built as it is alone. A rim's band runs over its
  arc and ends in the meridian planes through the arc's ends. Where it runs
  past the side at a cusp, the side's plane trims it, and its fill ends
  flush with the side. A straight edge's band ends in the plane square to
  the edge through the cusp, as a band ends where a round of another kind
  continues it tangentially, and at its far end in its end face's plane.
  The request builds the same solid in any order.
- **Built** as the one boolean's tools (amendment "M5a"). A rim's tool is
  its section (the corner between the cap and the wall less the ball's
  disc, or the ball's fill) turned about the drum's axis over its arc. A
  straight edge's tool is its section swept along the edge. Only the
  round's own face (the torus or the cylinder) is named for its edge. The
  tools are applied in the solid's order of its edges, so the request's
  order does not change the result.
- **Refused by name:** a rim and a straight edge leaving a cusp the same way
  on one side of it while the top runs on past the cusp with its edges
  there left sharp. The S has nothing to close against on the cusp's other
  side, and the request is refused with a request to round the chain on
  both sides (the trimmed assembly's S), as before. Superseded where the
  trimmed assembly builds it: amendment "a cusp's S on one side while the
  top runs on".
- **Not this route:** a rim and a straight edge leaving a cusp the same way
  where the top ends at the cusp (the one-sided S); a rim on a cone; any
  other edge in the request.

### Amendment: a rim arc ending obliquely on a flat (#902)

A drum cut flat across its rim leaves the cap's rim an arc whose ends meet
the flat obliquely: the flat's plane neither holds the rim nor stands
square to it there. Rounded alone, where does the round end?

- **Rule.** The ball rolls on past each end through the corner it still
  finds in the solid, while its contact on the cap stands on the cap, and
  the flat's plane trims it. No sliver of sharp crease is left between the
  cap and the flat. This is the trim by the adjacent face that SOLIDWORKS
  documents as its "Keep surface" overflow [help.solidworks.com, Example:
  Overflow Type]. It matches how a straight round ends on its end faces'
  planes (amendment "M5f").
- **Built** as one tool of the one boolean (amendment "M5a"): the rim's
  section (the corner between the cap and the wall less the ball's disc,
  or the ball's fill) turned a whole turn about the drum's axis, cut by the
  plane of each face the arc ends on. Only its torus is named for the edge.
- **Also built so:** an arc in pieces of one circle asked together (a
  wall's seam splits a drum's foot), a concave arc (a drum standing on a
  block's top over its side, its foot filled and the fill cut back to the
  side), and an end on a plane square to the arc (a meridian plane: the
  whole turn cut there ends the round flush in it, as §4.3 ends a strip in
  its third face). The cap and the wall must share no edge but the arc's
  pieces: the whole turn rounds every edge they share.
- **Not this route:** an arc ending on more than one face at an end or on
  a curved face, and an arc asked with other edges.

### Amendment: contour chains on swept walls use engine-owned tools

A contour here is a chain of constant-radius creases on a planar host
whose other faces are straight walls, circular walls or walls swept from
fitted profiles. A tangent continuation may cross a split between planar
hosts. The rules are the existing end-trimming and marched-junction
rules, which mainstream systems share: a straight run under an
overhanging planar roof ends on that roof's plane; a split rim that stops
is restricted to its own arc; a fitted crease meeting its contour at a
re-entrant turn includes the junction between the two bands. No new
vertex type is introduced and nothing beyond the request is rounded.

The route is the last the engine tries. A rim whose band runs on past a
free end (amendment "corners where a round's rail misses") and every
other route keep the requests they take, so the route only receives
requests no engine route closed before it.

The engine builds every section on the original solid, using the whole
request for tangent stops and junction planes. The exact planar section,
revolved arc and marched band supply unapplied tools; no contour adapter
calls the old edge dispatcher, intercepts an old constructor, or applies
one edge before constructing the next. The tools retain their source-edge
credits when their adjoining caps are joined. Topology, fitted-contact
error, orientation and generated history remain independently checked.

The rib's straight and curved sections already have independent section
area, centroid and turn-moment oracles. For the radius-one fillet on the
closed pocket floor, the added volume is `L A - 2 pi M`, where `L = 18 pi`
is the perimeter, `A = 1 - pi/4` the corner-section area and
`M = 5/6 - pi/4` its first moment. The pocket's unrounded volume is
`30*30*10 - 5*(45 pi)`: a five-deep 120-degree annular sector of radii
6 and 12, with two radius-three semicircular ends. This independently
derived oracle replaces the recorded rounded-volume figure without
changing its tolerance.

### Amendment: rolling corners (#902)

`fillet_edges` closes a vertex where convex and concave rounds meet on one
G1 patch between the bands' ends (§13, questions 2 and 7). Users of other
systems also expect the corner the rolling ball itself makes, with no
patch. This amendment adds it as an option, `FilletCorner::Rolling`, through
`fillet_edges_corner`. `FilletCorner::Setback`, the default, is
`fillet_edges` unchanged.

**What the systems offer.**

| System | Behaviour |
|---|---|
| Autodesk Fusion [P44] | A corner type per fillet feature, **rolling ball** or **setback**. Asymmetric fillets always take the setback. |
| SOLIDWORKS [P45] | **Round corners**: constant-size fillets with "a smooth transition between edges", the ball's own corner. Setbacks are a separate, opt-in section. |
| SMLib [14, S1] | Equal radii at a vertex: the ball's sphere where it touches all three faces, and rolling-ball surfaces where a ball rolls round an edge's end. Coons patches only where neither fits. |

**Rule.** With rolling corners asked, every vertex where two or more of the
rounds meet is closed by the ball rolling through it. Every new face is a
plane, a cylinder, a torus or a sphere, and every contact is G1. The result
does not depend on the request's order.

- **Three rounds of one sense:** the ball's sphere, as by default.
- **Two convex and one concave, or two concave and one convex** (the lone
  round), at a vertex of three planes with every edge rounded at one
  radius `r`:
  - The lone band runs unbroken to where its ball touches the face across
    from it (the cap).
  - The other two bands end square at their corners on their other faces.
    Between those ends their ball rolls round the lone band's end, touching
    the cap and the lone band, so its centre runs on the circle of radius
    `2 r` about the lone band's axis, in the plane of the bands' centres.
  - The corner face is the torus that ball sweeps: major radius `2 r`,
    minor radius `r`, its axis the lone band's axis.
  - It meets the lone band along that band's end arc (the torus's inner
    equator), each straight band along a meridian, and the cap along its
    contact circle of radius `2 r`.
  - This is a torus only where the lone round's edge stands square to the
    cap, so the bands' centres stay `2 r` off its axis on a circle.
    Elsewhere the sweep is a canal surface on an ellipse, which is not
    analytic, and the vertex is refused.
- **Two rounds of one sense beside an edge left sharp:** the bands meet
  along the curve where they cross, as by default.
- **Refused by name, naming the vertex**, never closed by setbacks
  instead:
  - a vertex where convex and concave rounds meet that touches a curved
    face, or joins other than three edges (where the rounds are of one
    sense the default's corner is taken: amendment "rolling corners on
    curved faces");
  - radii that differ;
  - a convex and a concave round beside an edge left sharp;
  - a lone round's edge not square to the cap;
  - any vertex the default would close on a patch.

  A request that lies outside the engine's rounds between planes (a round
  reaching a face beyond its ends, a face with a hole) is refused too.

**Oracle.** With every edge of a prism rounded, the rolled body is the
opening of the prism's closing. Its volume follows from the cap's offset
region by slices. The printCAD stairs (three blocks, ten faces, 24 edges,
`r = 1`) give `4966 + 157π/3 + π²`. That equals the blocks less the
straight spandrels, `(1 − π/4)·r²·L` each, the twelve octant corners at
`r³(1 − π/6)` each, and the four elbows, each giving back
`r³(1 − 7π/6 + π²/4)`.
A boss 6 square and 6 high on a plate 20 square and 4 thick, rounded at
`r = 1` on its upright edges and its foot, gives `1808 + 17π/3 − π²`. There
each elbow takes `r³(1 − 7π/6 + π²/4)`, the same term with its sign turned,
as rounding the complement turns it.

- [P44] Autodesk Fusion API, `FilletFeature.isRollingBallCorner`:
  https://help.autodesk.com/cloudhelp/ENU/Fusion-360-API/files/FilletFeature_isRollingBallCorner.htm
- [P45] SOLIDWORKS, Constant Size Fillets (Round corners):
  https://help.solidworks.com/2019/english/SolidWorks/sldworks/r_constant_size_fillets.htm
- [S1] NVIDIA, SMLib Fillets: https://docs.nvidia.com/smlib/manual/smlib/fillets/

### Amendment: a cusp's S, one side or two (#902)

Where a flat side runs tangent into a drum's wall, the top between them
narrows to a cusp at the tangent point T. A round of the top's edge with the
side (convex) and of its rim on the wall (concave), leaving T on one side,
cannot stand side by side within two radii of T: the round there is an S,
off the side in a convex arc and onto the wall in a concave one. Both sides
of T may be asked (a plate's side tangent to a boss: two S's) or one, where
the top ends at T (an arm's side running into a hub: one S).

- **Lobes.** Each side of T asked with both rounds is one lobe. The two
  bands end in the meridian plane through the point where their contacts
  cross on the top, each cut by it, so the lobe turns smoothly from the one
  band into the other.
- **Two sides: one patch over T.** The hole is the two lobes and, between
  them, the rim's contact on the wall and the straight strip's contact on
  the side; one patch closes it, standing clear of T, where the side and
  the wall are tangent. On a plate's side tangent to a boss the rounded
  solid is tangent within 3e-5 at every radius from a twenty-fourth of the
  wall's to a quarter, against the 2e-4 target. A patch per lobe
  over the S's feet (below) has its pole at T, where at small radii the two
  feet leave T almost parallel: there it held 1e-4 and missed the target at
  a twenty-fourth (2.2e-4), finer nets diverging. The two-sided cusp keeps
  the patch over T (user, 2026-11-17).
- **One side: the S's feet.** Where the top ends at T there is no second
  lobe to close against, and the hole is the lobe and the S's feet down to
  T. Below the lobe the S is read section by section square to the
  straight edge: where the top stands open `u` of what it does at the
  lobe's end on a face (open: from the side's line out to the wall's foot
  along the radius, which reads at every station, also where the lobe's
  end on the side lies past the foot's reach along the top, as it does
  for the wider rounds), the S turns `t` off each face, `cos t = 1 - u`,
  its foot on the side `sin t` of the lobe's depth below the top and on
  the wall `sin t` of its height above. Both feet run from the lobe to T,
  each meeting its band's contact tangentially at the lobe. The patch over
  the lobe and its feet has its pole at T; the lobe runs from a convex band
  into a concave one, which one polynomial patch cannot follow, so the
  patch is fitted on each side of the line from that join to the pole and
  joined there C1 on a knot line, one face. Edges the S meets at T (the
  side's edge with the wall below the top, the wall's seam) run on to T.
- **Held** to the bands' fit target, gap and tangency, either way.
- **Refused by name:** one lobe asked while the top runs on past T with its
  edges left sharp (amendment "rounds either side of a cusp"). Superseded:
  that case is built, see amendment "a cusp's S on one side while the top
  runs on".

### Amendment: a round rolling round a vertex (#902)

Where the narrow face is a triangle over the seam, the ball resting on its
far edge reaches the triangle's corner before the seam's end. Past it the
ball touches the corner alone: its centre runs on the circle where the
plane a radius over the touched face meets the sphere a radius about the
corner, tangent there to the path it had on the edge (§13 question 8: the
edge's pipe and the corner's sphere are tangent where the edge's touch
reaches the corner).

- **The case.** Amendment "M5d"'s run, where the touch on the far edge
  reaches an end of the edge inside the seam. The ball on the edge, then
  on the corner, are two pieces meeting on the arc's plane where the touch
  is at the corner; the corner's piece runs on to the seam's end.
- **Built** as M5d's tools, one a piece. The corner's band is the ball's
  arc from the touched face to the corner turned about the line through
  the corner square to that face, exact: a spindle torus, every arc ending
  on the corner. The edge's piece is fitted as M5d's are. Each fill is
  closed by the touched face's plane, the narrow face's plane from the
  crease to the support, and the arcs' planes at its ends. Where a
  piece runs to the seam's end, the arc's plane caps it where the plane
  crosses the seam's end.
- **Ends.** Not cut to the end faces' planes. Both end faces carry the
  support (the corner lies on each, and the far edge in one), so their
  planes pass through the balls' touch on it and would cut into the fill
  the support itself holds. The fill is each piece's sections whole, to
  the arc's plane through the seam's end, and the oracle sums those
  sections unclipped. Sampled point by point against the rolling balls
  (each point in its nearest centre's normal plane, outside that ball,
  between the seam's two faces and on the corner side of the chord from
  the touch to the support), the material past the end faces' planes is
  that region and nothing else. Cutting flush (the boolean does it on the
  exact band) takes 28% of the fill, the edge piece's whole touch on the
  touched face with it, which lies past the plane of the end face holding
  the far edge: its round would no longer meet that face tangentially.
- **Proof** as M5d's, the corner standing for the far edge, at every
  station whether or not the end faces' planes keep it: each ball's touch
  on its face lies on the face, no face of the solid stands nearer than
  the radius (the end faces and the narrow face asked too, so the ball
  meets the faces at the corner only at the corner), and the arcs' planes
  cross the seam in order. The edge's piece is proved on the edge, within
  its normal cone, up to the switch.
- **Not this route:** a corner where the circle and the edge's path do not
  meet within the seam, a corner on a curved face, a piece whose fill
  shell does not close.

### Amendment: inside-corner meetings on a solid holding seams

Rounds meeting at an inside corner (amendment "M5e") read the face both
creases share. On an imported solid with curved faces elsewhere, a seam
can split that face at the corner into two faces on one plane, and the
creases then share none. As for a solid bounded by planes holding seams,
the request is read again on a copy of the solid with its seams dissolved,
and the history is told of the solid as given. Nothing else changes:
where the solid holds no seam, or the dissolved solid does not close the
corner either, the request goes on as before.

### Amendment: edges split into collinear pieces (#902)

Three blocks stepping up, fused and unified, keep each straight side
edge of the step face as three collinear pieces: the unification merged
the faces, not the edges. The engine read such a solid as given only
where a seam ran into it, so with no seam left it declined, and the
current path rounded each piece as its own band and closed each inside
corner with its signed-corner box. Each side edge became three
cylinders with seams between them, and each inside corner kept the
parts of three faces inside its box as faces of their own, their
outline standing as edges.

- **The rule.** The engine reads a solid bounded by planes with its
  seams dissolved **and** its runs of edges on one curve joined, whether
  or not any face merged: a vertex only two pieces of one straight edge
  between the same two faces meet at is no vertex of the part's form.
  A round or chamfer of any piece is one of the whole edge, as tangent
  propagation carries a pick along its smooth continuations by default
  (SOLIDWORKS in the M6 refusals amendment's survey, Onshape in the
  amendment on varying rounds on curved seats). Pieces asked at
  different radii are declined there, as for pieces a seam split.
- **What it fixes.** The stepped block's every-edge round goes through
  the engine: one cylinder per straight edge, one sphere per convex
  corner, one G1 patch per inside corner (§13 question 7), and every
  original face one face. The setback assembly leaves no same-carrier
  pieces beside their parent faces on this class. The signed-corner box
  of the current path is not reached.
- **Unification joins edges by default.** `ogeom_heal::unify_same_domain`
  now joins each run of edges on one curve after merging faces, as
  `merge_edges` does, and counts both; its history modifies each face
  of a merged set into the merged face. Public documentation of
  mainstream refinement does the same: in CadQuery's `clean`, "adjacent
  faces and edges lying on the same underlying geometry are merged", and
  it runs after every boolean unless `clean=False` [P46]; build123d's
  `clean` returns the shape "with extraneous internal edges removed"
  [P47]. `unify_same_domain_with` and `UnifyOptions { edges: false }`
  keep the edges split; `unify_same_domain_around` merges faces only,
  since an edge it joined is another face's boundary too.
- **Pieces on one carrier merged after any route.** Other routes do
  leave them: strips applied by a boolean, rolling bevels meeting in
  corner pieces, a slab rounded at half its thickness whose top and
  bottom bands are one half cylinder. Every round and chamfer the engine
  builds is followed by one merge (`engine::merged`): the faces the
  blend made, unified with their neighbours on one carrier
  (`unify_same_domain_around`), the result checked, and the merge's
  history composed after the blend's. A face taken into another is
  modified into the merged face, so an edge still names the face its
  blend became. A cluster the unification leaves split (a band closing
  round its axis with no seam) stays split; splits between faces the
  blend left alone stay as the input held them.
- **A cusp is a pinch, not a smooth turn.** A concave band's end lies on
  the face it ends on, so the merge joins the two: an L-bracket's end
  face then holds the band's arc between the leg's top edge and the
  wall's edge, each tangent to the arc where they meet. The route for
  rounds either side of a cusp read those tangent points as cusps (a
  flat side tangent to a drum's wall) and built the rim's round as its
  separate tools, cut off at the "cusp". A cusp is where the top narrows
  to nothing: the straight edge and the rim's arc leave the tangent point
  the same way. Where they leave it opposite ways, the top runs on
  straight across it, the chain turns smoothly, and the route does not
  apply.
- **A cusp edge is not a fold.** P5 refused two faces meeting along an
  edge facing opposite ways as folded flat. Two elementary surfaces of
  different kinds (a rim's torus and a drum's cone) cannot lie on one
  another; meeting so, they are tangent along the edge and the material
  narrows to it, as a rim's band carried round past a crease's flush end
  leaves it (the crease's cone piece beside the torus). The fold guard
  holds where the two can coincide: one surface kind, or a B-spline or
  other non-elementary surface. The merge above joins the drum's cone
  that the crease's end split, and the rim's round carried past the
  corner, refused before for that edge and built by the current path,
  goes through the engine's exact route to the same volume.

### Amendment: rounds standing apart (#902)

A request can hold edges that no route closes together, though each group
of them closes alone: a boss's base rim and a rib's base creases on one
plate. The rib's creases end flush in its end faces' planes, beyond the
plate's top, and the curved assembly does not carry an end face over a
fill's end section, while the rounds asked alone are built.

- **Rule.** The edges fall into groups, two edges in one group where a
  chain of edges sharing vertices joins them. Where no route closes the
  request whole, each group is rounded as it is alone on the solid as
  given. What its rounds take off and what they add are its tools. Where
  every tool stands clear of every other group's (a gap wider than the
  confusion tolerance), the tools are applied as the one-boolean
  fallback's (`engine/apart.rs`).
- **Why it is the rounding asked.** A band built alone meets only the
  faces within its tool's reach. Where no other group's tool comes within
  that reach, no band of one group meets a band of another or a face
  another group changes, so the rounds asked together are the rounds
  alone. This is amendment "M5e"'s convex-turn rule (each crease rounded
  as it is alone, the cuts taking their union) for groups that do not
  meet at all, with fills as well as cuts. No new behaviour is defined:
  each group's result is the one the engine already builds for it.
- **In turn.** A group's tools are read by booleans between the solid
  and its rounded copy, which share every face but the round's; on a part
  whose faces the booleans cannot match (two bores crossing), they refuse.
  The groups are then rounded in turn, each on the result of those
  before, those holding a marched seat first (the march reads the solid
  as given more surely than a rebuilt one), and a group the engine does
  not round there is moved first, once. The rounds are those standing
  apart only where each later group's edges come through the earlier
  rounds as one edge on the same ends, and each face an earlier round
  made comes through the later ones with its edges as they were (a band
  that met it would cut a new edge into it); otherwise the route
  declines.
- **Not this route:** groups whose tools touch or overlap (their bands
  meet and the meeting is the routes' to close), a group the engine does
  not round alone, or a request of one group.

### Amendment: a straight crease run out across a curved face (#902)

Amendment "M5e" carries a straight crease's band out past an end that
runs out of the solid, the band ending flush on the faces it meets there.
Its seat check excused a contact in air only past a planar end face. A
curved end face can stand over the crease's faces before the crease
ends: a drum cut across a block's bottom edge overhangs the back wall, so
the edge's contact on the wall runs off the wall onto the drum short of
the edge's end.

- **Rule.** Per section 13, question 8, the contact that runs off its
  face onto the curved face the crease ends at ends there: the band runs
  on straight and the curved face trims it. The seat check excuses a
  contact in air past such a face as it excuses one past a planar end
  face, by the test amendment "a crease's seat past an earlier blend"
  makes: straight in from the crease face's plane, the first surface met
  is the end face's own cylinder or sphere, the point met on the solid's
  boundary and the way to it in air. A groove beside the crease is met
  first over a groove, and is not excused.
- **Past the gap.** Carried on across the gap the curved face leaves (a
  groove), the band can meet the solid again farther along its line.
  That piece is no ball's of this crease, which ended in the gap. The
  tool is kept to a box about the piece holding the crease's own
  spandrel, a radius past it and clear of every other piece's box;
  where the boxes meet, the route declines as before.
- **Oracle.** The grooved block's bottom edge and groove crease rounded
  together (`blends.rs`) are held to the rolling ball: the groove band's
  opening of the floor less the drum and the wall band's spandrel
  carried on until the drum trims it, each section bounded by
  1-Lipschitz cells.

### Amendment: a corner one ball rounds, from its bands alone (#902)

A corner of three convex creases asked together, where one ball a radius
in from all three faces fits, closes on that ball (the corner sphere, the
default of section 13, question 2). On a solid with a curved face the
trimmed assembly can build the corner but not every band's far end (a
straight band running out across a drum at its other end), and the
request went to the current path.

- **Rule.** No new behaviour: the corner is the rolling ball's. Where the
  three faces are planes or walls holding the material inside them, the
  solid is locally convex there, and the round is what no ball within the
  solid reaches. A point of a crease's spandrel is reached by no ball
  within that crease's two faces, so by none within the solid: each
  crease's round alone, carried past its ends as it is alone, is part of
  the corner's round. What is left of the round is the corner's cell: the
  points whose nearest centre of a ball within the solid is the corner
  ball's centre (the cone of the three faces' normals there), farther than
  a radius from it.
- **Build** (`engine/ball_corner.rs`). The ball's centre is solved a radius
  in from the three faces (Newton's method on their signed distances). It
  is proved to touch each face within its trim and to stand at least a
  radius from every face of the solid. Each crease's round alone is a
  tool. So is the cell: the normal cone, held within the three faces'
  tangent planes at the ball's touches moved a radius out (the solid stands
  within them there), less the ball. They are applied as the one-boolean
  fallback's.
- **Oracle.** The shaved cube's corner (`blends.rs`) takes the union of its
  three bands alone, in closed form, plus the cell outside it, integrated
  by rays from the ball's centre over the normal cone.
- **Not this route:** a face of another kind or a wall holding the
  material outside it, a vertex with a fourth edge, radii that differ, a
  crease whose round alone the engine does not build, or a ball that does
  not fit.

### Amendment: convex bands crossing apart (#902)

Amendment "rounds standing apart" rounds each group of a request alone
and applies the tools where every tool stands clear of every other
group's. Two convex edges that share no vertex can still cut into each
other: the shaved cube's top edges along two sides, with the rim arc
between them left sharp, each run out across the drum at both ends, and
near the shaved corner the two bands' cuts overlap.

- **Rule.** Two groups' cuts may touch or overlap; the cuts take their
  union, where the request passes the overlap guards every request of
  two edges or more is held to (no two bands overlap across a face they
  share, as section 11's rounds that no ball fits). A request failing
  them leaves this route and is refused there by name: two edges of a
  face narrower than their two setbacks still are. A fill must still
  stand clear of every other group's tool, cut or fill.
- **Why it is the rounding asked.** Where two convex rounds meet across
  material both take off, each band removes what no ball within its two
  faces reaches, and a point that one band removes is reached by no ball
  of the solid there. This is amendment "M5e"'s convex-turn rule (each
  crease rounded as it is alone, the cuts taking their union), applied to
  groups whose cuts meet away from any shared vertex. A fill meeting
  another group's tool would change material that tool's round was
  proved on, so that case still declines.
- **Oracle.** The shaved cube's two top edges in one call (`blends.rs`)
  take twice one band alone, in closed form, less the square where the
  two overlap within the drum, integrated column by column.
- **Not this route:** a fill touching another group's tool, or two bands
  overlapping across a face they share.

### Amendment: a far edge leaving the seam (#902)

Amendment "M5d"'s run caps its fill on the arcs' planes past the seam's
ends, stepping out until an arc clears the end face's plane, and cuts the
fill back to that plane. Where the far edge the ball rests on leaves the
touched face's plane past the seam's end, faster than the ball can follow
(beyond some station no ball a radius over the touched face reaches the
edge's line), the balls run out before any arc clears the end face's
plane: the stress case `blend/boxes/48`, a short concave crease in a 20
degree wedge between two oblique boxes, where the first box's bottom
strip is narrower than the setback and its far edge runs off at 23
degrees to the crease.

- **Rule.** No new kind of end. At an end where the ball runs out before
  its arc clears the end face's plane, the run ends where its arc's plane
  crosses the seam's end, as a run rolling round a vertex does (amendment
  "a round rolling round a vertex"). The other end is capped and cut back
  as before. The current path ends the same band on the same plane at
  that end.
- **What is left.** Between that arc's plane and the end face's plane the
  wedge stays sharp, as the current path leaves it.
- **Proof** as M5d's, over the stations the fill spans.
- **Oracle.** `blend/boxes/48`'s two edges in one call (`blends.rs`): the
  fill is the ball's sections square to its path (bounded by the crease,
  the strip, the arc and the touched face, clipped to the end faces'
  half-spaces) moved along the path from where they leave the first end
  face's plane to the arc's plane through the seam's other end, by
  Simpson's rule; the other edge is a plain round between square faces.
- **Not this route:** a ball that runs out within the seam.

### Amendment: a bevel stops at a concave floor (#1145)

A downstream report: "a chamfer wider than the face or wall beside it
builds instead of being refused". Two shapes of it:

- **A bevel past a convex far side.** A 12 mm chamfer on the top edge of
  a block 10 high runs past the front's foot, a convex edge. That is the
  notch of amendment "M5b", deliberate: the bevel carries on and trims
  what it reaches (`1400` mm³ off a 20 × 20 × 10 block, the bevel's
  triangle less the part below the bottom).
- **A bevel past a concave far side.** A 2 mm chamfer on a step's back
  edge, the wall under it 1 high and meeting the block's top at a concave
  edge. M5b built it by its tools, cutting a groove 1 deep into the floor
  beyond the wall.

**Rule.** A convex bevel on a straight edge between planes whose leg on a
face crosses a concave edge of that face's outer ring, between the
edge's ends, is refused by name, naming its spec: its bevel would cut
below the floor or ledge the face ends on. The maintainer decided this in
#1145. It matches the rounds: a round whose ball would come down on the
floor beyond its wall does not seat and is refused (section 13, question
4: what does not fit fails, naming the edges). A notch trims the element
a blend overflows onto; a floor beyond a concave foot is material the
bevel would pass under, a change to a face the request does not name.

- A convex far side still rolls over (M5b).
- A bevel reaching under a feature standing within its face (a hole of
  the face, such as a block on the top) is non-local, as amendment "M5a"
  builds it, and not this.
- A leg that reaches the concave edge exactly builds.
- A rim's seat (a cylinder, a cone or a sphere under a planar cap) is
  read the same way (amendment "a rim's bevel stops at a concave
  floor"), and so is a crease on a marched seat (amendment "a marched
  seat's bevel stops at a concave floor").

**Oracle.** The step: at 1.5 the request of three edges is refused naming
the step's spec alone; at 1 the step's bevel alone takes `20 · 1²/2`.

### Amendment: tangent chains stop at a kink (#1145)

A downstream report: a fillet picked on one edge with tangent
propagation "runs past a 4 degree kink instead of stopping there". Two
questions: which edges a tangent chain takes, and what a round that
stops at a kink does at its end.

**Which edges a chain takes.** The default tangency tolerance, from
public documentation:

| Source | Default |
|---|---|
| Surface modeler A [P48] | The document's angle tolerance, 1 degree by default, decides when two curves or surfaces are tangent. |
| Parametric system B [P49] | The modeling angle tolerance, the largest angle allowed between tangents at matching points: 0.5 degree by default (a user forum gives the value; the vendor's API reference defines the setting). |
| Parametric system C [P50] | Tangency propagation runs "up to the first edge that is not continuous in tangency". No angle given. Its tangency check works down to 0.5 degree. |
| Parametric systems D, E and F [P51] | Tangent propagation or tangent chain is on by default. No angle given. |
| Parametric systems G and H | Not found. |

No source takes a turn of more than 1 degree as tangent. A 4 degree
kink stops tangent propagation in every documented system.

**Rule.**

- `tangent_chain` takes the picked edges and every edge reached through
  vertices where the way one edge leaves is the way the next arrives,
  within an angle. `TANGENT_CHAIN_ANGLE` is 1 degree, the widest
  documented default. Edges meeting back to back (a cusp) are not a
  chain.
- The angle chooses edges only. The engine runs bands on into each other
  with no corner only where the edges are tangent (`setbacks::TANGENT`,
  a millionth of a radian). It closes every other vertex of a chain as a
  corner, as it would if the edges had been picked one by one.

**A round that stops at a kink.** The reported part is a block whose
front turns 4 degrees into the material at x = 20, its first top edge
rounded alone. That is one rounded edge at a vertex of three convex
edges. Under amendment "corners where a round's rail misses" the band
runs on past the vertex while the turned front trims it, until the trim
reaches its contact on the top, `r cot(kink)` past the vertex. The
report's figure, `20 r²(1 − π/4)`, caps the band square at the vertex,
which leaves a cliff the ball does not leave.

The engine built that run-on only where it ended before the part's far
end. Where the far end came first, it treated the band as reaching a
face beyond its ends, and the cut-off route ended its tool square at the
vertex. So a 0.9 degree kink on a 20 mm face was capped and a 4 degree
one was not. Now the cut-off route's tool runs on too: at a vertex of
three convex edges it ends square where the band's last contact crosses
the end face's plane (the farthest setback corner past the vertex), and
the cut takes only what of that the solid holds.

**Oracle.** At `x` past the vertex the band takes what of its spandrel
lies at `x tan(kink)` or more from the front. With `A(s)` that area over
`y ≥ s` for the unit round, the round takes

`20 (1 − π/4) + cot(kink) ∫₀^m A(s) ds`, `m = min(1, (far − 20) tan(kink))`,

and `∫₀^1 A(s) ds = 5/6 − π/4`. For the reported part (4 degrees,
far = 40) that is `4.977542`, against the report's `4.292037`.

- [P48] Vendor help page, units and tolerances, and the vendor's
  tolerance FAQ.
- [P49] Vendor API reference, part modeling angle tolerance; the default
  from a user forum thread on modeling preferences.
- [P50] Vendor user guide, edge fillet propagation modes; vendor support
  record on the tangency check's lower limit.
- [P51] Vendor help pages on tangent propagation and tangent chain
  selection.

### Amendment: rolling corners on curved faces (#1145)

Amendment "rolling corners" refused every vertex where rounds meet that
touches a curved face or joins other than three edges. `Setback` builds
many of those: pocket floors and rims whose outlines have arcs, tangent
chains of lines and arcs, a pyramid's apex of four rounds, a round on an
edge already rounded (rounds meeting on its cylinder or torus), a kidney
pocket's floor of arcs alone.

**What the default builds there.** At each of those vertices the rounds
are all of one sense, and the default's corner is already the rolling
ball's (section 13, question 2):

- rounds meeting tangentially: the band runs on, with no corner;
- three rounds of one sense where one ball touches every face, on planes
  or a curved face: the ball's sphere (the corner tool of the curved
  route);
- an apex of four rounds or more: the envelope of the ball, its sphere
  where one ball touches every face, spheres and cylinders where none
  does;
- two rounds of one sense beside an edge left sharp: the bands cross.

None is a patch. The patch, and the setback, come only where convex and
concave rounds meet.

**Rule.**

- Where the rounds at every vertex of the request are of one sense, a
  rolling corner is the default's corner, built by `fillet_edges`.
- Where a vertex mixes convex and concave rounds, amendment "rolling
  corners" stands: between planes at a vertex of three edges, the torus
  round the lone band's end; elsewhere refused by name, naming the
  vertex.
- Mixed rounds on a curved face are refused by `Setback` as well ("a
  convex and a concave round meeting at one vertex close as one corner
  only where every face of the vertex is planar"), so every request
  `Setback` builds, `Rolling` builds too. `Rolling` now closes some of
  them (amendment "mixed rolling corners on curved faces").
- The order of the edges changes nothing, as for the default.

Making `Rolling` the default corner is a separate decision and is not
taken here.

**Oracle** (`tests/rolling_curved.rs`), each with `Rolling` in three
orders, sound and with no patch face:

- the obround pocket's floor and rim, the kidney pocket's floor, the
  square pocket's floor, a rounded rectangle's top loop and a rounded
  block's top loop: the section swept round the loop by Pappus, within
  1e-6;
- a pyramid's apex of four rounds, a notched corner of three convex
  rounds on a concave cylinder, uprights rounded into an earlier rounded
  top: as `Setback`;
- a D-shaped boss's foot with its uprights (mixed, on its cylinder):
  refused by name by both at the time; `Rolling` now builds it (amendment
  "mixed rolling corners on curved faces").

### Amendment: chamfers the engine refused on 0.19.0 (#1161)

A downstream report: chamfers that built before M6 were refused with
"the blend engine builds no chamfer of this kind here". Three causes.

**`Point` where no three bevels meet.** `ChamferCorner::Point` and the
default facet differ only at a vertex where three or more bevels meet
(section 13, question 1). The `Point` call took the mitre's own route for
every request, and that route lacked the default's routes for rims and
for solids holding seams. So a cone's or a boss's foot rim on a plate,
or an edge of a block holding seams, was refused with `Point` and built
by default.

- A `Point` request with no vertex of three bevels takes the default's
  route and builds exactly what the default builds.
- Where three or more do meet, the mitre now reads a solid holding seams
  with the seams dissolved, as the facet does.

**Faces a boolean imprinted inside a flat.** Where a fuse joins two parts
whose sides are coplanar (a bar fused between flanges, faces flush in
front), it can leave a piece of one face whose every side is a seam.
Dissolving the seams of a solid bounded by planes now matches such a
piece to the face that holds its neighbour across a seam, as the curved
reading already did. Before, the whole solid was read as one the engine
cannot chamfer. A concave edge where the bar meets a flange was then
refused, though the same part built as one prism chamfers.

**Concave bevels round a tangent loop.** In mainstream CAD, a chamfer
picked on a pocket's floor edge with tangent propagation takes the whole
loop of lines and arcs and builds one continuous bevel, as a fillet on
the same loop does. The engine's curved rings build convex bevels only,
and the one-boolean route for creases joined only pieces of one curve
between the same two faces.

- **Joints.** Asked edges meeting at a vertex are a joint when their
  curves leave it back to back within `setbacks::TANGENT`, they share one
  face (the floor), and their other faces stand tangent there. Every other
  edge at the vertex must be a seam or the smooth edge between those two
  faces.
- **Tools.** Each bevel's tool is its section swept along its edge (a
  prism along a line, the triangle turned about a wall's axis along an
  arc). It is capped on the plane square to the chain at each joint, cut
  where the edge is convex and fused where it is concave. Both sides of a
  joint share the plane and the section, so the tools meet face to face.
- **Unchanged.** A kink is no joint: amendment "tangent chains stop at a
  kink" stands, and such a request is handled as before. Amendment "a
  bevel stops at a concave floor" concerns convex bevels and is not
  touched.

**Oracle** (`tests/chamfer_point_and_pockets.rs`), sound, within 1e-6:

- the obround and kidney pocket floors at 1: the pocket plus `d²/2` along
  the lines and, by Pappus, round each arc, its centroid `d/3` from each
  face;
- the slanted bar's four concave edges at 1: `base + 20 cos(a)`, with
  `a = atan(13/36)`;
- the cone and boss rims, and edges of a block holding seams: `Point`
  equal to the default to the bit, both to Pappus or the prism's section;
- a corner of three bevels on the seamed block: `Point` the mitre, and the
  default `d³/12` less.

### Amendment: a cusp's S on one side while the top runs on (#900)

Amendments "rounds either side of a cusp" and "a cusp's S, one side or two"
refused one case by name: a plate's top edge and the boss's foot arc asked
on one side of the tangent point T, the top running on past T with its
edges there left sharp. They held that the S had nothing to close against
on T's other side. It does not need anything there. The one-sided S of an
arm running into a hub closes down to T on its own feet, and nothing past T
enters its hole.

- **Rule.** The lobe and its two feet close on one patch with its pole at
  T, as the one-sided S does where the top ends at T. The top's edges past
  T, with the side and with the wall, stay sharp and run on to T, where the
  S's feet meet them. The rim's arc past T is not asked and stays sharp.
  The asked arc's band ends at its far end in the meridian plane there,
  capped where it runs on into the rest of the circle.
- **The arc as asked.** A closed rim asked in part is rounded the whole
  turn where that builds (the curved routes' `whole_rims`). Here the whole
  turn is a different request: the foot's band past T would run beyond the
  plate's side next to the S. Where the whole turn does not build and the
  request as asked makes a cusp corner, the trimmed assembly takes the arc
  as asked.
- **A cap at a concave band's end.** Where a rim's arc runs on into the
  rest of its circle left sharp, the only faces at the vertex are the cap
  and the wall, both square to the cap, so their outward normals do not
  say which way the cap faces. The cap faces into the band where the edge
  is convex and out of it where the edge is concave. The cord between the
  band's contacts tells which: it lies behind the faces' outward normals
  at a convex edge and in front of them at a concave one.
- **Held** as the one-sided S is: the bands' fit target, G1 within 2e-4.
  On printCAD's plate it is G1 within 6.6e-5 at every radius from 0.25 to
  1.5, and its volume is within the patch's own distance of the S oracle
  (0.029 at r 1).
- **Still refused by name** (`engine::cusp_sides`): the same request
  wherever the trimmed assembly declines it.

### Amendment: what #877 left of asymmetric rounds (#907)

#907 listed five asymmetric-round cases still refused by name. The rules
below come from surveys this record already holds; none is new
behaviour.

- **Concave and mixed corners.** Amendment "engine-owned mixed asymmetric
  requests" holds that a vertex whose edges are all blended is blended. The
  engine's planar request already closes these corners on one patch,
  tangent to each round's band: a pocket's inside corner, and the mixed
  corner where a concave upright meets two convex top edges. They now have
  tests of their own.
- **Two of three blended, the third left sharp.** The same amendment's
  survey: two blended edges of one convexity meet by intersection, and the
  vertex is not rounded. That amendment refused this for an asymmetric
  band. It is now built:
  - each band's tool comes from its original seat and runs to its end
    faces, the circular round's ends;
  - the tools are applied together, so the solid loses (or gains) their
    union and the two band surfaces meet along their intersection.
  - Still refused by name: bands of opposite convexity, and a curved edge
    or a bevel at such a vertex.
- **Varying distances.** Amendment "varying rounds on curved seats" lays a
  law by arc length from the edge's start, held level past an open edge's
  ends and closed round a closed one. Amendment "engine-owned variable-law
  and conic tools" holds that every cross-section type takes variable
  values. So:
  - Round a rim and along any crease the march seats, the asymmetric band
    is marched under its two laws, each leg walked its station's distance.
  - On a straight convex edge an end not square to the edge takes the
    constant band's ends (M5f), the laws held level past it: run on to an
    end face leaning back, trimmed by one overhanging.
  - A concave edge stays capped square, and an end not square to it is
    refused.
  - No request takes a law over a chain of asymmetric rounds, so a chain
    stays out of reach.
- **Along a drum's ruling.** The section follows amendment "round
  distances on curved seats": each leg runs its distance along its face's
  section square to the crease (round the drum's arc), and the conic is
  tangent to both sections at the contacts. Along a ruling, with a flat or
  a second drum of a parallel axis, that section is the same at every
  station. The band is therefore the shoulder conic extruded along the
  edge, exact, capped square at end faces square to the edge. A ruling
  ending at a face not square to it is left to the march.
- **The largest pair on other seats.** On a seat that is neither a
  straight edge between planes nor a rim, each leg is walked as the band
  walks it, at stations over the edge. The scale is stepped up until a
  contact first leaves its face, then bisected, and the answer stands back
  by the reach the faces are read within. Building the band is the only
  check on such a seat, so the answer is the elliptical round built on a
  copy of the model; a refusal there is `NotDone`.

**Oracles** (`tests/asymmetric.rs`):

- the ruling, convex and concave: the drum less (or plus) the drum's
  height times the shoulder conic's corner, less (or plus) the sliver
  between the tangents and the arc, `r² (tan β − β) / 2`; to 1e-8;
- two of three: the box less each band's corner along its edge, plus once
  what both take. That is the column of the conic's depth times the
  round's width, integrated over the left face, or the shallower depth
  over the top. Sampled independently, to 1e-7 and 1e-5;
- the pocket's two fills: the same, with the shallower fill;
- the concave and mixed corners: the bands' own sections past the corner's
  cube, the remainder within the cube, and the patch tangent within
  2e-4;
- varying laws: the integral of the conic corner over the edge, plus a
  quarter of the end section's first moment where the end leans back by a
  quarter, or the part of each section under the overhang, to 1e-6. Round
  a rim, Pappus station by station; along the spline crease, the corner
  less the crease's curvature times its first moment. Both within 1e-4
  of what the band takes, the march's fit;
- the largest pair: 10 along the half drum's ruling (the flat is 10
  wide) and 2/0.7 on the slab's spline crease (the wall is 2 deep), each
  to 1e-4. The answer builds valid, and 1.01 times it is refused.

### Amendment: the rest of asymmetric rounds (#907)

The amendment before this left four asymmetric-round cases refused by
name. The public user documentation of mainstream CAD was surveyed on
2026-11-20 for the two of them this record had no rule for:

| Question | What mainstream CAD documents |
|---|---|
| A varying size along a chain of tangent edges | Tangent propagation carries a variable fillet across the tangent edges; sizes are set at vertices and control points along the edges, and the value at a shared vertex is one value for both edges. Asymmetric (two-distance) fillets take variable sizes too, with control points at vertices where filleted edges meet. |
| A concave fillet whose edge ends at a face not square to it | Not documented. Forum answers describe the fillet stopping on that face, and split faces or lofts as workarounds where it does not. No source asks for a square cap past the face. |

Rules:

- **A chain of varying distances.** `fillet_edges_distances_variable`
  takes a tangent chain, as `fillet_edges_variable` takes one for a
  radius. Both laws are laid by arc length along the walk from the first
  edge's start, closed round a closed chain, so a junction reads one
  value for both edges. The named face lies beside every edge.
  - Each straight edge between planes takes its exact band under its
    stretch of the laws. Any other edge takes the marched band, which now
    stops flush in its end section where a chain mate continues the
    crease. The bands are joined across each junction into one tool.
  - Measured limit: where both laws vary at a junction between a straight
    edge and an arc, the two bands share the end section but their
    tangent planes across it differ by up to a few hundredths of a
    radian. The sections stand in planes square to the crease, and the
    marched band's along-derivative at its trimmed end is held only to
    the fit. The circular law's chain shows the same (about 0.016 rad on
    the same block). Volume and contacts are within the march's fit.
- **A concave edge ending at a face not square to it.** The fill ends in
  the plane of each face at its end. Where one is not square to the
  crease, the fill runs on with its laws held level and is trimmed by
  each such plane, kept on the side of the material that face bounds. A
  fill capped square stood past a face leaning in (the constant fill
  failed there with folded faces) and left a notch beside one leaning
  out. The constant fill and the varying one now end alike, so level laws
  are still the constant round to the bit. An end meeting a face that is
  not a plane stays refused by name.
- **Two bands at a vertex whose third edge is left sharp.**
  - One cuts and the other fills: question 6's rule. The bands end
    through their faces' setback corners, and the signed patch closes
    the corner. That patch reads only the strips' end sections and the
    setback corners, so it now takes conic sections as it takes circular
    ones. A request with a bevel, a concave sharp edge, or setback corners
    that do not stand either side of the vertex stays refused by name.
  - A bevel and a round of one sense: they meet by intersection, as two
    rounds do. This is what question 3's "bevels first" gives here: the
    round's edge, shortened by the bevel, runs on to the bevel's plane,
    which leans back from it. The bevel's tool is its prism between its
    ends, so each of its ends must meet only faces square to it or
    running along it, or it is refused by name.
  - A curved edge at such a vertex stays refused by name.
- **The largest pair on curved seats.** The seat reading also stops where
  a station's faces, read in its section plane, leave their contacts
  along tangents that no longer meet in front of the crease: the shoulder
  the conic stands on. On a convex ruling that is a leg a quarter round
  its drum (the shoulder at infinity). On a concave one it is the drum's
  tangent crossing the other face's line past that face's contact.
  - The concave bound is where the build itself refuses, so the confirm
    step brackets it at once.
  - On a convex ruling the section gives out before the shoulder does.
    As the drum's leg turns far past the flat's, the conic standing on
    the shoulder bulges toward the drum, and from some scale it crosses
    the drum between the contacts: no band of that section is tangent to
    the drum there without leaving it. With legs 0.1 and 1 on the half
    drum that is at 1.69143, well short of the quarter turn (5π/2 =
    7.854). The ruled section is refused by name there ("crosses the
    drum"), read exactly: the arc against the drum's circle is a quartic
    in the arc's parameter, the contact's double root divided out, and
    the rest a quadratic whose sign over the arc is read in closed form.
    The seat reading on a ruled seat applies the same test, so seat and
    build agree: the confirmed answer is 1.69140, refused at 1.69150.
  - Past that scale the build used to go on: the boolean took the band
    bulging out of the drum and returned a solid with material outside
    the original, valid by `check`, up to about 1.9, and refused only
    from there. So the old answer of 1.9085 was already past the
    geometry's limit.
  - Just below it, from about 1.675, the boolean refused bands that clear
    the drum by under 1e-4. The band's tool shared the drum's end planes,
    and in each the band's end conic, hugging the drum's rim circle, was
    read as running along it. Where the solid lies wholly behind an end's
    plane (read on its guaranteed bound) a cut's prism now runs on past
    that end, so the end plane cuts the band in its own exact section.
    A fill does not run on: past the solid it would add material.

**Oracles** (`tests/asymmetric.rs`):

- the chain round a block's rounded upright: each station's conic corner
  integrated along the straight edges, and round the arc times
  `(R - c) / R`; within 1e-4 of what the band takes, the march's fit.
  Valid, each band tangent to its faces within 1e-3;
- the concave fill on the slab on the leaning block, both leans: constant
  distances `A (L - k c)` to 1e-8; varying laws, the corners' integral plus
  the end's run-on moment, or each section under the overhang, to 1e-6.
  The same solid whichever face is named, and tangent to both faces
  within 1e-7;
- the L-bracket's signed corner with an asymmetric round on either band
  or both: the bracket less the convex corner along 12, plus the fill
  along 20, the remainder within the cube of the widest leg. Every point
  near the corner lies within the bracket or the fill (2e-4), and the
  patch is tangent to the bands but along the sharp edge (2e-4);
- a bevel beside an asymmetric round (box) or fill (pocket): each section
  along its edge, less (or plus) once the shallower depth over the face
  under both; to 1e-5;
- the largest pair on the concave ruling: the root of
  `r (1 / cos(b s / r) - 1) = a s`, to 1e-4 relative. The answer builds
  valid, and 1.01 times it is refused.

### Amendment: overflow handling

A round often runs into geometry that is not one of its own two faces: an
existing round, a sharp edge, another feature. Amendment "M5c, non-local
rounds keep their surfaces" settled one behaviour for all of it: the round
keeps its own surface and the feature it reaches is trimmed or extended.
§13, question 8 records the other behaviour practice offers there, the ball
resting on what it meets. This amendment makes the choice a request
option, as mainstream CAD does.

**The survey.** The public user documentation of mainstream CAD, read in
general terms:

- Every system with an overflow option offers more than one behaviour and
  resolves them in a fixed order, with every behaviour on by default.
- The first is a **smooth** overflow: where a round reaches a face joined
  tangentially to its own (most often an earlier round), the ball rolls on
  across the tangent edge onto that face, and the new round's surface
  bends to meet it with no step and no crease.
- The second keeps a **sharp edge** the ball meets: the round there is
  swept by a ball touching one of its original faces and resting on the
  edge, which is kept on its curve. Several systems also apply the rule at
  a band's end alone, as a separate switch.
- The last keeps the round's own surface and **trims or extends** the
  faces it reaches (the notch, or "keep surface").
- A system that offers no option builds the last, or refuses.
- A tight concave zone, narrower than the ball (a narrow passage), is a
  separate setting in the same option set where it exists: a choice of how
  general (and how slow) the solver may be.

**Rule.** A fillet request carries `Overflow { blend_into_round,
hold_edge, hold_edge_at_end, trim_through }`, every field on by default.
Where a round reaches geometry other than its own two faces, the
behaviours are tried in that order and the first that builds wins:

1. **Blend into round**: the round reaches an existing round of the body.
   In the transition zone the ball rests on the existing round instead of
   on the face it left, so the new round meets it G1. This is §13
   question 8's smooth overflow, the seat moving across the tangent edge.
2. **Hold edge**: the ball would roll over a sharp edge of the body. The
   edge is kept, and the round there is swept by a ball touching one
   original face and resting on the edge (§13 question 8's edge seat).
   **Hold edge at the end** applies the rule where only a band's end runs
   over such an edge, and matters only with hold edge on.
3. **Trim through**: amendment "M5c"'s rule, today's construction.

A behaviour switched off passes the request to the next. Where every
behaviour that is on fails to build a round that reaches other geometry,
the request is refused by name; trim through stays the construction where
it is on, so the default request builds what it built before.

**What is built in this step.**

- **Blend into round**, for a boss's foot: the closed rim between a
  planar floor and a cylindrical wall standing square on it, whose round's
  contact circle crosses one straight edge of the floor along which a
  convex cylindrical round meets the floor tangentially.
  - Over the floor the round is the torus its ball sweeps.
  - Past the tangent line the ball rests on the wall and on the existing
    round. Its centre runs where their offsets meet: the cylinder about
    the wall's axis and the cylinder about the round's axis. The sections
    are exact at every station, and the piece is fitted through them
    within a confusion.
  - The two pieces meet where the ball touches the tangent line. The
    spine's tangent is level there, so the section lies in the plane
    through the wall's axis on both sides: one quarter circle, G1.
  - The fill between the round, the wall, the floor and the existing round
    is one closed tool, fused into the solid by amendment "M5a"'s boolean.
    It is built only where nothing of the solid but those three faces
    stands within the box the fill can take up, which proves the fill
    touches nothing else. The fuse must add the tool's whole volume, and
    the result must check valid.
  - Any other round reaching a round (two rounds at once, one not tangent
    to the floor, a concave one, a wall not square to the floor, an open
    arc of a foot) is not blended in this step and goes on to the next
    behaviour. The part "Blend into round, more feet" below builds the
    first three of those.
- **Hold edge** and **hold edge at the end** were not built in this
  step; the amendment's "Hold edge" part below builds them.
- **Trim through** is the engine's construction unchanged. With it off,
  a result that changes any face other than the edges' own and those at
  their ends is refused by name. A face counts as unchanged where the
  result holds a face of its area and bounds on its surface. A foot that
  blend into round would take is refused by name when both are off.
- **Narrow passages** were not part of this step. The part "Narrow
  passages" below builds them under the same option set.

**Oracle** (`tests/overflow.rs`). The case is a block 40 × 30 × 15 whose
long top edge is rounded r 5, with a boss r 6, 15 tall, standing at
`y = 12`. Its foot round, r 3, reaches 2 past the edge round's tangent
line.

- The result is valid and closed, and every edge of the fitted piece
  meets its neighbours within 1e-5 rad of tangent.
- The reference volume is the morphological closing, computed without
  the construction. A point above the floor is in the fill where its
  distance from the centres clear of the solid exceeds the radius. That
  distance comes from projections onto the two offset boundaries and a
  search along the spine. The columns are integrated by Gauss–Legendre
  quadrature, split where the seat changes. The reference reproduces
  Pappus for a boss clear of the edge round to 1e-6. Its own error bound
  is the change from 48 to 96 points per range, 6e-6. The fill the round
  adds matches the reference within that bound plus 1e-6.
- Blend into round off: trim through's refusal, as before. All off: a
  refusal naming blend into round. Trim through off on a rim trimming
  through a pocket's walls: a refusal by name. A local round with every
  behaviour off builds the same closed-form solid as the default.

**Hold edge.** The public documentation surveyed above defines the rule
by the ball: it touches one original face and rests on the kept edge. It
gives the end case as the one where only the ending part of the round
overflows a sharp edge, as a sub-option of hold edge, and says that where
hold edge and trim through could both apply, hold edge is applied. It
describes no geometry past that, so the following is settled here.

- **Where the ball rests.** Along the stretch (the *run*) where the band's
  contact on face `a` would lie past a sharp convex edge of `a`, the ball
  touches the other face `b` and rests on the edge: its centre stays on
  `b`'s offset at the radius from the edge's curve, found from the band's
  centre toward `a`'s plane. The round there is the canal surface those
  balls sweep. Each ball's arc runs from its contact on `b` to its foot on
  the edge; both lie on its characteristic circle, in the plane through
  its centre square to the spine. The edge lies on the canal, so it is
  kept on its curve as the edge between the round and the face past it.
- **The seat change.** Where the contact crosses the edge, the ball
  touches `a`, `b` and the edge at once: the offset of `a` and the pipe
  about the edge are tangent there, the spine's tangent runs along the
  edge, and the section is the band's own. The canal meets the band G1
  on that section.
- **At the end.** A run that reaches one end of the band and not the other
  is *at the end*, and is held only with both hold edge and hold edge at
  the end on. A run inside the band is held with hold edge on. A run over
  the whole band is not held: the round would take face `a` whole, so it
  is larger than its seat allows, and the largest-radius readings stand. With the option it needs off, the request goes on to
  trim through, and with trim through off too it is refused by name,
  naming the option.
- **Past an end face.** The centres run on as they are past the band's
  end, the held edge's curve carried on past its vertex there, and the end
  face's plane cuts the round (amendment "M5f"), as it cuts a band. No
  vertex seat is built: a run that would need the ball to rest on a
  vertex of the held edge (one ending inside the band) is not held.

**What is built for hold edge in this step.**

- A round of one straight edge between two planes meeting convexly at any
  angle, ending flush on an end face square to the edge at each end. The
  held edge is a line or a circle (a whole one or an arc) in the plane of
  `a`, convex, its face falling away from `a`. There is one run, with a
  seat change at one end at least, and the balls' feet over it stay on
  the held edge but within the radius of an end face the run reaches.
- The canal is fitted through its sections by the planes square to the
  edge, each point found exactly on a ball's arc, within a confusion. Its
  ends lie in an end face's plane or on a band's section.
- The region between the edge and the round (the bands, the canal, the
  planes of `a`, `b` and the end faces) is one closed tool cut from the
  solid by amendment "M5a"'s boolean. It is built only where nothing of
  the solid but those faces and the face past the held edge stands within
  the box the tool can take up. The cut must take the tool's whole volume,
  the result must check valid, and its edges on the held edge's curve must
  still hold every ball's foot over the run. Otherwise the request goes
  on to trim through.
- Anything else (a chain of edges, a curved face, a held edge that is
  neither a line nor a circle, a run needing a vertex seat, two runs, an
  oblique end face) is not held in this step and goes on to trim through.
  The part "Hold edge, more seats" below builds vertex seats, several
  runs, chains along one line and oblique end faces.

**Oracle** (`tests/hold_edge.rs`). The reference is the definition,
sampled without the construction: the ball's centres in closed form (at
the radius under the top while the contact lies on it, else at the radius
from the held edge's curve), and a point of the corner cut away where it
lies farther than the radius from every centre, on the corner's side of
its nearest one. The cut is integrated by Gauss-Legendre quadrature in
`x` and down the front in `s`, `z = -s^2`, split at the seat changes. It
reproduces the band's closed form to 1e-9. Its bound is the change from
24 to 48 points per range.

- Hold edge: a block 30 × 20 × 15, its front top edge rounded r 4, a
  pocket r 5 sunk 8 into its top about `(15, 7)`. The rim comes within 2
  of the edge, so between `x = 11` and `x = 19` the ball rests on it.
- Hold edge at the end: a block 20 long whose top narrows from 6 to 2
  deep, a cliff 6 high falling from its back edge, the front top edge
  rounded r 4. Past `x = 10` the ball rests on the cliff edge to the end.
- Each result is valid and closed, the held edge is kept on its curve over
  the run, the volume cut away matches the reference within its bound
  plus 1e-6, and the canal meets the band and the front within 1e-5 rad of
  tangent.
- Hold edge off: trim through's notch, the rim cut away. Hold edge at the
  end off: the same at the end, with no change to a run inside the band.
  Either off with trim through off: a refusal naming it. Trim through off
  alone: hold edge builds.

**Chamfers.** A chamfer takes the same option set, limited to the
behaviours that apply to a flat bevel. Mainstream CAD offers the same
choice for chamfers where it offers one for rounds, without the edge seat.

1. **Blend into round**: the bevel reaches an existing round of the body.
   In the transition zone the bevel's surface is adjusted to meet the
   round with no step: its leg on the face it left runs on over the round.
2. **Trim through**: amendment "M5b"'s rule, today's construction. The
   bevel keeps its own surface, and the faces it reaches are trimmed or
   extended to meet it.

**Hold edge** and **hold edge at the end** do not apply: a bevel has no
ball to rest on an edge. They are ignored for chamfers. Both behaviours
are on by default and tried in that order. With both off, a chamfer that
reaches other geometry is refused by name; it never yields a wrong solid.

**What is built for chamfers in this step.**

- **Blend into round**, for the same boss's foot as the round's, every
  edge of the rim asking the same legs (any of the three spellings),
  whose leg on the floor crosses the tangent line of a straight convex
  round meeting the floor tangentially, and stays on that round.
  - Over the floor the bevel is the cone between its legs.
  - Past the tangent line the legs are read as on any curved seat: the leg
    on the round is where the wall's offset (the cylinder `out` from the
    wall) meets the round, straight below the circle it draws on the
    floor, and the leg on the wall stays where the floor's offset meets
    the wall. Each section is the straight ruling between the legs in the
    plane through the wall's axis. The piece is fitted through its exact
    sections within a confusion.
  - The two pieces share the ruling where the leg crosses the tangent
    line. The round leaves the floor tangentially, so the leg's tangent
    is the same on both sides there, and the bevel meets the cone G1.
  - The fill is one closed tool fused into the solid, under the round's
    proofs: nothing else of the solid within its box, the fuse adding its
    whole volume, the result valid.
- **Trim through** is the engine's construction unchanged. Before this
  step the foot's bevel was built by it: the cone runs on, and the floor's
  plane is carried on under it as a ledge over the round.
- With trim through off, a result changing any face other than the edges'
  own and those at their ends is refused by name, as for rounds. A foot
  that blend into round would take is refused by name when both are off.
- `chamfer_edges`, `chamfer_edges_with` and `chamfer_edges_cornered` are
  the default set, so the foot bevel above now blends into the round
  where it used to leave the ledge.

**Oracle** (`tests/chamfer_overflow.rs`). The block and boss are the
round's, the foot bevelled 3 along the top and 3, then 2, up the wall: its
leg on the top reaches 2 past the tangent line.

- The result is valid and closed. The bevel meets the round along edges
  lying on it, with no ledge, and the two edges where the fitted piece
  meets the cone are tangent within 1e-5 rad.
- The reference volume integrates the bevel's definition without the
  construction: at each angle the straight line from the wall's leg to the
  floor below the leg on the floor, less the floor, by Gauss-Legendre
  quadrature split where the leg leaves the top and where a column crosses
  the tangent line. It reproduces Pappus for a boss clear of the round to
  1e-9. Its bound is the change from 48 to 96 points per range. The fill
  matches it within that bound plus 1e-6.
- Blend into round off: the solid built before this step, face for face
  and of the cone's whole fill. Both off: a refusal naming blend into
  round. A notch with trim through off: a refusal by name, and with it on
  the notch as before.

**Blend into round, more feet.** The public documentation surveyed above
describes the smooth overflow by its result alone: the round runs on over
the tangent face it reaches with no step. It says nothing of a round that
reaches several such faces, of a wall leaning from square, or of a concave
face. The papers define a constant-radius round by the rolling ball, and
its fill as the closing of the solid by that ball [3]; that is the rule
settled here, and it needs no choice beyond the seat.

- **Several rounds.** Where the contact circle crosses more than one
  tangent line, the ball rests on each round over its own arc and on the
  floor between, so the round is a torus piece between each two fitted
  pieces, meeting each G1 at its seat change as in the first case. It is
  built only where the arcs stand apart, so that no ball rests on two
  rounds, or on the corner where they meet, at once: at a corner the top's
  corner stands farther from the axis than the contact circle. Arcs that
  meet or overlap are not blended.
- **A leaning wall.** The wall is a cone about an axis square to the floor,
  narrowing away from it. The ball touches the cone and the floor or the
  round: its centre runs on the cone's offset (a cone about the same axis)
  where it meets the floor's. Over the floor that is a circle, so the
  round is a torus of that circle; past the tangent line the centre's
  height and radius solve together. The seat change keeps its section in
  the plane through the axis, since the round leaves the floor
  tangentially and the cone's offset depends on height alone, so the
  pieces meet G1 there. A cone widening away from the floor (an undercut),
  or leaning more than 45 degrees, is not blended.
- **A concave round.** The round curves up from the floor: its axis
  stands above the floor at its radius. The ball climbs it, its centre on
  the cylinder of the round's radius less the ball's about that axis. It
  is built only where the ball's radius is less than the round's and the
  centre stays inside that cylinder short of its edge, so every ball lies
  inside the round's cylinder, clear of the faces beyond it.
- **Proof.** The fill's room is the union of boxes between each two of
  some hundred stations per arc, each box holding both stations' sections,
  their feet and the floor under them, padded by the chord's sag. Nothing
  of the solid but the wall, the floor and the rounds crossed may stand in
  it, and no edge of those faces but the rim and the tangent lines. Every
  other face whose box meets a ball's is measured from each centre: it
  must stand farther than the radius plus the stations' spacing. The fuse
  must add the tool's whole volume and the result must check valid, as
  before.
- **Not blended** (these go on to trim through): an open arc of a foot
  (the rim must turn once round the axis), a round not tangent to the
  floor, a curved floor, and a bevel reaching anything but one convex
  round from a cylindrical wall.

**Oracle** (`tests/blend_into_round.rs`). The reference is the closing,
sampled without the construction, as for the first case, generalised: the
free centres are bounded by the wall's offset (a line in each plane
through the axis) and the floor's offset (the plane at the radius, each
round's offset cylinder past its tangent line). At each angle about the
axis the columns are integrated from where the fill starts on the wall to
where it ends on the floor, both found by bisection, split at the wall's
foot circle on a cone, and the angles are split at every seat change.
Gauss-Legendre throughout; the bound is the change from 32 to 64 points
per range. On a conical boss clear of any round it reproduces Pappus
(the section is the quadrilateral of the foot, the touches and the centre,
less the ball's sector) to 1e-6.

| case | fill built | reference | error | bound | worst G1 |
|---|---|---|---|---|---|
| boss r 6 at (12, 12) by a corner of two rounds r 5 | 75.108708877 | 75.108707894 | 9.8e-7 | 1.8e-5 | 9.7e-7 rad |
| cone r 6 to 4 over 15 at (20, 12), edge round r 5 | 55.909500322 | 55.909500788 | 4.7e-7 | 3.0e-6 | 5.6e-7 rad |
| boss r 6 at (20, 21), concave round r 6 meeting the floor 7 off | 92.217220722 | 92.217224303 | 3.6e-6 | 7.3e-6 | 2.1e-6 rad |

Each result is valid and closed, with one fitted piece per round reached.
With blend into round off each foot gets trim through's refusal, as
before, and with every behaviour off a refusal naming blend into round.
An open arc of a foot reaching the edge round gets what trim through
gives.

**Narrow passages.** A concave zone is a narrow passage where it is
tighter than the radius, so that the ball cannot roll through touching its
seats on both sides: two walls of one floor closer than twice the radius,
or a corner face whose curvature is tighter than the radius. The public
documentation surveyed above says, in general terms:

- The setting is a choice of how general the solver may be. The cheapest
  assumes the tight zone is made only by constant-radius rounds and
  refuses anything else; the next assumes the zone runs through one whole
  face of a wall; the most general lets the zone cover part of a wall
  face, the round being split to suit. Turning the general ones off keeps
  rebuilds fast. Nothing there describes the resulting geometry.
- Where two rounds of one operation intersect each other, the
  intersection is computed and the rounds are trimmed to it, so that the
  result is a valid solid [14].

The rolling ball gives no answer here: no ball touches both seats. The
rule settled here is the second statement's, carried across the zone:

- **Every round keeps its own constant section.** Each wall's band runs on
  past its edge's ends through the tight zone, to the plane of the next
  wall (a square corner) or to the plane square to the edge through the far
  side of the face closing the zone (a half cylinder between facing
  walls). The fill is the union of the bands' fills, cut to the solid's
  air: where two rounds cross they meet along their intersection, a G0
  crease, and the faces of the tight zone (a floor narrower than twice the
  radius, a corner face) are buried where the fill covers them. Each round
  still meets its own wall G1.
- At a square corner whose corner face is a cylinder tangent to both
  walls and tighter than the radius, the ball at the bands' meeting stands
  clear of that face, so the union is the rolling ball's own fill there,
  the two bands mitred as at a sharp corner. Across a narrow floor no ball
  is clear of the far wall, and the union is the trimmed intersection
  above.
- **The strategies** are the three switches `narrow_blend_made`,
  `narrow_whole_face` and `narrow_partial_face`, on by default and tried
  cheapest first. They say what the tight zone may be made of, and each
  takes what the ones before it take:
  - **blend-made only**: constant-radius rounds only. These are the
    request's own rounds crossing over a narrow floor, and corner faces
    that are cylinders tangent to both their walls (a pocket's corner
    rounds, a stadium's ends);
  - **whole face**: also one whole corner face that no constant-radius
    round made, such as a cylinder meeting its walls at kinks;
  - **partial face**: also a zone covering part of a face, the band split
    where the zone starts and ends.

  A passage that the strategies switched on do not take is refused by
  name, naming them. These switches stand apart from blend into round,
  hold edge and trim through, whose order and option-off behaviour they
  leave unchanged.

**What is built for narrow passages in this step.**

- A planar floor's rim of concave edges, with planar walls square to the
  floor. Each pair of neighbouring walls meets in one of these ways:
  - at a square corner, sharp or through a cylindrical corner face square
    to the floor and tighter than the radius, which stands within the
    radius of both walls;
  - across a half cylinder tangent to two facing walls.

  Open ends must be flush on end faces square to the edge. There must be a
  tight corner face, or two facing walls no farther apart than twice the
  radius but farther apart than the radius. Exactly twice the radius
  apart, the two contacts meet on one line of the floor and the two rounds
  make one half circle, the trough's section.
- Each band's fill is its spandrel swept along the band, one prism, and
  the prisms are fused into the solid one boolean each. They are built
  only where nothing of the solid but the floor, the walls, the corner
  faces and the end faces stands within the prisms' boxes, and none of
  those faces' edges but the request's and those between two of them.
  The fuses must add material, no more than the prisms hold, and the
  result must check valid. Otherwise the request goes on to the engine.
- **Not built:**
  - a zone covering part of a face. Nothing needs partial face alone
    yet; such a request goes on to the engine, which refuses it by name;
  - walls not square to the floor or to each other;
  - a curved floor;
  - a corner face that is not a cylinder;
  - facing walls no farther apart than the radius (a band would cross the
    wall it faces);
  - walls lower than the radius.

  These go on to the engine, which refuses them by name.

**Oracle** (`tests/narrow_passage.rs`). The reference is the union's
definition, sampled without the construction. Over a point of the floor
at distance `s` from the nearest wall plane, the fill stands
`H(s) = r - sqrt(r^2 - (r - s)^2)` tall. It is integrated over the air's
plan in one variable: the strip `s` in from the walls of a rectangle is
`2 (w + l) - 8 s` long, a stadium's chord at `s` is known in closed form,
and each corner face's material is taken off by its own integral. That
integral is checked against a direct integral over both halves of the
corner. Gauss-Legendre throughout, with `s = t^2` where a square root
starts. The bound is the change from 24 to 48 points, below 1e-15 in
every case. A through slot has a closed form.

| case (r 4, floor 7 deep) | fill built | reference | error | worst G1 (round to wall) |
|---|---|---|---|---|
| through slot 6 wide, 20 long | 135.662523995 | 135.662523995 | 1.1e-12 | 0 rad |
| stadium pocket 20 by 6, ends r 3 | 119.916295109 | 119.916295109 | 1.5e-12 | 0 rad |
| pocket 20 by 14, corner rounds r 2 | 198.927647132 | 198.927647132 | 6.1e-12 | 0 rad |
| pocket 20 by 14, corners r 3 at kinks 2 in | 194.342573041 | 194.342573041 | 2.9e-12 | 0 rad |
| through slot 8 wide (twice the radius), 20 long | 137.345175426 | 137.345175426 | below 1e-6 | 0 rad |

- Each result is valid and closed. The narrow floor is taken whole.
- With blend-made off, a blend-made passage gets the same solid from
  whole face, and with whole face off too, from partial face. With all
  three off it is refused by name. The kinked corners with whole face and
  partial face off are refused by name.
- A slot twice the radius wide or wider is no narrow passage: it builds as
  before with every strategy off.
- A slot no wider than the radius, and one shallower than the radius, go
  on to the engine's refusal by name.
- A slot exactly twice the radius wide answers within a time bound: the
  two rounds give the half-circle trough, and one round alone gives its
  own band. The ball of that one round stands tangent to the far wall all
  along, so the far-edge proof that every face clears the balls has no
  margin anywhere. That proof halves each interval between its stations,
  and it now stops after 2^17 distance measures per interval, the face
  then not shown clear, instead of walking the whole tree of 24 halvings.
  The most a passing proof uses is about 46,000.

**Hold edge, more seats.** The public documentation surveyed above
gives the edge seat by the ball (it touches one original face and is
supported by the kept edge) and adds that the ball may roll onto
adjacent edges; it gives no geometry for where one kept edge stops and
the next begins, for several stretches on one round, or for a round's
end on an end face that leans. §13, question 8 already rests the ball on
whatever of the boundary it meets, a vertex among the seats, its offset
the sphere about it [3]. The papers on rolling-ball blends define the
round as the envelope of a ball kept in contact with the boundary [2].
What follows is settled from those.

- **Where the ball rests.** The ball touches `b` and stands as near
  `a`'s plane as the region of `a` (its faces in that plane, taken as
  one) lets it: its centre on `b`'s offset at the radius from the nearest
  point of that region. That point is the ball's contact on `a` where the
  contact lies on `a` (the band), a point inside a sharp edge of `a` (the
  edge seat of the part above), or a vertex where two sharp edges of `a`
  meet (a vertex seat). This reads the single run of the part above as
  before.
- **A vertex seat.** Where the foot on a held edge reaches the vertex
  where that edge stops and the next goes on, the ball rests on the
  vertex: its centre stays the radius from it while it touches `b`. Every
  ball of that piece passes through the vertex, so the piece closes to it
  there. The nearest point moves from an edge to a vertex only where the
  vertex is a convex corner of `a`'s region; the centre's path is then
  smooth (C1) across both changes, the two offsets are tangent, and the
  pieces share the ball's arc there, which lies in its plane square to the
  spine. The pieces meet G1 on it.
- **A ball on two edges at once.** Where `a`'s region turns the other
  way at a vertex (a pocket's corner, a notch), the nearest point jumps
  from one edge to the other without passing the vertex, the ball resting
  on both at once. The centre's path turns a corner there, and the balls
  next to the corner reach into the pocket past both edges, so the round
  would pass through the pocket's walls and the edges could not both be
  kept. That is not held; the request goes on to trim through, as before.
- **Several stretches.** Each stretch where the ball leaves the band is
  held on its own terms, and the band runs between them. With hold edge at
  the end off, a round any of whose stretches reaches an end of the band
  goes on to trim through as a whole.
- **A chain.** Edges end to end along one line, each between a face of
  the same two planes, are one round over the line from the chain's first
  vertex to its last; the seams between faces of one plane are not edges
  of `a`'s region. Inside the chain no other face may meet it.
- **An oblique end face.** An end face whose normal stands within 60
  degrees of the edge cuts the round in its own plane, as amendment
  "M5f" cuts a band: the band in the ellipse its plane cuts from the
  band's cylinder, a fitted piece through sections by planes that turn
  from the piece's other end to the end face's plane. Where no ball's
  foot reaches the end face's plane (a held edge leaving the end too
  steeply, its feet trailing the balls), the piece has no end there and
  the request goes on to trim through.
- **The boss's underside (amendment "M5c").** There the contact on the
  wall runs into the boss, not past a sharp edge into air: the edge where
  the wall meets the boss's underside is concave, and the vertex where it
  meets the wall's cut is not a corner of the wall's region the ball could
  rest on from outside. It does not reduce to a vertex seat, and stays
  with trim through.

**What is built for hold edge in this step.** A round of a chain of
straight edges along one line between two planes meeting convexly, ending
on end faces whose normals stand within 60 degrees of the edge. The edges
of `a`'s region within reach are lines or circles in its plane. The
ball's seat is read at a thousand stations and each change found by
bisection; the pieces between are kept only where every change is the
band onto an edge, or an edge onto a vertex where it stops, and no piece
is shorter than a ten-thousandth of the radius. A piece on an edge or a
vertex is fitted through its balls' arcs, or, reaching an end face,
through its sections by planes turning to that face's plane. A piece on a
vertex closes to the vertex along one side of its patch; its pieces' ends
there are that patch's own boundary iso-curves, whose pcurves are the
chart's straight sides, read exactly. The tool, its proofs and the
options are those of the part above, every face past a held edge or at a
held vertex counted with the faces past the held edge.

**Oracle** (`tests/hold_edge_seats.rs`). The reference is the definition,
sampled without the construction: the ball's centre stands the radius
behind the front, at depth `sqrt(r^2 - h^2)` under the top, `h` the
distance in the top's plane from its contact to the top's region, worked
from the region's polygon or pockets directly. The cut is integrated as
in the part above, split at the seat changes. On the leaning ends: the
band's part near each end is its spandrel over the length between the end
planes, a closed form from the spandrel's first moment
`r^3 (5/6 - pi/4)`; at a canal's leaning end the columns are integrated up
to the end plane, in `q` with `x = 20 - q^2 / 4` where the plane clips
them.

| case (r 4) | cut away (built) | reference | error | bound | worst G1 |
|---|---|---|---|---|---|
| vertex seat: top narrowing behind a cliff edge with a kink at `(14, 2.5)`, held from `x = 4` to the end | 53.507792555 | 53.507792541 | 1.4e-8 | 6.0e-8 | 2.1e-6 rad |
| two runs: two pockets r 5 beside a band 40 long | 124.215057101 | 124.215056840 | 2.6e-7 | 1.3e-8 | 1.7e-6 rad |
| a chain of three edges over a split top and front, pocket r 5 | 96.443822236 | 96.443822276 | 4.0e-8 | 6.5e-9 | 1.7e-6 rad |
| both end faces leaning (`x = y / 4`, `x = 30 + 3 z / 10`), pocket r 5 | 94.756504255 | 94.756504295 | 4.0e-8 | 6.5e-9 | 1.7e-6 rad |
| a canal cut by a leaning end face (`x = 20 + z / 4`) | 63.109947456 | 63.109947560 | 1.0e-7 | 4.6e-9 | 1.1e-6 rad |

The volume comparisons are in the heavy tier; the builds, joins and kept
edges are in the default one.

- Each result is valid and closed, every fitted piece meets its
  neighbours, the band and the front within 1e-5 rad of tangent, and the
  held edges are kept on their curves.
- The vertex seat reaches its part's end: with hold edge at the end off,
  or hold edge off, the request goes on to trim through, which refuses
  this part by name; with trim through off too it is refused naming hold
  edge at the end.
- A square pocket beside the edge (the ball on two rims at once) is not
  held: the request gets trim through's solid, the same as with hold edge
  off.

### Amendment: largest-size readings confirmed by the build (#1173)

Apps clamp a size slider with `largest_fillet_radius`,
`largest_fillet_radius_all`, `largest_fillet_distances` and
`largest_chamfer_distance`, and show the limit before the user tries a
size. Until now each read only the seat: where a contact first leaves
its face. §3.5 had them become the predicates' exact boundary, but the
predicates are not all the refusals. Below the seat a request can be
refused by a corner that does not close, a bevel reaching a concave
floor, an overflow the engine does not build, or a result failing its
check. Above it a request can build: a bevel rolls over a convex far
side as a notch (amendment "M5b"), and a contact can run off a face
between the stations the check reads. So the reading offered sizes the
call refused, and hid sizes it built.

**Rule.** A reading answers what its blend call does.

1. **The seat starts the search.** The closed-form seat reading stays,
   unchanged, and gives the first size probed. Where the seat is not read
   (a marched chamfer's crease, a walk whose side depends on the size),
   the search starts from the bound asked.
2. **The build confirms it.** The call the request makes, with the same
   options, is asked on a copy of the model: at the start, then a
   relative `LIMIT_TOLERANCE` (1e-4) past it if that built, then doubling
   up or halving down until one size builds and one is refused, then
   bisecting between the largest size that built and the smallest refused
   until they are within the tolerance. The answer is the size that
   built. The probes are `fillet_edges_corner` with the corner asked,
   `fillet_edge_distances` with the section asked, and
   `chamfer_edges_cornered` with the corner asked: the public calls
   themselves, so no refusal is left out.
3. **A range that is not one interval is reported, not hidden.**
   `LIMIT_SAMPLES` (3) sizes evenly below the answer, and as many spread
   geometrically from its refusal to the bound, are probed too. A refusal
   below moves the answer down to the boundary under that refusal. A
   build above is reported. Either way the answer is where the sizes that
   build first break, with the size that builds again beside it.
4. **The cost has a cap.** At most `LIMIT_PROBES` (64) builds. A bracket
   and bisection from any start in a range of a billion take about 45.
   Where the seat is the limit the reading takes `2 + 2 · LIMIT_SAMPLES`
   builds (8); where it is bisected from a factor of two, about 14 more.
   The cap keeps a reading's cost bounded on a seat the bracket does not
   find. If it stops the search, the answer is still a size that built,
   and the report's refused size shows the wider gap.

**API.** The four readings keep their signatures and change their
contract: each answers the confirmed size, `None` where nothing above the
confusion distance builds. Beside them, `blend_limit` takes a
`BlendSizing` (rounds with their `FilletCorner`, an asymmetric round
with its `RoundProfile`, a chamfer with its `ChamferCorner`) and returns
a `BlendLimit`: the size, the smallest size refused, a size found
building again, the seat reading and the probe count.

- The readings change in place rather than gaining `_exact` twins. An app
  clamping a slider needs the number that builds, and four seat-only
  readings left beside exact ones would be the wrong call for every
  caller.
- The report is one function, not four `_report` variants. Its request
  enum is `#[non_exhaustive]`, so a request of several bevels (a chamfer
  reading over many edges) joins as a variant without breaking callers.
  It is also the only place the call's options reach the probes.
- Errors: a usable bound and distances are still checked first. A seat
  the reading refuses (a seam, a seat no measure fits) no longer ends the
  reading: the search runs from the bound, and the seat's refusal is
  returned only where no size builds.

No reading is kept from one call to the next. A `Model` carries no
revision, so a kept answer cannot tell that its solid changed. An app
reads once per selection.

**What the search found.** Pinned in `tests/exact_limits.rs`:

- A 2 mm plate's edge bevels as a notch until its setbacks cross on the
  end faces. That happens where the bevel's plane passes the plate's far
  bottom edge, at 22 for the symmetric bevel and 20/3 + 2 for a scale of
  3 along the top and 1 down the side. `tests/largest.rs` pinned the
  seat's 2 for both. Those tests now pin the new answers, since the notch
  is what the call builds.
- The cone rim's round fails to build within 2e-4 of the exact seat
  (3.6178 against 3.6180). `tests/largest.rs` had stepped under it by
  hand; the reading now answers it.
- Three rounds at a 2 mm plate's corner, radius 2, made the one-boolean
  fallback recurse without end. The polytope it rounds handed its own
  rounds to a polytope again, the same one. A polytope's rounds are no
  longer handed on. Those rounds now build up to the seat.

**Oracles** (`tests/exact_limits.rs`). For every case the answer builds
valid a tolerance under it and is refused a tolerance over it. The cases
are boxes square and oblique (one edge, all edges, a chamfer), a drum's
rims, a bore's rim, a cone's rim, the corpus's L and U channels and
stepped block, a chopped cube's corner, and a boss's foot. At least one
case per refusal:

- setback: the plate's edge rounds to its 2 mm;
- a corner: the plate's rolling corner, refused below the setback; and
  the chopped cube's corner, stopped at 2.366, short of what any of its
  edges rounds to alone;
- a concave floor: the step's riser top bevels to its 10 mm height;
- self-intersection: the plate's bevel at 22 and 20/3 + 2, as above;
- validity: both rims of a drum 3 high, to 1.5, where the bands meet;
- overflow: the foot of a boss of radius 2 standing 2 from the block's
  edge rounds to 2, where its contact runs over that edge.

**Cost** (`tools/ogeom-bench`, one thread, release):

- `limit_box_all`: all twelve edges of a 10 mm box, where the seat is
  the limit. 131 ms for 8 builds, against 17 ms for the one build of
  `fillet_box_all`.
- `limit_drum_rim`: a drum's rim. The round meets the cap's axis below
  the seat, so the limit is bisected: 56 ms for 21 builds.

The test profile is optimised less. There the cases above take 20 to
400 ms on planes and 0.1 to 4 s on curved faces.

### Amendment: mixed rolling corners on curved faces (#1145)

Amendment "rolling corners on curved faces" left a vertex mixing convex
and concave rounds refused by name wherever a face there is curved: a
D-shaped boss's foot with its uprights, a D-shaped pocket's rim with its
uprights. Mainstream CAD closes these corners with the ball, as it does
between planes.

**What the ball does there.** Take the D-shaped boss: at a foot corner
the upright's round is convex, the foot's two rounds (along the chord's
plane and round the boss's cylinder) concave, all of radius `r`. The
upright's edge is a line square to the plate, and its faces (the chord's
plane and the cylinder) are ruled along it, so its band is a cylinder of
radius `r` about an axis square to the plate.

- The ball touching the plate and that band has its centre `r` above
  the plate and `2 r` from the band's axis: the torus of amendment
  "rolling corners", major radius `2 r`, minor `r`.
- It meets the chord's foot band where the centre circle touches that
  band's centre line, on the meridian square to the chord's plane.
- It meets the cylinder's foot band (a torus about the boss's axis) where
  the two centre circles touch. The band's axis stands `R - r` from the
  boss's (`R` the boss's radius), so the circle of radius `2 r` about it
  touches the circle of radius `R + r` about the boss's axis from inside,
  on the meridian through both axes. Where centre curves touch, the two
  sweeps share the ball's circle there, so the joint is G1.
- The corner face thus spans more than a quarter turn of the torus
  (`π/2 + asin(r / (R - r))` on the half disc), but it is still the ball's
  sweep, and every new face is a plane, a cylinder or a torus.

**How it is built.** That torus is exactly the round of the upright
band's foot on the plate: a concave edge between the band (a convex
cylinder of radius `r`) and the plate, rounded `r`. And the foot's other
two rounds run on tangent into it: the band's foot arc is tangent to
each of their edges where the band touches the face beside it. So the
lone bands are built first, then the rest of the request as tangent
chains through each lone band's foot. Between planes the same order
gives the torus of amendment "rolling corners", which the engine still
builds there directly.

**Rule.** With `FilletCorner::Rolling`, where a vertex mixing convex and
concave rounds touches a curved face:

- it must join three edges, all three in the request (the lone round
  and the two of the other sense);
- the face across from the lone round (the cap) is a plane;
- the lone round's edge is a line square to the cap;
- each face beside it is a plane containing its direction or a cylinder
  whose axis is parallel to it.

Then the lone bands go first, with every round reached from them
through vertices of one sense (whose corners are the default's: a D
boss's top loop with its uprights closes on spheres), and the rest
after, each lone band's foot on its cap added to the second pass. Any
other mixed vertex on a curved face is refused by name, naming the
vertex and the condition that fails, as is a round meeting a lone round
that is itself reached from another lone band first. The order of the
edges changes nothing: both passes are order independent. `Setback`
still refuses these corners.

Not built: a lone round whose edge is an arc or a line not square to the
cap (the sweep round its end is then a canal surface on an ellipse, not
analytic), and a lone round beside a cone.

**Oracle** (`tests/rolling_curved.rs`), sound in three orders, no patch
face, every edge between two curved faces G1 within 1e-6 rad, volume
within 1e-6:

- At `t` from the plate within the foot round, the boss's section is the
  rounded D (the D with each corner rounded `r` in plan) grown by the
  round's offset `w = r - sqrt(r² - (r - t)²)`, since the torus's section
  there is the circle of radius `r + w` about the band's axis. The rounded
  D is convex, so by Steiner's formula its area grows by `P w + π w²`.
  With `∫w dt = r²(1 - π/4)` and `∫w² dt = r³(5/3 - π/2)`, the D boss
  (radius 4, 6 high, on a 20 × 20 × 4 plate, `r = 1`) comes to
  `1600 + 6 A + P r²(1 - π/4) + π r³(5/3 - π/2)`, with `A` and `P` the
  rounded D's area and perimeter in closed form.
- With the top loop too, the top's sphere corners shrink the section the
  same way: less `P r²(1 - π/4) - π r³(5/3 - π/2)`.
- The D pocket (radius 5, 5 deep in a 30 × 30 × 10 block) with its rim
  and uprights: the void's section grows the same way within the rim
  round, `9000 - 5 A - P r²(1 - π/4) - π r³(5/3 - π/2)`.
- The foot line and one upright without the foot arc: refused by name,
  naming the vertex ("beside an edge left sharp").

### Amendment: a rim's bevel stops at a concave floor (#1145)

Amendment "a bevel stops at a concave floor" read straight edges between
planes only. A rim's bevel running down a curved wall past a concave
foot still went to M5b's tools and cut a groove into the floor: a boss
on a plate, its top rim chamfered farther down the wall than the wall is
high.

**Rule.** The same rule, on a rim's seat. A convex rim bevel (a planar
cap and a coaxial cylinder, cone or sphere wall) is refused by name,
naming its spec, where a concave edge bounding its wall or its cap has a
point strictly inside that leg's strip, within the rim's turn:

- on a cylinder or cone wall, the strip between the rim's plane and the
  contact's depth along the axis (the leg read along the wall's ruling);
- on a sphere wall, between the rim's latitude and the contact's (the
  leg read along the great circle);
- on the cap, between the rim's circle and the contact's.

The wall's outer ring is read, and on a wall in one piece round its axis
any ring of circles about the axis (its foot has a ring of its own). A
ring standing within a face is a feature the bevel reaches under, the
non-local case of amendment "M5a", as before. A leg that reaches the
concave edge exactly builds. A convex far side still rolls over (M5b): a
drum's rim chamfered past its short wall is a notch.

**Oracle** (`tests/curved_chamfer.rs`): a boss of radius 6, 1 high, on a
plate, its top rim chamfered 1.5 is refused by name; at 0.8 it builds
the boss less the bevel's triangle turned about the axis by Pappus,
within 1e-6; with 1 on the wall exactly it builds to its Pappus volume.
A dome (a ball of radius 5 standing out of a block, cut flat at 4.5)
chamfered 0.5 across the top and 4 down the sphere, past its foot on the
block, is refused by name; symmetric at 1 it builds, sound. The
quarter-disc block on a slab now refuses its rim alone as well as its
side.

### Amendment: largest-size readings for laws (#1175)

An app editing a variable round's law (its start radius, its end radius,
an interior point, or the whole law at once) needs the same limit a
constant round gets: the largest value that builds. The amendment above
gives that for one size. A law has several values, and only one is free
at a time.

**Rule.** The law reading is the constant reading with the free value as
the size. Everything else in the amendment above holds: the contract
(the size builds, a relative `LIMIT_TOLERANCE` above it is refused), the
search, the cap, and the report of a range that is not one interval.

1. **What is free.** One value of the law, the rest held as given: its
   first point, its last, or the point at an index (a linear law's points
   are its start and its end). Or the scale `k` of the whole law, every
   value times `k`. A two-distance round frees one value of either law,
   or the scale of both.
2. **The probe is the whole call.** Each size is built under the whole
   law so set: `fillet_edges_variable`, `fillet_edge_variable` or
   `fillet_edge_distances_variable`. The law is interpolated smoothly
   between its points, so the free value moves the radius along the whole
   stretch either side of it, and on a curved seat the march checks the
   chain everywhere. Building the whole call covers that; a reading of
   the free station alone would not.
3. **The seat is read along the chain.** The walks the constant seat
   reads, at both ends of the edge and every sixteenth between, give each
   station's reach: how far its contact can go before it leaves the face.
   The seat is the largest size at which the law's measure at every
   station stays within reach, bisected on the law itself without
   building. The ends are read too, since a linear law's free end is
   bound at the end. On a chain of several edges, where the stations'
   places along the chain are not read, the constant round's seat is held
   at every sixteenth of the chain.
4. **A closed chain's start is its end.** On a closed chain a law with a
   point at fraction 0 and one at fraction 1 states one point twice, and
   the build refuses the two apart. Freeing either moves both.

**API.** `BlendSizing` gains three variants, so `blend_limit` reads laws
and returns the same `BlendLimit`: `RoundLaw` (a chain, a law and a
`LawSize`), `RoundLinear` (one edge, start, end and a `LawSize`) and
`DistanceLaws` (an edge, a face, two laws, the section and a
`DistanceLawSize`). `LawSize` is `At(LawStation)` or `Scale`;
`DistanceLawSize` is `OnFace(LawStation)`, `OnOther(LawStation)` or
`Scale`; `LawStation` is `Start`, `End` or `Point(index)`. Three thin
readings return the size alone: `largest_law_radius` (a point),
`largest_law_scale` and `largest_law_distances`.

- Variants rather than a second report function. `BlendSizing` is
  `#[non_exhaustive]`, so the variants break no caller, and the search,
  the report and the cap stay in one place.
- One `LawStation` for every law, rather than start and end arguments:
  a law of points has as many values as points, and a linear law is the
  two-point case of it.
- `fillet_edge_variable` has its own variant. It is its own call, which
  falls back to the linear round, so a reading of it probes that call.
- A station the law lacks is refused before any build, and so is a law
  that states no usable radius at the bound.

**What the search found.** `fillet_edge_distances_variable` checks its
contacts at the interior sixteenths of the edge only, leaving the ends
to a face that may overhang the corner. On a cap whose depth at the
edge's end is 6, a distance law along the cap freed at that end builds,
valid, up to 6.2: where the contact at the last interior station, 15/16
of the way, leaves the cap. The reading answers 6.2, since that is what
the call builds. The round radius law on the same edge is refused at 6.
That check is the engine's; this amendment does not change it.

**Oracles** (`tests/law_limits.rs`). Each answer builds valid a tolerance
under it, is refused a tolerance over it, and stands at the closed form:

- a linear law on a straight edge between planes, the end radius free,
  on a cap whose depth runs from 4 to 6 along the edge: the depth at the
  edge's end, through `fillet_edges_variable` and `fillet_edge_variable`;
- a law of three points over a cap notched to a depth of 2 under its
  middle, the middle point free: 2, and the same for a distance law
  along the cap;
- a drum's top rim, 3 high, under a closed law whose start and end are
  one point, the start free: 3, the wall's height;
- a law of 0.5 on a 2 mm plate's edge, scaled: 4, the constant round's
  limit over 0.5; two constant distance laws of 3 and 1 scale to the
  constant pair's limit.

Every case's seat is the limit, so each reading takes 8 builds.

### Amendment: a concave band ends against its end face (#1145)

A concave straight crease between planes whose ends meet a curved wall
at a sharp corner was refused: a D pocket's floor line rounded alone,
its band ending at each end against the pocket's cylinder. The arcs
alone and every rim edge alone already built. Amendment "M5e" ends a
convex crease's band flush on the face it meets at a free end, whatever
that face's surface; nothing ended a concave one.

**Survey.** The practice already recorded for free ends (amendment
"M5e": the end face extended out to the fillet surface; section 13,
question 8: a contact that runs off its face onto the face the crease
ends at ends there) is the same for both senses in mainstream CAD's
public documentation: a fillet whose end meets a face left unrounded at
a sharp corner runs on along its spine and that face trims it. No source
caps a concave fillet in a plane square to its edge, or stops it at its
contact's end, where the wall is curved.

**Rule.** Where every asked edge is a straight concave crease between
planes, no two meeting or sharing a face, and each end meets exactly one
face past the crease's own (a wall across its way, or the solid's side
it runs out to), each fill runs on a few radii past both ends and the
end face's surface, carried on past the face, trims it, the crease's
side kept (`engine/meeting.rs`, `concave_ends`). The fill is the void's
spandrel of the crease: the wedge of its faces' planes, turned to the
void, less that wedge rounded. Its contacts must seat on both faces
between the ends; a contact running off its face past the end face is
the end face's to trim. The fill must meet the void in one piece. The
band meets the end face along their intersection, a sharp edge, as it
does where the wall is a plane.

**Crossing bands.** Concave creases meeting at a vertex of three edges
whose third is left sharp, across a turn convex in the void on the face
they share (a pocket floor's corner, its upright not asked), are each
rounded as they are alone, and the fills take their union
(`engine/alone.rs`, which did this for convex creases). This is the
crossing of amendment "rolling corners on curved faces" (a square
pocket's floor) where a band is a torus: the D pocket's line band
crosses its arc's torus band along their intersection. `Rolling` and
`Setback` agree there, as the rounds at the vertex share their sense.

**Oracle** (`tests/rolling_curved.rs`), sound, alike in three orders and
under `Setback`, no patch face:

- The D pocket (radius 5 on its chord, 5 deep in a 30 × 30 × 10 block),
  its floor line rounded 1 alone: at `s` in from the chord the fill's
  section is `h(s) = 1 - sqrt(2s - s²)` high and runs across the D,
  `2 sqrt(25 - s²)` long, integrated after `s = u²` (2e-16 relative).
- Its whole floor loop, the line and the arc: the line's fill, plus the
  arc's half torus fill by Pappus (`π (5 - x̄) (1 - π/4)`, `x̄` the
  section's centroid from the wall), less what both fill at each corner:
  in `(s, t)`, `t` in from the cylinder, both fill up to `h(max(s, t))`
  over the plan's area element `(5 - t) / sqrt((5 - t)² - s²)` (2e-16
  relative).

### Amendment: a marched seat's bevel stops at a concave floor (#1145)

Amendments "a bevel stops at a concave floor" and "a rim's bevel stops
at a concave floor" read straight edges between planes and rims. A
crease on a marched seat (a straight edge between a plane and a
cylinder along it, a spline's crease along its wall) still went to M5b's
tools: a step whose riser is a cylinder, its top edge chamfered past the
riser's foot, cut a groove into the ledge.

**Rule.** The same rule, on a marched seat. A convex bevel there is
refused by name, naming its spec, where a concave edge bounding either
face's outer ring has a point strictly inside that face's strip: between
the crease's ends, and in its section nearer the crease than the leg's
contact. Each leg is walked as the marched bevel walks it, along the
face's section square to the crease, read at 64 stations and along the
chords between them. A leg that runs off its face's surface (a swept
wall that ends at its face's edge) takes in the face's whole section:
the longest leg that still seats is read, and the rest of the leg is
added to its reach. A spec whose contacts do not read at all is passed
to the build, which reports it.

**Oracle** (`tests/curved_chamfer.rs`, `tests/chamfer_marched_open.rs`):
the cylinder step (an arc about (9, 3) of radius `sqrt(10)` from the
ledge at y = 4 up to the top at y = 6, 20 long) chamfered 3.5 is refused
by name, the foot being `sqrt(10) acos(3/5)`, about 2.93, down the arc;
at 2 it builds the section less the triangle to the riser's contact
`2 / sqrt(10)` round the arc and the arc's segment beyond that chord,
within 1e-6. The open spline slab with a block fused behind it, its wall
ending 2 from the crease at a concave foot, chamfered 2.5 is refused by
name; at 1.5 the rule passes it to the build.

### Amendment: rolling corners on arcs, and concave corners from the void's side (#1145)

Amendment "mixed rolling corners on curved faces" built the torus round
a lone band's end where the lone edge is a line square to the cap. What
it left was read again in this order: corners on outlines with arcs,
a lone edge not square to its cap, cones, curved caps.

**Corners on outlines with arcs.** Bosses and pockets whose outlines
mix lines and arcs that are not full circles, or are arcs alone, have
their mixed corners at the uprights. There the lone round is the upright:
a line square to the cap, its faces planes or cylinders along it. So the
existing rule already covers a line meeting an arc in a convex or a
concave bow, two arcs meeting (a lens), and reflex corners (two discs
fused into a figure eight, at its waist), on bosses and pockets. Tests
now pin each one.

One case did not build: a D pocket rounded all over. Its floor corners
join three concave rounds on the pocket's cylinder, the floor's line and
arc and the upright. The lone-band rule rounds those first, and the
engine had no route for three concave rounds meeting on a curved face,
though it closes the convex case (the boss's top) on the ball's sphere.
`Setback` refused them too, so this was a gap in the engine, not in the
rolling rule.

**Rule: concave rounds from the void's side.** Rounding a concave edge
with a ball is rounding the same edge of the void, where it is convex:
the closing of a solid is the complement of its void's opening. Where
every edge of a `fillet_edges` request is concave and the engine's own
routes decline, the request is taken from the void's side
(`engine/void_side.rs`):

- A box stands about the asked edges, `8.371 r` clear of each, and the
  solid is cut from it. What is left is the void beside the edges, where
  each asked edge is convex.
- The void's images of the edges are rounded by the engine. What that
  takes off the void is the fill, which is fused onto the solid.
- It is refused, and the request declined as before, where any of these
  holds: a fill reaches a face the box made (a band running that far
  would be the box's round, not the solid's), the void's rounds do not
  build, or the result is not one sound solid with one shell.
- Each asked edge's band faces are recorded as generated from it, so a
  later pass (the lone-band rule's second) finds them.

Named refusals of the request stand before this route is tried. Every
corner it builds has rounds of one sense, so `Setback` and `Rolling` build
the same solid. With it the D pocket's floor and uprights close on the
ball's sphere, and the pocket rounded all over rolls.

**A lone edge not square to its cap, cones, curved caps: still refused.**

- A lone round whose edge is an arc (a half cylinder lying on a plate,
  its end arcs convex over a concave foot) has a torus band. The ball
  rolling round that band's end keeps its centre a radius off the cap
  and twice the radius off the torus's centre circle. That path is a
  section of a torus, not a circle.
- A lone edge leaning to its cap (a frustum halved through its axis, its
  uprights along the cone's rulings) gives an ellipse.
- A curved cap (a rib below the top of a drum it is fused to) bends the
  path in turn.

In every one of these the sweep is a pipe of the ball's radius about a
non-circular centre curve. It is a canal surface, generally not rational
and not analytic.

**Survey.** Mainstream CAD's public documentation describes the
rolling-ball corner by the ball, rolled along each edge and round the
corner, and offers setbacks (one patch over the opening) as the other
corner type. Nothing restricts the ball's corner to analytic surfaces.
The literature agrees that a constant-radius rolling-ball blend is in
general irrational and is approximated where its spine is not a conic
[S2, S3]. So the corner these cases ask for is the fitted pipe about
the centre curve, G1 to the bands, not a patch. A setback patch would
be a different corner, so it is not substituted.

The engine cannot build that pipe in the lone-band order yet. The
second pass would round a tangent chain whose middle link is marched
(a plane against an inclined cylinder or a torus, with the two lines
beside it). The engine declines that chain from either side, and on a
leaning block it even declines some of the lone marched links alone.
Until that chain builds, these corners stay refused by name, naming the
vertex and the condition. `Rolling`'s promise that every new face is a
plane, a cylinder, a cone, a sphere or a torus holds for everything it
builds today. Admitting the fitted pipe will relax it for these corners
only.

**Oracle** (`tests/rolling_curved.rs`), in three orders, `check` clean,
no patch face, every edge between two curved faces G1 within 1e-6 rad.
At `t` within a round, the section of the boss or of the pocket's void
is the outline rounded `r` in plan, grown or shrunk by the round's
offset `w`. Steiner gives `A ± P w + π w²` for it, with `A` and `P`
computed independently in the test: each corner's circle is centred
where the two pieces' carriers, offset `r`, meet, and the area follows
by Green's theorem piece by piece (checked against the rounded D's
closed form to 1e-12). Measured relative errors:

- lens boss foot and uprights 5e-16, all round 5e-13;
- boss with an inward arc, foot and uprights 3e-16;
- lens pocket rim and uprights 2e-16, all over 2e-14;
- pocket with an inward arc, rim and uprights 2e-16;
- figure-eight boss top and waist 4e-16, figure-eight pocket floor and
  waist 7e-16;
- D pocket floor and uprights 1e-13, the same under `Setback`;
- D pocket all over 1e-13.

The half cylinder's foot, the half frustum's foot and uprights, and the
rib on the drum are refused by name. `Setback` refuses them too.

- [S2] Vilnius University, blending natural quadrics by a rational
  variable rolling ball (the general rolling-ball blend is not rational,
  so a rational one is built instead):
  https://epublications.vu.lt/object/elaba:3210101/3210101.pdf
- [S3] SINTEF, SAGA WP3, Improved blends between primitive surfaces
  (approximate fixed-radius rolling-ball blends where the spine is
  irrational):
  https://sintef.no/projectweb/computational-geometry/pro/saga/wp3-algebraic-geometry-for-cad-applications/improved-blends-between-primitive-surfaces

### Amendment: a marched band leaning into its end face (#1145)

Amendment "rolling corners on arcs" left a lone edge leaning to its cap
refused. Read again from the lone band's own round: a half frustum on a
plate, its uprights on the cone's rulings, rounded alone, was declined.
The upright is a marched seat (a plane and a cone along a ruling, the
ball's centre on a hyperbola), and its band stopped at the foot vertex:
the plate continues past it, the crease runs on into the material, and
the band caps in its section's own plane there. That plane leans off the
plate as the ruling does, so the cut would notch the plate on one side
of the ruling and leave a sliver standing on the other. A domain error
(the march reading the bounded ruling and cone past their ends) was all
that stopped it.

**Rule.** Where the engine's marched routes decline a request of
straight creases, it is tried once more with the bands carried past
their ends (`marched::carried_past_ends`):

- the crease is read on its whole line, and a cylinder or cone host on
  past its height range by two radii, a cone short of its apex (further,
  a cone's spine runs over the top of its hyperbola onto the next
  ruling);
- at an end where the crease's continuation runs into the material
  across one planar face, the band runs on until the ball's centre
  stands one and a half radii past the vertex, and its wedge is clipped
  by that face's plane, so the face trims the band and keeps all it
  holds;
- the end section must lie wholly past the face before the clip, or the
  request is declined as before;
- the wedge's face on the clipping plane is credited to the end face,
  not to the round: the band is the round's only face.

**Oracle** (`tests/leaning_end.rs`). The half frustum run on below the
plate is convex, so rounding its uprights is its opening: the balls of
radius `r` inside it, whose centres fill `E = {y ≥ 10 + r, ρ ≤ R(z) −
r / cos α}`. A point is kept where it lies within `r` of `E`, the
distance read in closed form to `E`'s plane, its cone and the two
hyperbolas where they meet. The volume taken off is integrated section
by section in polar form, 4096 rays by 8 Gauss panels (2048 rays agree
within 3e-8). Measured: relative error 1.8e-6, both orders alike to
1e-9, `Rolling` the same, `check` clean, the plate whole below its top,
the band G1 to the plane and the cone within 5.3e-4 rad along its rails
(the march holds its fit in distance, the bands' fit target of 2e-4).

**The rolling corner beside it stays refused.** With the uprights
built, the lone-band order rounds the foot chain next: the plate meeting
the cut face, the band and the cone. It builds, `check` clean, and the
volume is wrong, by 0.016 of fill at the two corners against a sampled
reference (the closing of the rounded frustum by balls clear of it and
the plate). The concave ball changes walls where its contact crosses the
band's rail, not where the creases meet: the rail leans, so at the
height of the contact it stands 0.35 from the crease's vertex. The chain
caps its straight and marched links at the vertex's section, and the
straight band's end face is left exposed beside a gap. For an upright
square to its cap the two places are one, which is why the D boss
closes. The corner here is the pipe of radius `r` about the curve `{z
= cap + r, dist(c, spine) = 2r}` between the sections where the contact
crosses each rail, with the straight and rim rounds ending on those
sections. Until that is built, the mixed vertex is refused by name as
before ("does not stand square"). The lone arc (a half cylinder lying on
a plate) and the curved cap (a rib under a drum's top) meet the same
transition and stay refused too.

### Amendment: a ball changing walls off the vertex (#1145)

Amendment "a marched band leaning into its end face" left the half
frustum's foot, rounded once its uprights are, 0.016 of fill short with
`check` clean: a wrong solid reachable in two calls. The cause is
general. Where two rounds of a request meet tangent at a vertex on one
shared face (the cap) and their other faces (the walls) meet G1 along a
seam there, the ball rolling along the chain changes walls where its
contact crosses the seam at its own height, at the radius from the cap.
The chain's rounds are joined in the section square to the chain at the
vertex, so where the crossing stands off that section they leave a gap
or a step beside the seam. An upright seam (the D boss's) crosses in the
section; a leaning one (a band's rail on the frustum's cut face, 0.35
off) or a curved one upright at its foot (an arc's torus band meeting
the end disc along a circle, 0.17 off) does not.

**Detection** (`engine::switch_over`). For each two rounds of one sense
and radius meeting tangent at a vertex on one shared face, with an edge
of their two walls at the vertex: walk the seam from the vertex to where
the ball off the walls' mean normal stands the radius from the cap (on
the material's side for convex rounds), and measure how far along the
chain that centre stands from the vertex's section. Beyond ten times
the bands' fit target the crossing is off the vertex. `fillet_edges`
and `fillet_edges_with` read this before any route, so no route joins
such a chain at the vertex: the request is built as below or refused by
name, naming the vertex and the distance.

**The pipe** (`engine::pipe_corner`, through the contour route). Where a
concave link's crossings at both its ends stand past their vertices over
its neighbours, the ball between them touches the cap and the link's
wall alone, its centre on the curve at the radius from both. The round
there is the pipe of the radius about that curve:

- each end's ball is solved again off the neighbour's own wall at the
  seam's point, so it stands the radius from the cap and from the seam
  exactly; its section is the plane square to both contacts, which is
  the neighbour's own section there;
- the centre is marched between the two sections in steps of a
  two-hundredth of the radius on the cap's and the wall's offsets, then
  64 stations are settled on it, closer toward the ends, each the ball's
  exact arc from its cap contact to its wall contact;
- the sections are fitted within half the confusion, and the fit is
  held to the bands' fit target of the ball about the centre's run
  halfway between them;
- the tool is the fill itself: the pipe, legs on the cap and on the
  link's wall bounded by the link's foot, the seams and the pipe's rails,
  a leg on each neighbour's wall past its vertex, and the two section
  planes;
- each neighbour's band (a straight round, a rim's torus) is trimmed
  back to its section and joined to the pipe there, the whole loop one
  fill.

A joined pair of wedges keeps its bands' sense (fill or cut): the flush
join had always made a cut, which no fill reached before.

**Rolling.** `Rolling`'s lone-band order now takes a lone line leaning
to its plane cap, beside planes, cylinders along it or cones it is a
ruling of: the lone bands go first, and the foot then rolls round them,
through the pipe.

**Oracle** (`tests/switch_over.rs`). Run on below the plate, the solid
is a plate under the opening `S = E ⊕ B_r` of a convex body. The free
centres are `F = {z ≥ cap + r} ∖ int D`, `D = E ⊕ B_2r`, and a point
outside `S` below `cap + r` is filled exactly where none of three holds:
its lift to `cap + r` lies outside `D`; it moved `2r − d` along the
normal from its nearest point of `E` reaches `cap + r`; it lies within
`r` of `∂D ∩ {z = cap + r}`. For `E` convex every other candidate
nearest point of `F` is farther than `r`. The fill is integrated in
polar sections, `z = cap + r s²`, the panels split where the cone's
contact stands. On an upright D (no slope) it gives the closed form
`p r²(1 − π/4) + π r³(5/3 − π/2)` within 2.9e-6. Measured on the half
frustum: fill within 2.3e-6 of the reference (3.438297, converged to
3e-8), the same in either order, `check` clean, every join G1 within
2.4e-4 rad, and `Rolling` in one request the same solid in either
order.

**Still refused, by name.** A convex chain (the half frustum's top loop
once its uprights are rounded, 0.19 off), and a link with one crossing
off and one in its vertex's section (a half cylinder lying on a plate,
its end arcs rounded: the arc band's foot hands over in the section on
the cylinder's side and 0.17 off it on the end disc's). The curved cap
(a rib under a drum's top) is declined by every route as before.

## 12. Acceptance corpus

From #902:

- The stepped bracket, chamfer and fillet on every edge, matching the
  reference render.
- Every edge of a box, an L-bracket, a T and a U channel, a stepped block,
  and a plate with a boss (tangent and not), at several radii and distances.
- Mixed requests: rounds with bevels, asymmetric distances (#877), variable
  radii.

On every case:

- every permutation of the request gives the bit-identical solid;
- no refusal on planar polyhedral parts inside the §3.5 domain;
- exact oracles: polyhedral volume for planar chamfers, closed-form section
  areas for prismatic rounds, the morphological reference on pure-convex and
  pure-concave regions; plus the region and no-overhang guards;
- the §7 speed targets;
- the real-world check unchanged or better.

## 13. Settled by survey

Each behavioural question is settled by what mainstream applications and
kernels do. The evidence comes from three kinds of source:

- the public user documentation of SOLIDWORKS, Autodesk Fusion, Onshape,
  PTC Creo, CATIA V5, Siemens NX, FreeCAD and Blender [P1–P16];
- the Apache-2.0 SMLib sources in usd-brep and their solidcore port [14];
- the papers.

No kernel's implementation source other than SMLib's was read. Where a
document was not found, the table says so rather than guessing.

### Question 1: a chamfer corner where several bevels meet

| Source | What it does |
|---|---|
| Fusion [P1, P2] | **Default:** the "chamfer" corner type, "a small patch … created at the corner". **Options:** miter ("the chamfer surfaces extend and intersect one another"), and blend. |
| Creo [P3] | Picks a transition by context. **Corner Plane:** three pieces only. **Patch:** three or four. **Intersect:** a sharp boundary. |
| SMLib [14] | Linear sections at a fully chamfered corner (`MakeAnalyticCornerSurface`, `sm_CreatePlanarSurface`). **Three points:** an exact plane. **Four:** a bilinear patch through them, so not flat where they aren't coplanar. **Five or more:** a Coons or n-sided patch. It never builds several flat triangles. **Separately:** its chamfer API forces mitred corners (`SetBevel`) and refuses more than three edges. |
| Blender [P4, P5] | **Default:** "grid fill", where the profile continues through the intersection. **Option:** "cutoff", a flat cap at the end of each bevelled edge. |
| SOLIDWORKS [P6] | Vertex chamfer is a separate feature. What edge chamfers make at a vertex: not found. |
| Onshape [P7]; NX; CATIA; FreeCAD | Not found. |
| #902's reference render | A triangular facet at every vertex. |

**Practice agrees:**

- Where three bevels meet, every source closes the corner with a flat face
  (Fusion's default corner, Creo's Corner Plane, SMLib's plane, the
  reference render) rather than mitring the bevels to a point. The mitre is
  an option everywhere, never the default.
- Where four or more meet, or the corner points are not coplanar, the
  corner is a patch through them (Fusion, Creo Patch, SMLib). No source
  triangulates.

**Rule:**

- The default chamfer corner is the **facet** through the setback corners
  wherever they are coplanar, always so at three bevels.
- Otherwise it is **one patch** through them, built in M4: G0, the bilinear
  patch at four corners and their Wachspress blend at five or more (see the
  amendment in §11).
- The mitre stays available as an option (`ChamferCorner::Point`).
- **API change:** today `chamfer_edges` mitres; under this rule it takes the
  facet. M1 records the change in the changelog as a change of default.

### Question 2: the default round corner and its setback

| Source | What it does |
|---|---|
| SOLIDWORKS [P8, P9] | Setback parameters are an optional section with distances the user types. No default distance. FilletXpert: where all three edges are convex or concave, or have the same radius, "only one fillet corner type exists". |
| Fusion [P10, P11] | Corner type: rolling ball or setback. Asymmetric fillets support only setback. The dialog's default: not found. |
| Creo [P12, P13] | Default, context-sensitive transitions. **Corner Sphere** for three pieces, at the largest round's radius by default, with adjustable setbacks. **Patch** for three or four. |
| NX [P14] | Corner setback and stop-short options. No default documented in what could be read. |
| CATIA [P15] | Blend corners are an opt-in reshape. Created by vertex they take "default setback distance", its value not stated. |
| Onshape [P7] | No setback option. A tech tip fixes three-edge corners by hand with a fill surface. |
| SMLib [14] | **Equal radii:** the end arcs lie on one sphere, so a sphere. **Unequal radii:** a rolling-ball blend or a surface of revolution where it fits, otherwise a Coons patch. **Setback:** off by default (`m_dSetBackDist = 0`), opt-in, and not exposed by its public API. |
| Várady–Rockwood [5, 6] | Setbacks as a quality option for vertex patches, with distances chosen per corner. |

**Practice agrees:**

- The default corner needs no setback: a rolling-ball corner, which is the
  sphere where one ball fits, and a blend patch otherwise.
- A setback is opt-in, with distances the user gives. No source documents a
  default distance.

**Rule:**

- The default round corner is the **sphere** wherever one ball touches all
  three faces (equal radii at a convex or concave trihedral vertex), as
  SMLib and Creo's Corner Sphere build it.
- Elsewhere it is **one G1 patch** bounded by the bands' natural ends, with
  no setback.
- **Setbacks** are applied only when asked, with the request's distances,
  as SOLIDWORKS, NX, CATIA and Fusion offer them.
- Creo's sphere at the largest radius for unequal rounds is a minority
  default and is not adopted; SMLib and Fusion patch there.

**Built (M4, setbacks asked).** P2 (the chamfer side's, #1038) stops each
strip square to its edge at its setback and trims each face round the
vertex on the conic tangent to both rails. P3 (the round side's) takes
those as the vertex's sides: each strip ends on its section's arc at its
setback, its corners there its rails; each face's ring passes the vertex
along its conic trim; and one transfinite patch closes the `2 n`-sided
opening, the end arcs on their strips and the trims on their faces' planes,
G1 to both. Five or more sides blend with their distances to the fourth
power: a trim's ribbon is its face's plane running on past the corner's
other faces, and cubed it entered beside a neighbouring trim at second
order (the patch stood 2e-4 past its face).

### Question 3: a chamfer and a round meeting at one vertex

| Source | What it does |
|---|---|
| Fusion [P10] | Rule fillet: fillets only. No mixed feature. |
| Creo [P3, P13] | Transitions cover pieces of the same kind only. |
| Onshape [P16] | Sheet-metal corner break cannot combine both in one feature. A second feature can round a chamfered edge. |
| SMLib [14] | Its per-edge API lets section types mix. Such a corner falls through to a Coons patch, with no test covering it. |
| SOLIDWORKS; NX; CATIA; FreeCAD | No mixed feature documented. |

**Practice agrees:** mixing takes two features applied in turn. No mainstream
application closes a mixed corner in one feature.

**Rule:**

- A request holding both kinds is two features: **the chamfers first, as one
  set, then the rounds, as one set, on that result**.
- That is the order in which `blend_edges` applies them today ("bevels
  first, then rounds"), and the order users are advised to take [P23].
- The result is the same for every order of the request. A user wanting the
  other order makes two calls.

### Question 4: edges that don't fit

| Source | What it does |
|---|---|
| SOLIDWORKS [P17] | The whole feature fails. "What's Wrong" names why. FeatureXpert can split the feature on request. |
| Onshape [P7, P18] | An error and a visualisation. The feature is marked in the list. |
| FreeCAD [P19, P20] | The operation fails ("reduce its radius"). |
| Creo [P21, P22] | Ordinary round: not found. **Auto Round** builds the rest and lists the unrounded edges, which "does not indicate that the model has an error". |
| NX [P14] | Opt-in "limit blend to avoid failure areas". |
| SMLib [14] | The executive is all-or-nothing, with per-corner and per-edge diagnostics (`SmFilletErrorInfo`). Two convenience APIs silently skip edges that can't be blended. |

**Practice agrees:**

- The default is all-or-nothing, with the failing edges named (SOLIDWORKS,
  Onshape, FreeCAD, SMLib's executive).
- Building the rest and reporting the skipped edges is an explicit,
  separate mode (Creo Auto Round, NX's opt-in limit).

**Rule:**

- Every blend call is **all-or-nothing**, its refusal naming the edges and
  the predicate that failed (as `*_checked` does).
- `*_partial` is the **opt-in mode**: build every requested edge that the
  §3.5 predicates admit, and return the skipped edges with their reasons.
- **Which edges are skipped** depends on the request as a set, not its
  order:
  1. take the failing predicate whose edges come first by stable key;
  2. drop those edges;
  3. evaluate again, until the rest fit.

### Question 5: rounds running on into each other with different sections

Two rounded edges meeting tangent at a vertex (one straight edge split by a
third edge, or a straight edge running into a rim) whose sections differ
there: different radii, or the same radius set back differently by a third
edge's round on one face.

| Source | What it does |
|---|---|
| Creo [P13] | **Blend** between two collinear round pieces: a fillet surface between them, then a stop. **Continue:** the round geometry extended across. The default is chosen by context. |
| SOLIDWORKS [P24, P9] | Variable size fillets: a **smooth transition**, the fillet changing from one radius to the other matching its neighbour's tangency, or a straight one. Setbacks make "a smooth transition ... into the fillet corner". |
| SMLib [14] | The Nx2 tangent case "might have to put another fillet surface in between the two that are there now" (`AdjustTangentCorner`). |
| Fusion; NX; Onshape; CATIA | Not found. |

**Practice agrees:** two such rounds are joined by a transition surface
tangent to both, never left stepped or cut square.

**Rule:**

- Where the two sections agree at the vertex, one strip runs on through
  it (today's case).
- Otherwise the two strips end at their natural ends and **one G1 patch**
  joins them, bounded by their end sections and by each face's setback
  curves between the two contacts, through whatever the third edge's
  round or sharp edge leaves there.

### Question 6: a convex and a concave round beside an edge left sharp

| Source | What it does |
|---|---|
| SMLib [14] | `SmFillet3x2MixedConvexityCorner`: the two rails do not meet; a **three-sided patch** goes between them, its third side along the sharp edge's face. |
| SOLIDWORKS [P25, P26] | FilletXpert corners need three filleted edges; what two fillets beside a sharp edge make: not found. |
| Creo [P13] | Stop cases chosen by context; nothing specific to mixed convexity found. |
| Fusion; NX; Onshape; CATIA | Not found. |

**Rule:** the two strips end at their natural ends and **one G1 patch**
closes the corner between them, bounded by both ends and the sharp edge's
faces' setback curves, as SMLib builds it. It is proved against the
material set out under question 7.

### Question 7: mixed convexity where three or more rounds meet (saddle corners)

| Source | What it does |
|---|---|
| SOLIDWORKS [P25, P26] | Three fillets of one convexity, or of one radius: "only one fillet corner type exists". Mixed: alternatives, a **triangular patch**, a **quadrilateral patch**, or two triangles splitting a quadrilateral, chosen by which two edges share convexity and radius. |
| Creo [P13] | **Patch** for three or four round pieces, whatever their convexity. **Round Only 1 / 2** have variants for three pieces of mixed convexity. Corner Sphere only where three pieces overlap. |
| SMLib [14] | Every fully filleted vertex with a concave edge is a `SmFilletConcaveNxNCorner`: one patch through the rails, with an opt-in setback. |
| Fusion; NX; Onshape; CATIA | Not found. |

**Practice agrees:** no source puts a sphere at a mixed corner (no one ball
touches every face from one side). Each closes it with one patch through
the strips' ends; the split or quadrilateral layouts are options.

**Rule:**

- A vertex where rounds of mixed convexity meet closes on **one G1
  n-sided patch** through the strips' natural ends, as unequal radii do.
- No class of mixed convexity is refused by name: every source above
  builds them.

**The material a mixed corner is proved against (P5).** A corner's patch
may stand only where the solid was or where a concave round fills:

- the original solid's material at the vertex: around the vertex's cycle
  of faces, the half-spaces of faces joined by convex edges intersect and
  the runs joined by concave edges unite. Three faces give the table
  `corner_material` already proves against: all convex, every face; one
  concave edge, the opposite face and either of the edge's; two, the
  opposite face or both of the convex edge's; three, any face;
- or a concave round's fill: outside both of its faces, within its
  contacts.

Where a face is tangent to the next at the vertex (a boss running into a
plate's side, the cusp of #900), the two count as joined by an edge of
the convexity the faces either side give them, so the cusp is proved
within the plate's corner, the boss above the top, or the foot's fill.

### Question 8: a round whose contact runs off its face (roll-over)

| Source | What it does |
|---|---|
| Parasolid family: T-FLEX [P27], Solid Edge [P28], NX [P29] | Overflow is resolved in a fixed order, every option on by default: **smooth overflow** (the round rolls across a tangent edge onto the next face), then **retain cliff edge** (the round is "formed by rolling a ball that is tangent to one of the original faces and supported by the retained edge"), then **notch** (the round is kept and the element it overflows is trimmed or extended to meet it). Smooth edges are resolved before sharp ones. |
| SOLIDWORKS [P30] | Overflow type **Default** picks **Keep edge** (the model's edge stays, the fillet surface splits into pieces and its edge may dip) or **Keep surface** (the fillet stays smooth and the adjacent edge is changed) from the convexity of the edges met. |
| Onshape [P7] | **Allow edge overflow** (on by default) lets the fillet change the edges it meets; off, or with edges picked to keep, the edges stay and the fillet adapts. |
| SMLib [14, P31] | Its filleting framework's default rolls the fillet over onto the adjacent surface; a cliff rollover is trimmed by the cliff face's extension (`CreateCliffRollOverSideFG`). Its handler can instead cut the fillet, **roll a ball along the edge**, or quit. "The ball may actually roll onto an adjacent surface or onto adjacent edges." |
| Fusion; Creo; CATIA | Not found. |

**Practice agrees** that a round meeting an edge it was not seated on
either rolls on that edge (the ball resting on it, the edge kept) or keeps
its own surface and changes the neighbour. Every source offers the first;
the Parasolid family tries it before notching, and SOLIDWORKS picks it for
the convexities a boss's foot or a shelf's edge give. The second changes a
face the request did not name, which is the non-local family (§4.5).

**Rule: the rolling ball rests on whatever of the boundary it meets.** The
ball of a round touches the solid's boundary on each side at a *seat*: a
face, an edge of that face, or a vertex. Its centre lies where the two
seats' offsets meet (the face's offset surface, the pipe about the edge,
the sphere about the vertex), and the round is the envelope of those
balls, piece by piece. A seat changes where the ball touches three seats
at once, a point found by one intersection:

- the contact on a face reaching the face's boundary edge: the seat moves
  to the edge (a cliff edge, or a shelf's rim);
- the ball on an edge coming down on the face beyond: the seat moves to
  that face (the edge tangent, the piece on the edge has no length: smooth
  overflow);
- the contact on an edge reaching the edge's end: the seat moves to the
  vertex.

At each change the two offsets the centre leaves and joins are tangent, so
the pieces meet G1 on one section. This is the one end condition of §4.3
for a contact that leaves its face: P2 finds the crossing and labels it
`RollOver` (face, edge, point); the strip runs on, re-seated, instead of
ending.

**Proof (P5).** Each piece is proved as the strips are, within the material
it stands in or the fill a concave round makes, and the run is proved
**clear**: every ball of it touches the solid only at its seats (none
holds anything nearer than its radius), checked on the spine with the
distance's Lipschitz bound between stations. A run that fails is refused
by name, as is a seat change this milestone does not build (a vertex seat,
a ball holding an edge on both sides at once).

### Decided by the maintainer

- **Non-local blends** (a strip reaching a face that is not its neighbour)
  are refused by name until M5.
- **Special-case work on the current path** is limited to regressions until
  M2. #877 continues on it.

## References

1. J. Vida, R. R. Martin, T. Várady. "A survey of blending methods that use
   parametric surfaces." *Computer-Aided Design* 26(5):341–365, 1994.
   doi:10.1016/0010-4485(94)90023-X.
2. B. K. Choi, S. Y. Ju. "Constant-radius blending in surface modelling."
   *Computer-Aided Design* 21(4):213–220, 1989.
   doi:10.1016/0010-4485(89)90046-8.
3. J. R. Rossignac, A. A. G. Requicha. "Offsetting operations in solid
   modelling." *Computer Aided Geometric Design* 3:129–148, 1986. Rounding
   and filleting as compositions of offsets.
4. S. Y. Ju. *Construction of blend surfaces in surface modelling.*
   Doctoral thesis, KAIST, 1989. Constant- and variable-radius edge blends.
5. T. Várady, A. Rockwood. "Vertex blending based on the setback split." In
   *Mathematical Methods for Curves and Surfaces*, Vanderbilt University
   Press, 1995, pp. 1–16.
6. T. Várady, A. Rockwood. "Geometric construction for setback vertex
   blending." *Computer-Aided Design* 29(6):413–425, 1997.
7. T. Várady, C. M. Hoffmann. "Vertex blending: problems and solutions." In
   *Mathematical Methods for Curves and Surfaces II*, 1998, pp. 501–527.
8. P. Salvi, T. Várady, A. Rockwood. "Ribbon-based transfinite surfaces."
   *Computer Aided Geometric Design*, 2014.
9. P. Salvi, T. Várady. "A multi-sided generalization of the C⁰ Coons
   patch." arXiv:2002.11347, 2020. T. Várady, P. Salvi, Gy. Karikó. "A
   multi-sided Bézier patch with a simple control structure." *Computer
   Graphics Forum*, 2016, doi:10.1111/cgf.12833.
10. J. A. Gregory. "N-sided surface patches." In *The Mathematics of
    Surfaces*, Oxford University Press, 1986.
11. P. Charrot, J. A. Gregory. "A pentagonal surface patch for computer
    aided geometric design." *Computer Aided Geometric Design* 1(1), 1984.
12. I. C. Braid. "Non-local blending of boundary models." *Computer-Aided
    Design* 29(2):89–100, 1997.
13. I. Stroud. *Boundary Representation Modelling Techniques.* Springer,
    2006. Blending in b-rep modellers.
14. SMLib fillet design, read from the Apache-2.0 sources in `usd-brep`
    (`source/SMLib`: `SmFilletExecutive`, `SmFilletCorner`,
    `SmFilletSolver`, `SmFilletIntersector`) and their Rust port in
    `solidcore` (`crates/solidcore/src/blend/`). Concepts taken: a whole-model
    pseudo-topology built before geometry; phased corner and edge passes;
    the offset-surface intersection as the section solver; pluggable
    cross-section generators; N × M vertex classes with a bevel-or-blend
    choice; "way long" sheets trimmed by intersection; one merge at the end.
15. A. P. Rockwood, J. C. Owen. "Blending surfaces in solid modeling." In
    *Geometric Modeling: Algorithms and New Trends*, SIAM, 1987.
16. C. M. Hoffmann, J. Hopcroft. "The potential method for blending
    surfaces and corners." In *Geometric Modeling: Algorithms and New
    Trends*, SIAM, 1987.

### Practice sources (public user documentation)

- [P1] Autodesk Fusion, Chamfer reference:
  https://help.autodesk.com/cloudhelp/ENU/Fusion-Model/files/SLD-REF-CHAMFER.htm
- [P2] Autodesk Fusion API, ChamferCornerTypes:
  https://help.autodesk.com/cloudhelp/ENU/Fusion-360-API/files/fusion_ChamferCornerTypes.htm
- [P3] PTC Creo 9, About Chamfer Transitions:
  https://support.ptc.com/help/creo/creo_pma/r9.0/usascii/part_modeling/part_modeling/About_Chamfer_Transitions.html
- [P4] Blender manual, Bevel modifier:
  https://docs.blender.org/manual/en/latest/modeling/modifiers/generate/bevel.html
- [P5] Blender manual, Bevel edges:
  https://docs.blender.org/manual/en/latest/modeling/meshes/editing/edge/bevel.html
- [P6] SOLIDWORKS 2024, Chamfer PropertyManager:
  https://help.solidworks.com/2024/english/SolidWorks/sldworks/r_chamfer_pm.htm
- [P7] Onshape, Fillet: https://cad.onshape.com/help/Content/fillet.htm;
  Chamfer: https://cad.onshape.com/help/Content/chamfer.htm; tech tip on
  blending fillets with Fill Surface:
  https://www.onshape.com/en/resource-center/tech-tips/how-to-blend-fillets-using-the-fill-surface-feature
- [P8] SOLIDWORKS 2024, Constant Size Fillets (setback parameters):
  https://help.solidworks.com/2024/english/SolidWorks/sldworks/r_constant_size_fillets.htm
- [P9] SOLIDWORKS 2024, FilletXpert overview:
  https://help.solidworks.com/2024/english/SolidWorks/sldworks/c_FilletXpert_Overview.htm
- [P10] Autodesk Fusion, Fillet reference:
  https://help.autodesk.com/cloudhelp/ENU/Fusion-Model/files/SLD-REF-FILLET.htm
- [P11] Autodesk Fusion API, FilletFeatureInput.isRollingBallCorner:
  https://help.autodesk.com/cloudhelp/ENU/Fusion-360-API/files/fusion_FilletFeatureInput_isRollingBallCorner.htm
- [P12] PTC Creo 9, About Round Transitions:
  https://support.ptc.com/help/creo/creo_pma/r9.0/usascii/part_modeling/part_modeling/About_Round_Transitions.html
- [P13] PTC Creo 9, Round Transition Types:
  https://support.ptc.com/help/creo/creo_pma/r9.0/usascii/part_modeling/part_modeling/Round_Transition_Types.html
- [P14] Siemens NX, Edge Blend: help copy at
  http://www2.me.rochester.edu/courses/ME204/nx_help/en_US/tdocExt/content/a/create_features_feat_edgeblend_dia.xml
  (read through search snippets; the site is offline), and
  https://learnnx.com/lesson/siemens-nx-edge-blend-stop-short-corner-setback/
- [P15] CATIA V5, Corner reshape (help copy):
  https://catiahelp.azurewebsites.net/English/Fm1UserMap/pdg-c-CornerReshape.htm
- [P16] Onshape, Corner break (sheet metal):
  https://cad.onshape.com/help/Content/sheetmetalcornerbreak.htm
- [P17] SOLIDWORKS 2024, What's Wrong / FeatureXpert:
  https://help.solidworks.com/2024/english/SolidWorks/sldworks/HIDD_NEW_WHATS_WRONG_AUTOFIX.htm
- [P18] Onshape, Error indicators:
  https://cad.onshape.com/help/Content/errorindicators.htm
- [P19] FreeCAD wiki, PartDesign Fillet: https://wiki.freecad.org/PartDesign_Fillet
- [P20] FreeCAD wiki, Part Fillet: https://wiki.freecad.org/Part_Fillet
- [P21] PTC Creo 9, About the Auto Round Feature:
  https://support.ptc.com/help/creo/creo_pma/r9.0/usascii/part_modeling/part_modeling/About_the_Auto_Round_Feature.html
- [P22] PTC Creo 9, Using Troubleshooter to View Unrounded Edges:
  https://support.ptc.com/help/creo/creo_pma/r9.0/usascii/part_modeling/part_modeling/Using_Troubleshooter_to_View_Unrounded_Edges.html
- [P23] Onshape forum, chamfer before fillet (user advice, not
  documentation):
  https://forum.onshape.com/discussion/21699/how-can-i-put-a-fillet-tail-on-my-chamfer-edge-like-this

- [P24] SOLIDWORKS 2024, Variable Size Fillets:
  https://help.solidworks.com/2024/English/SolidWorks/Sldworks/r_Variable_Size_Fillet.htm
- [P25] SOLIDWORKS 2021, FilletXpert Corner PropertyManager:
  https://help.solidworks.com/2021/english/SolidWorks/sldworks/hidd_fillet_mgr_corner.htm
- [P26] SOLIDWORKS API, IFeatureManager::FilletXpertMakeCorner (the
  corner types):
  https://help.solidworks.com/2018/English/api/sldworksapi/SolidWorks.Interop.sldworks~SolidWorks.Interop.sldworks.IFeatureManager~FilletXpertMakeCorner.html
- [P27] T-FLEX CAD 17, Overflow Processing of Edge Blend:
  https://tflex.com/help/eng/T-FLEX%20CAD/17/blendedges_overflowprocessing.htm
- [P28] Solid Edge (Parasolid), rounding overflow processing order:
  https://soliddna.com/SEHelp/ST5/EN/i_v/roundingoverflow1a.htm
- [P29] Siemens NX Open, EdgeBlendBuilder (the overflow resolution
  options: roll over smooth edges, roll onto edges, maintain blend and move
  sharp edges):
  https://docs.plm.automation.siemens.com/data_services/resources/nx/12/nx_api/custom/en_US/nxopen_python_ref/NXOpen.Features.EdgeBlendBuilder.html
- [P30] SOLIDWORKS 2022, Example: Overflow Type:
  https://help.solidworks.com/2022/english/SolidWorks/sldworks/c_example_overflow_type.htm
- [P31] SMLib manual, Fillets: https://docs.nvidia.com/smlib/manual/smlib/fillets
- [P46] CadQuery, class reference (`Workplane.clean`, the `clean`
  parameter of `union` and `cut`, and `occ_impl.shapes.clean`):
  https://cadquery.readthedocs.io/en/latest/classreference.html
- [P47] build123d, direct API reference (`Shape.clean`):
  https://build123d.readthedocs.io/en/latest/direct_api_reference.html

No implementation source of any other kernel was read for this record.
