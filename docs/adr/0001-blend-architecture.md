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
  application refuses a fillet whose faces cannot hold it. A request of two
  edges or more that no route builds is checked as the current path checks
  it before building any band (its corners read, each ball seated, no two
  bands overlapping), within the current path's cut-offs, so a concave
  band running past its face's side is kept for the cut-off at that side
  as there; one failing is refused there in the same words, without
  entering the current path. The check refuses nothing the current path
  builds, since that path refuses on it before building anything. A
  request of one edge is checked as its entry checks it (`request_fits`
  through `fillet_edges`, the ball on its two faces through
  `fillet_edge`), and where that fails, the strip or ring beside it is
  read as taken off whole, as the current path reads it; only where that
  reads nothing either is it refused, in the same words.

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
  Where the rings decline (a bevel rolling down a front whose foot meets
  the floor at a concave edge), the union of the tools is applied.
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
the wall's cut, a vertex seat this milestone refuses (§13, question 8);
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
  beside a rim, and needs no seated booleans.
- **Checked.** While the engine's marched seats run, reaching an old
  edge constructor (one that builds and applies its own wedge, or reads a
  strip off whole for the one-edge entries) is logged for the old-path
  count and panics in a debug build, so every marched test fails if the
  engine goes back to old code.

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
- **Domain.** Built where the earlier band is a cylinder seated on planes
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
  both sides (the trimmed assembly's S), as before.
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
  - a vertex where rounds meet that touches a curved face;
  - a vertex that joins other than three edges;
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

### Amendment: a cusp's S, one rule for one side or two (#902)

Where a flat side runs tangent into a drum's wall, the top between them
narrows to a cusp at the tangent point T. A round of the top's edge with the
side (convex) and of its rim on the wall (concave), leaving T on one side,
cannot stand side by side within two radii of T: the round there is an S,
off the side in a convex arc and onto the wall in a concave one, closing to
nothing at T. Both sides of T may be asked (a plate's side tangent to a
boss: two S's) or one, where the top ends at T (an arm's side running into
a hub: one S).

- **Rule.** Each side of T asked with both rounds is one lobe, built alike
  whether the other side is asked or not. The two bands end in the
  meridian plane through the point where their contacts cross on the top,
  each cut by it, so the lobe turns smoothly from the one band into the
  other. Below the lobe the S is read section by section square to the
  straight edge: where the top stands open `u` of what it does at the
  lobe's end on a face (open: from the side's line out to the wall's foot
  along the radius, which reads at every station, also where the lobe's
  end on the side lies past the foot's reach along the top, as it does
  for the wider rounds), the S turns `t` off each face, `cos t = 1 - u`,
  its foot on the side `sin t` of the lobe's depth below the top and on
  the wall `sin t` of its height above. Both feet run from the lobe to T, each meeting its
  band's contact tangentially at the lobe.
- **Built** as one patch per lobe over the hole of the lobe and its two
  feet, the feet meeting at the patch's pole, T. The lobe runs from a
  convex band into a concave one, which one polynomial patch cannot follow:
  the patch is fitted on each side of the line from that join to the pole
  and joined there C1 on a knot line, one face. Edges the S meets at T (the
  side's edge with the wall below a top ending there, the wall's seam) run
  on to T. Held to the bands' fit target, gap and tangency.
- **Refused by name:** one lobe asked while the top runs on past T with its
  edges left sharp (amendment "rounds either side of a cusp").

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
- **Not this route:** groups whose tools touch or overlap (their bands
  meet and the meeting is the routes' to close), a group the engine does
  not round alone, or a request of one group.

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
