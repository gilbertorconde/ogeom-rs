# The ogeom data model

**Status: normative.** Everything in `crates/` implements this document. A
change here is a design change and must be argued as one (see
`CONTRIBUTING.md`).

---

## Why this document exists

A B-rep kernel's public API is a thin layer over its data model. What the model
can express limits what the kernel can do. Two examples:

- Flatten a location chain into a 4×4 matrix, and assembly instancing becomes
  impossible.
- Give an edge one curve instead of a list, and boolean face splitting has
  nothing to split with.

Neither can be fixed later at the API layer.

The field has converged on the model below over thirty years, because the
alternatives do not work. It is described here in its own terms, not as a port
of any implementation. Each invariant is:

- **cheap now:** a few days of design attention;
- **effectively impossible to retrofit:** adding it later means reworking every
  algorithm in the kernel;
- **load-bearing:** its failure mode is written out, so the cost of dropping it
  is concrete.

In two places the conventional design is wrong and we deliberately differ:
stable identity (§8) and predicate abstraction (§9).

Notes marked *Elsewhere* give the conventional name for a concept, because that
is how the field talks about it. They are a glossary, not a dependency: ogeom
links against no existing CAD kernel. Where a note cites usage counts, they come
from `api_surface.json`, a profile of one large application. The counts show a
requirement is real; they do not define scope. See `SCOPE.md`.

---

## 1. A shape is a triple

```rust
pub struct Shape {
    tshape:      TShapeId,      // arena key: the shared, positionless topology node
    location:    Location,      // where this instance sits
    orientation: Orientation,   // which way its boundary faces
}
```

`Shape` is cheap to copy and is passed by value everywhere. The heavy data
(children, geometry, tolerances) is stored once, in an arena, behind
`TShapeId`.

This separation is why B-rep scales. The same `TShape` can appear at many
locations with many orientations without copying any geometry.

**Consequence: traversal composes.**

- A sub-shape's *effective* location is the product of every location from the
  root down to it.
- Its *effective* orientation is the composition of every orientation on that
  path.

An explorer that yields sub-shapes without composing both is wrong. It produces
plausible-looking garbage, not a crash.

> *Elsewhere:* `TopoDS_Shape` = `{Handle(TopoDS_TShape), TopLoc_Location, TopAbs_Orientation}`.

---

## 2. Location is a chain, not a matrix

```rust
pub struct Location {
    // (datum, power) pairs. Empty == identity.
    chain: SmallVec<[(DatumId, i32); 2]>,
    // composed transform, computed on demand
    cached: OnceCell<Trsf>,
}
```

A `Datum` is a reference-counted rigid transform. A `Location` is a *sequence*
of `(datum, integer power)` pairs.

**Why not a 4×4 matrix:**

- **Composition is concatenation.** Roughly O(1), no matrix multiply, and no
  drift from repeated floating-point composition.
- **Identity comparison is structural.** Two shapes are at "the same place" if
  their chains are equal. No need to compare 16 floats against a tolerance.
  This is what lets 10,000 identical bolts in an assembly share one piece of
  geometry *and* be recognized as instances of it.
- **Inverses are exact.** Negate the power instead of inverting a matrix.

The composed `Trsf` is computed lazily and cached. Nothing outside `ogeom-topo`
should ever need to read the chain itself.

> *Elsewhere:* `TopLoc_Location`, a linked list of `(datum, power)` pairs.

---

## 3. Orientation composes multiplicatively

```rust
pub enum Orientation { Forward, Reversed, Internal, External }
```

| Value | Meaning |
|---|---|
| `Forward` | The material is on the surface's default side. |
| `Reversed` | The material is on the other side. |
| `Internal` | The boundary lies *inside* the material (for example, a stiffener edge embedded in a face). |
| `External` | The boundary lies outside the material (reference geometry). |

Composition is a monoid, applied at **every** level of descent:

```
compose(Forward,  x) = x
compose(Reversed, Forward)  = Reversed
compose(Reversed, Reversed) = Forward
compose(Internal, _) = Internal
compose(External, _) = External
```

An edge's orientation *within a face* depends on that face's orientation
*within its shell*, which depends on the shell's orientation within the solid.
Reversing a solid must not require touching any child.

> *Elsewhere:* `TopAbs_Orientation` and `TopAbs::Compose`.

---

## 4. Identity has three levels

There are three distinct equivalences, with three distinct hashers. They are
**not** interchangeable:

| Predicate | Compares | Used for |
|---|---|---|
| `is_partner` | tshape only | "Is this the same underlying topology, anywhere, either way round?" |
| `is_same` | tshape + location | Set membership in most algorithms; the common case |
| `is_equal` (`==`) | tshape + location + orientation | Exact identity; ordered containers |

Every map and set type names the equivalence it uses and enforces it in its
hasher. A `HashMap` keyed on `is_equal` semantics that hashes only the tshape is
a silent correctness bug.

> *Elsewhere:* `IsPartner` / `IsSame` / `IsEqual` and the `ShapeMapHasher`
> family. Mixing them up is a common, well-documented source of bugs in
> applications built on kernels that expose all three.

---

## 5. Tolerances are per entity, and they only grow

Every vertex, edge and face carries its **own** tolerance: the radius of the
sphere, pipe or slab within which the entity is considered to lie.

**Containment rule** (an invariant), for entities in a boundary relationship:

```
tol(vertex) >= tol(edge) >= tol(face)
```

- Operations may only *increase* tolerances, never silently decrease them.
- An operation that cannot satisfy the rule has failed and must say so.

This is not a workaround for sloppy code. Exact arithmetic cannot represent the
intersection curve of two NURBS surfaces (it is transcendental). A
tolerance-carrying topology is the only known way to build a kernel whose
models close. Every production kernel works this way. See §9 for what exact
predicates *can* do.

> *Elsewhere:* a per-entity `Tolerance` on the vertex, edge and face records; a
> validity checker that enforces the rule; and a boolean post-pass that
> increases tolerances until the rule holds.

---

## 6. An edge carries a list of representations

An edge has a list of representations, not one curve:

```rust
pub enum EdgeRepr {
    Curve3d      { curve: CurveId, location: Location, range: (f64, f64) },
    PCurve       { curve: Curve2dId, surface: SurfaceId, location: Location },
    PCurveClosed { curve: Curve2dId, curve2: Curve2dId, surface: SurfaceId, location: Location },
    Polygon3d    { polygon: PolygonId, location: Location },
    PolygonOnTri { polygon: PolygonOnTriId, triangulation: TriangulationId, location: Location },
}
```

One edge holds, at the same time:

- a 3D curve;
- **one pcurve per adjacent face**;
- two pcurves where it is a seam on a closed surface;
- cached discretizations.

**Why it must be a list:**

- A boolean splits faces in 2D parameter space. Without a pcurve on each face,
  there is nothing to split *with*.
- Different surfaces have different parameterizations, so one 2D curve cannot
  serve two faces.
- A seam edge on a cylinder appears at both `u = 0` and `u = 2π`. One pcurve
  cannot express that.

**The `same_parameter` flag** asserts that all representations agree on the
parameterization: for the same `t`, `curve3d(t)` and `surface(pcurve(t))` are
the same point within tolerance. The claim can be false. A repair routine
re-establishes it, possibly by increasing the edge's tolerance.

> *Elsewhere:* `BRep_CurveRepresentation` and its subclasses; `SameParameter`.

---

## 7. Every operation emits history

```rust
pub trait Operation {
    fn generated(&self, input: Shape) -> &[Shape];  // new entities made *from* input
    fn modified(&self, input: Shape) -> &[Shape];   // what input *became*
    fn is_deleted(&self, input: Shape) -> bool;     // input has no image in the result
}
```

History is not optional, not deferred, and not "added when something needs
it". Every operation in `ogeom-algo`, `ogeom-bool`, `ogeom-fillet` and
`ogeom-offset` fills these in from the first commit that introduces it.

The reason is downstream. A parametric application records "fillet *that*
edge" and must still find that edge after the model is rebuilt with different
dimensions. It finds it again by walking history. This is the topological
naming problem. Every application built on a kernel that identifies topology by
pointer has had to solve it this way; in one well-known case, it took a decade
of work layered on the kernel's history maps.

Adding history later means revisiting every algorithm. Incomplete history is
worse than none: it fails silently and corrupts documents instead of raising an
error.

> *Elsewhere:* `Generated` / `Modified` / `IsDeleted` on the operation base
> class, plus a standalone history object.

---

## 8. Entity identity is stable (*deliberate divergence*)

The conventional design identifies topology by pointer. Every modelling
operation allocates new nodes, so every reference into a previous result is
lost. That *is* the topological naming problem. Every downstream fix tries to
reconstruct identity afterwards by walking history maps.

We record identity at creation instead:

```rust
pub struct EntityId(u64);         // stable for the lifetime of a document

pub enum Provenance {
    Primitive { op: OpId, role: PrimitiveRole },   // "the +Z face of box #3"
    Derived   { op: OpId, from: SmallVec<[EntityId; 2]>, role: DerivedRole },
    Imported  { file: FileId, external: ExternalRef },
}
```

An entity's identity is *what produced it, and from what*, not where it sits in
memory.

- When a boolean splits a face, each fragment knows it came from that face.
- A rebuild with different parameters produces entities with the same
  provenance, so a reference like "the fillet on this edge" survives.

History (§7) is still required. An embedding application consumes it, and it is
the honest answer where provenance alone cannot resolve a reference. But
provenance is the primary mechanism, and it only works if it is designed in
from the start.

---

## 9. Numerics are abstracted; the tolerance model is not negotiable

```rust
pub trait Predicates {
    fn orient3d(a: Point, b: Point, c: Point, d: Point) -> Sign;
    fn insphere(a: Point, b: Point, c: Point, d: Point, e: Point) -> Sign;
    // ...
}
```

Algorithms are written against this trait. Implementations may be
fast-filtered floating point, adaptive exact (Shewchuk), or interval-based, and
can change without touching any algorithm.

**What this does and does not buy:**

- Exact predicates solve the *polyhedral* robustness problem: the orientation of
  a point against a plane, in-sphere tests for Delaunay.
- They do **not** solve the CAD problem. The intersection curve of two curved
  surfaces has no exact representation to be exact *about*. That is why §5
  exists and cannot be traded away.
- Predicates decide exactly what can be decided exactly. Tolerances handle the
  rest.

Tolerance constants live in `ogeom-core`. The model's unit scale is
**explicit**, not assumed to be millimetres. Kernels often assume millimetres,
and then misbehave silently on models authored in metres or inches.

| Constant | Value at unit scale | Meaning |
|---|---|---|
| `CONFUSION` | `1e-7` | two points are the same point |
| `ANGULAR` | `1e-12` | two directions are parallel |
| `INTERSECTION` | `CONFUSION * 1e-2` | intersection convergence |
| `APPROXIMATION` | `CONFUSION * 1e1` | curve/surface fitting target |
| `P_CONFUSION` | `CONFUSION * 1e-2` | parametric-space confusion |

---

## 10. Geometry is reached through traits

```rust
pub trait Curve3d {
    fn range(&self) -> (f64, f64);
    fn value(&self, u: f64) -> Point;
    fn d1(&self, u: f64) -> (Point, Vector);
    fn d2(&self, u: f64) -> (Point, Vector, Vector);
    fn continuity(&self) -> Continuity;
    fn kind(&self) -> CurveKind;              // for analytic fast paths
    // ...
}
```

Every intersection, projection, extrema and tessellation algorithm is written
against `Curve3d` / `Curve2d` / `Surface`, never against a concrete type. To a
caller, a face adaptor (surface + location + trimming) and a bare analytic plane
are the same thing.

`kind()` lets algorithms *opt into* analytic fast paths (plane/plane
intersection should not go through a marching intersector). The general path
never needs to know what it is looking at.

> *Elsewhere:* the `Adaptor` family. This is the best idea in the conventional
> design, and we adopt it as is.

---

## 11. Memory: arenas, not reference counting

Topology lives in typed index arenas (`Vec<T>` plus a typed `u32` key), not
behind `Arc` or an intrusive reference count.

- No reference cycles to reason about. Intrusive reference counting has no
  cycle collection, so kernels built on it avoid cycles by convention only.
- Cache-friendly traversal.
- Keys are small, `Copy`, and serializable.
- It is what makes §8 possible at all.

The cost: a shape only has meaning relative to the arena that owns it. That is
the right trade, and the API makes it explicit instead of hiding it.

### Append-only arenas, and what relies on it

In practice the arenas are append-only: nothing in the kernel removes entries.
Two non-builder paths rely on this.

- **`Model::from_parts`** assembles a restored document by replaying the file's
  insertion order, which reproduces every handle.
- **`Model::absorb`** uses the same engine on a model that already contains
  data. Another document's parts are appended, with every handle shifted past
  what the target already holds, so an absorbed shape is indistinguishable
  from one built in place. This is how a serialized tool body meets a live one
  in a boolean.

Absorption:

- preserves the source's identities under a plain offset (the remap table says
  where each one landed);
- keeps the source's provenance verbatim, including source `OpId`s, which only
  have meaning in the source document's own rebuild;
- refuses, by name, a document at a different unit scale. Rescaling is a
  feature, not a default.

> *Elsewhere:* a transient base class with intrusive reference counting, plus a
> custom small-block allocator. We need none of it.

---

## 12. Errors are values

Operations return `Result<T, OgeomError>`. The variants cover the failures a
kernel actually needs: construction, domain, range, dimension mismatch, null
object, not done, numeric failure, invariant violation. They match the
categories applications already handle, so they translate cleanly into any
host's error model.

There are no exceptions, no `setjmp`, and no conversion of hardware signals into
throwable objects.

An algorithm that "did not converge" returns that fact. It does not return a
null shape and set a flag that the caller may forget to check.

> *Elsewhere:* a thrown `Failure` hierarchy and, in at least one kernel, a
> facility that turns SIGSEGV into a catchable exception. We do neither.

---

## Checklist for a new algorithm

1. Written against the geometry traits (§10), not concrete types.
2. Populates `generated` / `modified` / `is_deleted` (§7).
3. Assigns provenance to every entity it creates (§8).
4. Composes location and orientation correctly during traversal (§1, §2, §3).
5. Uses the right identity predicate, with a matching hasher (§4).
6. Keeps the tolerance containment rule, or fails loudly (§5).
7. Keeps edge representations consistent, or clears `same_parameter` (§6).
8. Returns `Result`; never a silently invalid shape (§12).
