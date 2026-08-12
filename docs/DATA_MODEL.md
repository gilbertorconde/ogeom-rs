# The ogeom data model

**Status: normative.** Everything in `crates/` implements this document. A change
here is a design change and needs to be argued as one — see `CONTRIBUTING.md`.

---

## Why this document exists

A B-rep kernel's public API is a thin skin over its data model. What that model
can express is what the kernel can do: flatten a location chain into a 4×4 matrix
and assembly instancing becomes impossible; give an edge one curve instead of a
list and boolean face splitting has nothing to split with. None of that can be
recovered at an API layer afterwards.

The model below is the one the field converged on over thirty years, because the
alternatives do not work. It is described here in its own terms rather than as a
port of anyone's implementation. Each invariant is:

- **cheap now** — a few days of design attention;
- **effectively impossible to retrofit** — it would mean reworking every algorithm
  in the kernel;
- **load-bearing** — with the failure mode spelled out, so the cost of dropping it
  is concrete rather than theoretical.

Two invariants — stable identity (§8) and predicate abstraction (§9) — are places
where the conventional design is wrong and we deliberately diverge.

Correspondence notes marked *Elsewhere* give the conventional name for a concept,
because that vocabulary is how the field talks about itself. They are a glossary,
not a dependency: ogeom links against no existing CAD kernel. Where a note
cites usage counts, they come from `api_surface.json` — a profile of one large
application, included as evidence that a requirement is real and not as a
specification of scope. See `SCOPE.md`.

---

## 1. A shape is a triple

```rust
pub struct Shape {
    tshape:      TShapeId,      // arena key: the shared, positionless topology node
    location:    Location,      // where this instance sits
    orientation: Orientation,   // which way its boundary faces
}
```

`Shape` is `Copy`-cheap and passed by value everywhere. The heavy data — children,
geometry, tolerances — lives once, in an arena, behind `TShapeId`.

This separation is the whole reason B-rep scales. The same `TShape` appears at many
locations with many orientations without a single byte of geometry being copied.

**Consequence — traversal composes.** A sub-shape's *effective* location is the
product of every location from the root down; its *effective* orientation is the
composition of every orientation on that path. An explorer that yields sub-shapes
without composing both is wrong, and wrong in a way that produces plausible-looking
garbage rather than a crash.

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

A `Datum` is a reference-counted rigid transform. A `Location` is a *sequence* of
`(datum, integer power)` pairs.

**Why not a 4×4 matrix:**

- **Composition is concatenation** — O(1)-ish, no matrix multiply, no drift from
  repeated floating-point composition.
- **Identity comparison is structural.** Two shapes are at "the same place" if
  their chains are equal, decided without comparing 16 floats against a tolerance.
  This is what makes 10,000 identical bolts in an assembly share one piece of
  geometry *and* be recognisable as instances of it.
- **Inverses are exact** — negate the power, rather than inverting a matrix.

The composed `Trsf` is computed lazily and cached. Nothing outside `ogeom-topo`
should ever need to look at the chain itself.

> *Elsewhere:* `TopLoc_Location`, a linked list of `(datum, power)` pairs.

---

## 3. Orientation composes multiplicatively

```rust
pub enum Orientation { Forward, Reversed, Internal, External }
```

- `Forward` — the material is on the surface's default side.
- `Reversed` — the material is on the other side.
- `Internal` — the boundary lies *inside* the material (a stiffener edge embedded
  in a face).
- `External` — the boundary lies outside it (reference geometry).

Composition is a monoid, and it is applied at **every** level of descent:

```
compose(Forward,  x) = x
compose(Reversed, Forward)  = Reversed
compose(Reversed, Reversed) = Forward
compose(Internal, _) = Internal
compose(External, _) = External
```

An edge's orientation *within a face* depends on that face's orientation *within
its shell*, which depends on the shell's orientation within the solid. Reversing a
solid must not require touching a single child.

> *Elsewhere:* `TopAbs_Orientation` and `TopAbs::Compose`.

---

## 4. Identity is a trichotomy

Three distinct equivalences, three distinct hashers, and they are **not**
interchangeable:

| Predicate | Compares | Used for |
|---|---|---|
| `is_partner` | tshape only | "is this the same underlying topology, anywhere, any way round?" |
| `is_same` | tshape + location | set membership in most algorithms; the common case |
| `is_equal` (`==`) | tshape + location + orientation | exact identity; ordered containers |

Every map and set type states which equivalence it uses in its name and enforces it
in its hasher. A `HashMap` keyed on `is_equal` semantics but hashing only the
tshape is a silent correctness bug.

> *Elsewhere:* `IsPartner` / `IsSame` / `IsEqual` and the `ShapeMapHasher`
> family. Conflating them is a recurring, well-documented source of bugs in
> applications built on kernels that expose all three.

---

## 5. Tolerances are per entity, and they grow

Every vertex, edge and face carries its **own** tolerance — the radius of the
sphere/pipe/slab within which the entity is considered to lie.

**Containment rule**, maintained as an invariant:

```
tol(vertex) >= tol(edge) >= tol(face)
```

for entities in a boundary relationship. Operations may only *increase* tolerances,
never silently decrease them, and an operation that cannot satisfy the rule has
failed and must say so.

This is not a workaround for sloppy code. Exact arithmetic cannot represent the
intersection curve of two NURBS surfaces — it is transcendental — so a
tolerance-carrying topology is the only known way to build a kernel that closes.
Every production kernel works this way. See §9 for what exact predicates *can*
do.

> *Elsewhere:* per-entity `Tolerance` on the vertex, edge and face records, a
> validity checker that enforces the rule, and a boolean post-pass that inflates
> tolerances until it holds.

---

## 6. An edge carries a list of representations

Not one curve. A list:

```rust
pub enum EdgeRepr {
    Curve3d      { curve: CurveId, location: Location, range: (f64, f64) },
    PCurve       { curve: Curve2dId, surface: SurfaceId, location: Location },
    PCurveClosed { curve: Curve2dId, curve2: Curve2dId, surface: SurfaceId, location: Location },
    Polygon3d    { polygon: PolygonId, location: Location },
    PolygonOnTri { polygon: PolygonOnTriId, triangulation: TriangulationId, location: Location },
}
```

A single edge simultaneously holds a 3D curve, **one pcurve per adjacent face**, two
pcurves where it is a seam on a closed surface, and cached discretizations.

**Why it must be a list:**

- Face splitting in a boolean happens in 2D parametric space. Without a pcurve on
  each face, there is nothing to split *with*.
- Surfaces are parameterized differently, so one 2D curve cannot serve two faces.
- A seam edge on a cylinder appears at both `u = 0` and `u = 2π`. One pcurve cannot
  express that.

**The `same_parameter` flag** asserts that all representations agree on the
parameterization — i.e. `curve3d(t)` and `surface(pcurve(t))` are the same point
within tolerance, for the same `t`. It is a claim that can be false, and there is a
repair routine that re-establishes it (possibly by growing the edge's tolerance).

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

Not optional, not deferred, not "added when something needs it." Every operation in
`ogeom-algo`, `ogeom-bool`, `ogeom-fillet` and `ogeom-offset` populates these from the first
commit that introduces it.

The reason is downstream: a parametric application records "fillet *that* edge"
and must still find that edge after the model is rebuilt with different
dimensions. It reconstructs the reference by walking history. This is the
topological naming problem, and every application built on a kernel that
identifies topology by pointer has had to solve it this way — in one well-known
case, a decade of work layered on top of the kernel's history maps.

Retrofitting history means revisiting every algorithm, and half-populated history is
worse than none: it fails silently and corrupts documents rather than erroring.

> *Elsewhere:* `Generated` / `Modified` / `IsDeleted` on the operation base class,
> plus a standalone history object.

---

## 8. Entity identity is stable — *deliberate divergence*

The conventional design identifies topology by pointer. Every modeling operation
allocates new nodes, so every reference into a previous result dies. That *is* the
topological naming problem, and every downstream fix is an attempt to reconstruct
identity after the fact by walking history maps.

We record it at creation instead:

```rust
pub struct EntityId(u64);         // stable for the lifetime of a document

pub enum Provenance {
    Primitive { op: OpId, role: PrimitiveRole },   // "the +Z face of box #3"
    Derived   { op: OpId, from: SmallVec<[EntityId; 2]>, role: DerivedRole },
    Imported  { file: FileId, external: ExternalRef },
}
```

An entity's identity is *what produced it and from what*, not where it happens to
live in memory. Splitting a face by a boolean yields fragments that each know they
came from that face; a rebuild with different parameters produces entities with the
same provenance, so a reference like "the fillet on this edge" survives.

History (§7) is still required: it is what an embedding application consumes, and
it is the honest answer where provenance cannot resolve a reference on its own.
But provenance is the primary mechanism, and it only works if it is designed in
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

Algorithms are written against the trait. Implementations may be fast-filtered
floating point, adaptive exact (Shewchuk), or interval-based, and can change
without touching a single algorithm.

**Be clear about what this buys.** Exact predicates solve the *polyhedral*
robustness problem: orientation of a point against a plane, in-sphere for Delaunay.
They do **not** solve the CAD problem, because the intersection curve of two curved
surfaces has no exact representation to be exact *about*. That is why §5 exists and
why it cannot be traded away. Predicates make the parts that *can* be decided
exactly decided exactly; tolerances handle the rest.

Tolerance constants live in `ogeom-core`, and the model's unit scale is **explicit**
rather than assumed to be millimetres — an assumption kernels commonly bake in,
which then misbehaves silently on models authored in metres or inches:

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
against `Curve3d` / `Curve2d` / `Surface`, never against a concrete type. A face
adaptor (surface + location + trimming) and a bare analytic plane are the same thing
to a caller.

`kind()` exists so algorithms can *opt into* analytic fast paths — plane/plane
intersection should not go through a marching intersector — without the general
path ever needing to know what it is looking at.

> *Elsewhere:* the `Adaptor` family. The best idea in the conventional design, and
> we take it wholesale.

---

## 11. Memory: arenas, not reference counting

Topology lives in typed index arenas — `Vec<T>` plus a typed `u32` key — not behind
`Arc` or an intrusive refcount.

- No reference cycles to reason about. Intrusive refcounting has no cycle
  collection, so kernels built on it avoid cycles by convention alone.
- Cache-friendly traversal.
- Keys are small, `Copy`, and serializable.
- It is what makes §8 possible at all.

The cost is that a shape is only meaningful relative to the arena that owns it.
That is the correct trade, and it is made explicit in the API rather than hidden.

The arenas are append-only in practice — nothing in the kernel removes — and two
non-builder paths lean on exactly that. `Model::from_parts` assembles a restored
document by replaying the file's insertion order, which reproduces every handle.
`Model::absorb` is the same engine pointed at a model that already has things in
it: another document's parts append with every handle shifted past what the
target holds, so an absorbed shape is indistinguishable from one built there.
That is how a serialized tool body meets a live one in a boolean. Absorption
preserves the source's identities under a plain offset (the remap table says
where each landed), keeps its provenance verbatim — including source `OpId`s,
which are meaningful only in the source document's own rebuild — and refuses,
by name, a document at another unit scale: rescaling is a feature, not a
default.

> *Elsewhere:* a transient base class with intrusive reference counting plus a
> custom small-block allocator. We need none of it.

---

## 12. Errors are values

`Result<T, OgeomError>`. The variants cover the failure vocabulary a kernel actually
needs — construction, domain, range, dimension mismatch, null object, not done,
numeric failure, invariant violation — chosen to line up with the categories
applications already handle, so they translate cleanly into any host's error
model. There are no exceptions, no `setjmp`, and emphatically no conversion of
hardware signals into throwable objects.

An algorithm that "did not converge" returns that fact. It does not return a null
shape and set a flag for the caller to forget to check.

> *Elsewhere:* a thrown `Failure` hierarchy, and — in at least one kernel — a
> facility that converts SIGSEGV into a catchable exception. We do neither.

---

## Checklist for a new algorithm

1. Written against the geometry traits (§10), not concrete types.
2. Populates `generated` / `modified` / `is_deleted` (§7).
3. Assigns provenance to every entity it creates (§8).
4. Composes location and orientation correctly on traversal (§1, §2, §3).
5. Uses the right identity predicate, with a matching hasher (§4).
6. Maintains the tolerance containment rule, or fails loudly (§5).
7. Keeps edge representations consistent, or clears `same_parameter` (§6).
8. Returns `Result`; never a silently-invalid shape (§12).
