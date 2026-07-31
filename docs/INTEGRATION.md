# Integration

**openGeometry's product is its own API.** Integration layers for other
languages and host applications are downstream, optional, and none of them
constrains the kernel's design.

This file exists for one narrow purpose: to record what embedding a kernel into
a real application demands, so that decisions taken now do not quietly foreclose
it. Nothing here should ever *drive* work. It should only ever veto a design that
would make embedding impossible.

---

## The design pressure worth taking seriously

Surveying a large application that embeds a B-rep kernel (~4,100 call sites)
turned up one finding that generalises beyond it:

**Consumers do not merely call a kernel. They extend it.** That survey found
seven classes deriving from kernel types and overriding their virtuals — a shape
subclass that intercepts every mutating member to keep an element map coherent,
several operation subclasses, custom message and progress sinks.

Two consequences:

1. **A pure C ABI is never sufficient** for a host that wants to specialise
   kernel behaviour. Any C++ integration layer will need real classes with real
   virtual dispatch on top of whatever FFI surface we expose.
2. **Extension points belong in the design, not in the shim.** Progress
   reporting, cancellation, diagnostics, custom tolerance policy and history
   observation should be traits in the Rust API. If a host has to subclass its
   way to them, we designed them wrong.

That second point is the actionable one, and it is a kernel-side task rather than
an integration-side one.

---

## What the kernel must already do

Satisfied by `DATA_MODEL.md`. Listed with the consequence of drifting.

| Requirement | §  | Consequence of getting it wrong |
|---|---|---|
| Shape is `(tshape, location, orientation)` and cheap to copy | §1 | Every by-value shape parameter in every host becomes an allocation |
| Location is a chain, not a flat matrix | §2 | Assembly instancing collapses; placement identity stops being decidable structurally |
| Orientation composes on descent | §3 | Face normals flip inconsistently — silent, and catastrophic downstream |
| `is_same` / `is_equal` / `is_partner` with matching hashers | §4 | Shape maps mis-key. Silent wrong answers |
| Per-entity tolerances with the containment rule | §5 | Imported geometry cannot be modeled with at all |
| Edges carry a representation list including per-face pcurves | §6 | Boolean face splitting has nothing to split with |
| `generated` / `modified` / `is_deleted` on every operation | §7 | **Downstream naming breaks silently and corrupts user documents** |
| Stable provenance | §8 | References into a rebuilt model cannot be resolved at all |
| `Result`, no exceptions, no signal conversion | §12 | Failures cannot be mapped cleanly into a host's error model |

The history row is the one to watch, because it is the only failure mode above
that is quiet. A parametric application records "fillet *that* edge" and resolves
it after a rebuild by walking history. Half-populated history does not error — it
reopens the document with the wrong faces filleted.

---

## Planned integration layers

None of these are scheduled. They are listed so their requirements are visible.

**C ABI** (`og-capi`) — the foundation for everything else. Opaque handles, POD
structs, explicit ownership. Straightforward once the native API is stable.

**Python** — the highest-value binding by a wide margin: it is how most people
would actually try this kernel. PyO3 over the native API, not over the C ABI.

**C++** — real classes with virtuals over the C ABI, for hosts that want to
specialise behaviour.

**Drop-in replacement for an existing kernel's headers** — technically possible:
a source-compatible façade exposing another kernel's class names and signatures,
built into libraries with the names that kernel's build-system probes expect, so
a consumer recompiles without being edited. Feasible, large, and firmly a
downstream project. Two things make it a poor thing to design *toward*:

- it drags in the other kernel's mistakes wholesale, including the pointer
  identity model that `DATA_MODEL.md` §8 exists to escape;
- some things cannot be supported at all. Any host API that hands a raw shape
  pointer to a third-party binding runtime requires binary layout compatibility,
  which is not a goal and would poison the design if it were.

**WASM** — the kernel is pure Rust with no C dependencies, so this is close to
free and worth keeping that way. Weigh it before adding any dependency that
would compromise it.

---

## Rules that protect this without constraining it

1. The native Rust API is designed for Rust. No parameter exists because a
   binding might want it.
2. Extension points are traits in the native API — progress, cancellation,
   diagnostics, tolerance policy, history observation.
3. No public type's design is compromised for C representability. The C ABI
   deals in handles; that is its job.
4. No dependency that would break WASM or introduce a C toolchain requirement,
   without an explicit decision recorded here.
