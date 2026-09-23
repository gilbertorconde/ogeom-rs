# Integration

**ogeom's product is its own API.** Integration layers for other languages and
host applications are downstream and optional. None of them constrains the
kernel's design.

This file has one purpose: to record what embedding a kernel in a real
application requires, so that today's decisions do not rule it out by accident.
Nothing here should *drive* work. It can only veto a design that would make
embedding impossible.

---

## The design pressure that matters

We surveyed a large application that embeds a B-rep kernel (about 4,100 call
sites). One finding applies beyond that application:

**Consumers do not just call a kernel. They extend it.** The survey found seven
classes that derive from kernel types and override their virtual methods:

- a shape subclass that intercepts every mutating member to keep an element map
  consistent;
- several operation subclasses;
- custom message and progress sinks.

Two consequences:

1. **A pure C ABI is never enough** for a host that wants to specialise kernel
   behaviour. A C++ integration layer will need real classes with real virtual
   dispatch on top of the FFI surface.
2. **Extension points belong in the kernel's design, not in the shim.** Progress
   reporting, cancellation, diagnostics, custom tolerance policy and history
   observation should be traits in the Rust API. If a host has to subclass to
   reach them, the design is wrong.

The second point is the actionable one. It is kernel work, not integration work.

---

## What the kernel must already do

`DATA_MODEL.md` covers all of these. Each row gives the consequence of getting
it wrong.

| Requirement | §  | Consequence of getting it wrong |
|---|---|---|
| Shape is `(tshape, location, orientation)` and cheap to copy | §1 | Every by-value shape parameter in every host becomes an allocation |
| Location is a chain, not a flat matrix | §2 | Assembly instancing breaks; placement identity can no longer be decided structurally |
| Orientation composes on descent | §3 | Face normals flip inconsistently. This is silent, and catastrophic downstream |
| `is_same` / `is_equal` / `is_partner` with matching hashers | §4 | Shape maps use the wrong keys. Silent wrong answers |
| Per-entity tolerances with the containment rule | §5 | Imported geometry cannot be modelled with at all |
| Edges carry a representation list including per-face pcurves | §6 | Boolean face splitting has nothing to split with |
| `generated` / `modified` / `is_deleted` on every operation | §7 | **Downstream naming breaks silently and corrupts user documents** |
| Stable provenance | §8 | References into a rebuilt model cannot be resolved at all |
| `Result`, no exceptions, no signal conversion | §12 | Failures cannot be mapped cleanly into a host's error model |

Watch the history row most closely: it is the only failure above that is
silent. A parametric application records "fillet *that* edge" and, after a
rebuild, finds the edge again by walking history. Incomplete history does not
raise an error. It reopens the document with the wrong faces filleted.

---

## Planned integration layers

None of these is scheduled. They are listed so their requirements stay visible.

| Layer | Notes |
|---|---|
| **C ABI** (`og-capi`) | The base for everything else. Opaque handles, POD structs, explicit ownership. Straightforward once the native API is stable. |
| **Python** | By far the most valuable binding: it is how most people would try the kernel. Built with PyO3 over the native API, not over the C ABI. |
| **C++** | Real classes with virtual methods over the C ABI, for hosts that want to specialise behaviour. |
| **WASM** | The kernel is pure Rust with no C dependencies, so this is nearly free. Keep it that way: weigh WASM before adding any dependency that could break it. |

**Drop-in replacement for another kernel's headers.** This is technically
possible: a source-compatible façade that exposes another kernel's class names
and signatures, built into libraries with the names that kernel's build-system
probes expect, so a consumer recompiles without code changes. It is feasible,
large, and firmly a downstream project. It is a poor thing to design *toward*,
for two reasons:

- it imports the other kernel's mistakes wholesale, including the
  pointer-identity model that `DATA_MODEL.md` §8 exists to avoid;
- some things cannot be supported at all. A host API that passes a raw shape
  pointer to a third-party binding runtime needs binary layout compatibility.
  That is not a goal, and pursuing it would damage the design.

---

## Rules that protect integration without constraining the kernel

1. The native Rust API is designed for Rust. No parameter exists just because a
   binding might want it.
2. Extension points are traits in the native API: progress, cancellation,
   diagnostics, tolerance policy, history observation.
3. No public type is compromised to make it representable in C. The C ABI deals
   in handles; that is its job.
4. No dependency that would break WASM or require a C toolchain, unless an
   explicit decision is recorded here.
