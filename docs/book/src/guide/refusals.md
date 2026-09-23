# Refusals

A plausible wrong answer is worse than an error: a boolean that leaves a
sliver of the tool inside, a blend that looks tangent but is not, an import
that silently drops a face. ogeom's rule: **when the kernel cannot produce a
correct result, it refuses, and the error names the reason.**

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:refused_by_name}}
```

Where this applies:

- **Degenerate inputs** are refused at construction (a zero radius, an
  empty wire, a face whose boundary does not close). The error names the
  offending parameter.
- **Restricted capabilities refuse outside their limits.** Example: the
  medial axis supports convex polygonal faces only. Given a reflex corner,
  a hole or an arc, it returns an error saying which one it found, instead
  of a wrong axis. Every `partial` row in
  [the parity ledger](../kernel/parity-ledger.md) states its restriction,
  and the code refuses at that same limit.
- **Exchange readers report what they skip.** An unsupported entity appears
  in the import report by name and number. Everything that was translated
  can be trusted.
- **Repairs report what they achieved.** Healing operations return measured
  deviations. A repair that cannot reach tolerance says so.

For callers: **treat every `Err` as information.** The message says which
input, which limit and which capability boundary you hit. The parity
ledger tells you whether that boundary is a known, scoped restriction.
