# Refusals

A CAD kernel's worst failure mode is not the error — it is the plausible
wrong answer: the boolean that returns a shape with a sliver of the tool
left inside, the blend that looks tangent until the toolpath gouges, the
import that silently dropped a face. ogeom's standing rule is that **when
the kernel cannot do something honestly, it refuses, and the refusal
names the reason**.

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:refused_by_name}}
```

This shows up everywhere:

- **Degenerate inputs** are refused at construction — a zero radius, an
  empty wire, a face whose boundary does not close — with the offending
  parameter named.
- **Restricted capabilities refuse outside their restriction.** The
  medial axis, for instance, is built for convex polygonal faces exactly;
  hand it a reflex corner, a hole or an arc and the error says which of
  those it met and why it matters, rather than returning an axis that is
  quietly wrong. Every `partial` row in
  [the parity ledger](../kernel/parity-ledger.md) states its restriction,
  and the code refuses at that same boundary.
- **Exchange readers report, not swallow.** An entity outside the
  supported set appears in the import report by name and number; the
  geometry that was translated is trustworthy precisely because the
  reader does not pretend about the rest.
- **Repairs report what they achieved.** Healing operations return
  measured deviations, and a repair that cannot reach tolerance says so.

For a caller, the practical consequence: **treat every `Err` as
information, not noise.** The message is written to tell you which input,
which limit, and which capability boundary you hit — and the parity
ledger tells you whether that boundary is a restriction someone has
already scoped.
