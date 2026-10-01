# Root cause tracing

A wrong shape shows up far from where it was made: the boolean refuses a
piece that a pcurve put in the wrong place, which a reader fitted from a
record the writer emitted with a reversed seam. The instinct is to fix where
the error appears. That treats a symptom.

Trace backward through the call chain until the original trigger, then fix
at the source.

## When

- The error is deep in a pipeline, not at the entry point.
- The backtrace is long.
- It is unclear where the bad value originated.
- A test or a stress case triggers it and the mechanism is unknown.

## The process

1. **Observe the symptom.** Quote it: the refusal text, the `check`
   diagnosis, the measured volume against the expected one.
2. **Find the immediate cause.** Which code produced that value?
3. **Ask what called it,** and with what. Walk up: the assembly of the
   result, the classification, the section, the surface pair, the reader.
4. **Keep tracing up.** At each step ask what value was passed and whether
   it was already wrong. A tolerance of the wrong scale, a parameter
   outside the domain, a direction not flipped by a mirror, a location
   dropped.
5. **Find the original trigger.** The first place a wrong value exists.
   That is where the fix goes.

## Instrumentation

When the chain cannot be read, print it. In a test or a repro program,
`eprintln!` before the suspect operation, with the inputs, their tolerances,
and the backtrace (`RUST_BACKTRACE=1`, or `std::backtrace::Backtrace::
force_capture()` at the spot). Print before the operation, not after it
fails. Include the parameter values and the entity ids, so two runs can be
compared.

Remove the instrumentation before committing. A diagnostic that earns its
place becomes a named diagnosis in `check` or a refusal with a reason, never
a print.

## Finding which case triggers it

When a property test or the stress harness turns up a failure, isolate it
first: the seed from `*.proptest-regressions`, or `ogeom-stress --case
<scenario>/<part>/<n>`, which replays one case and prints what happened.
Then shrink it by hand: fewer faces, a simpler placement, a round number.
The smallest reproduction is usually one placement of two primitives.

## Then

Having found the source, add validation where the invariant lives, so the
bad value cannot enter again: the tolerance containment check, a domain
assertion on the parameter, a refusal by name at the entry. One guard at
the owning boundary, not a check at every layer that copies it.
