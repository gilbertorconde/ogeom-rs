- **Added: blend into round for more boss feet** (#1171).
  - A boss's foot round now blends into every round its contact circle
    reaches, over arcs that stand apart (for example at a corner where two
    convex edge rounds of the floor meet), into a concave round tangent to
    the floor, and from a conical wall narrowing away from the floor. Each
    is fitted through exact sections and meets the torus G1 (ADR 0001,
    amendment "overflow handling", "Blend into round, more feet").
  - These feet used to be refused by trim through. With blend into round
    off they still are; option order and option-off behaviour are
    unchanged.
  - The proof that nothing else stands in the fill's way is now tighter:
    boxes along the sections in place of one box round the whole foot, and
    every other face near a ball measured from its centre.
