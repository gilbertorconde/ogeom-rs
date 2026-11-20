- **Fixed: asymmetric rounds on a convex drum ruling past where their
  section fits** (#907).
  - As the drum's leg turns far past the flat's, the section's conic
    bulges toward the drum and then crosses it between the contacts. That
    round is now refused by name ("crosses the drum"), read exactly. Before,
    it built a solid bulging out of the original, valid by `check`.
  - Just short of that scale the build no longer fails in the boolean: a
    cut's prism runs on past an end whose plane bounds the whole solid,
    instead of sharing that end plane with the solid.
  - `largest_fillet_distances` on a ruling stops there too. On the half
    drum with legs 0.1 and 1 it answers 1.6914 (the section meets the drum
    at 1.69143), where it answered 1.9085, already past that limit.
