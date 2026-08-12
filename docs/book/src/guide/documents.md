# Documents, assemblies and PMI

A `Model` holds geometry; a `ogeom::doc::Document` holds what a product
is made of: parts, the assembly tree, names, colours, annotations. The
exchange formats read and write documents, because that is what a STEP
file actually contains.

```rust,ignore
let mut document = ogeom::doc::Document::over(model);
let bolt = document.add_part("bolt", bolt_shape);
```

## Assemblies

Products form a tree: a part is a leaf, an assembly instances other
products (`Instance`) with placements. Instancing is real — two bolts in
an assembly share one shape node under different location chains, which
is what keeps a thousand-fastener assembly from being a thousand copies.
`roots()` gives the top products; `occurrences_of(root)` flattens the
tree into placed `Occurrence`s, each carrying its path string and its
placed shape.

## Attributes and PMI

Colours and named attributes attach to products and to individual faces.
PMI — the dimensions, geometric tolerances and datums of a manufacturing
drawing — attaches semantically (a dimension knows *which faces* it
measures, via the same stable references everything else uses) and
presentationally (`Callout` polylines, the drawn form). Both survive
STEP, as does the distinction.

## Views and notes

A `View` is a named camera with the callouts it presents — how annotated
models organise their PMI into readable sheets. A `Note` is authored
text, optionally pinned to a product. Both survive the native format and
STEP.

## Undo

The document records its own history: every structural change is a step,
and undo/redo walk them. `undo_depth()` reports how far each direction
goes.
