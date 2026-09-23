# Documents, assemblies and PMI

A `Model` holds geometry. An `ogeom::doc::Document` holds the product
information around it: parts, the assembly tree, names, colours and
annotations. The exchange formats read and write documents, because that
is what a STEP file contains.

```rust,ignore
let mut document = ogeom::doc::Document::over(model);
let bolt = document.add_part("bolt", bolt_shape);
```

## Assemblies

Products form a tree:

- A part is a leaf.
- An assembly places other products as `Instance`s.

Instances share geometry. Two bolts in an assembly point to one shape node
under different location chains, so a thousand fasteners are not a
thousand copies.

- `roots()` returns the top-level products.
- `occurrences_of(root)` flattens the tree into placed `Occurrence`s. Each
  has a path string and its placed shape.

## Attributes and PMI

Colours and named attributes attach to products and to individual faces.

PMI (the dimensions, geometric tolerances and datums of a manufacturing
drawing) is stored in two forms:

- **Semantic:** a dimension knows which faces it measures, through the same
  stable references used everywhere else.
- **Presentation:** `Callout` polylines, the drawn form.

Both forms, and the distinction between them, survive STEP.

## Views and notes

- A `View` is a named camera plus the callouts it shows. Annotated models
  use views to organise PMI into readable sheets.
- A `Note` is authored text, optionally attached to a product.

Both survive the native format and STEP.

## Undo

The document records every structural change as a step. Undo and redo walk
those steps. `undo_depth()` returns how many steps are available in each
direction.
