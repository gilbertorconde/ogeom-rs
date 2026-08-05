//! Product structure: parts, assemblies, instances with placements.
//!
//! A [`Document`] owns a [`Model`] and a table of products over it. A *part*
//! is a product with a shape; an *assembly* is a product with instances, each
//! naming a product and carrying a placement. Instancing leans directly on
//! the location chain (`docs/DATA_MODEL.md` §2): every instance's placement
//! is a datum in the model's own store, an occurrence of a part is the part's
//! shape *moved* by the chain of placements above it, and ten thousand
//! identical fasteners are ten thousand chains over one node.
//!
//! Appearance and naming ride alongside: a colour or a name attaches to a
//! product or to a topology node — a whole part or one face of it — and
//! resolution walks from the most specific to the least.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::Transform;
use ogeom_topo::{Location, Model, Shape, TShapeId};
use std::collections::HashMap;

/// A product in a document: a part or an assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ProductId(u32);

/// An RGBA colour, each channel in `[0, 1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Colour {
    /// Red.
    pub r: f64,
    /// Green.
    pub g: f64,
    /// Blue.
    pub b: f64,
    /// Opacity: 1 is opaque.
    pub a: f64,
}

impl Colour {
    /// An opaque colour from red, green and blue.
    #[must_use]
    pub const fn rgb(r: f64, g: f64, b: f64) -> Self {
        Self { r, g, b, a: 1.0 }
    }
}

/// One use of a product inside an assembly: which product, where, and what
/// this particular occurrence is called.
#[derive(Debug, Clone)]
pub struct Instance {
    /// The product this instance places.
    pub product: ProductId,
    /// The placement, as a location over the document's own datum store.
    pub location: Location,
    /// The occurrence's own name — "bolt-3", not the product's "bolt".
    pub name: Option<String>,
}

/// What a product is: geometry, or uses of other products.
#[derive(Debug, Clone)]
pub enum ProductKind {
    /// A part: a product that carries a shape.
    Part {
        /// The part's geometry, a shape in the document's model.
        shape: Shape,
    },
    /// An assembly: a product made of placed uses of other products.
    Assembly {
        /// The assembly's instances, in authoring order.
        children: Vec<Instance>,
    },
}

/// One product: a name, an optional colour, and what it is.
#[derive(Debug, Clone)]
pub struct Product {
    /// The product's name.
    pub name: String,
    /// The product's own colour, inherited by anything in it that has none.
    pub colour: Option<Colour>,
    /// Part or assembly.
    pub kind: ProductKind,
}

/// A placed part: the flattening of an assembly tree into shapes.
#[derive(Debug, Clone)]
pub struct Occurrence {
    /// The part this occurrence places.
    pub part: ProductId,
    /// The part's shape, moved by every placement above it.
    pub shape: Shape,
    /// The path of names from the root to here, `/`-separated, instance
    /// names where they exist and product names where they do not.
    pub path: String,
}

/// A model with product structure, appearance and names over it.
///
/// The model stays reachable — construction, booleans and measurement all
/// operate on it directly — and the document adds what a model alone does
/// not say: which shapes are products, how they assemble, what they are
/// called and what colour they are.
#[derive(Debug, Default)]
pub struct Document {
    model: Model,
    products: Vec<Product>,
    colours: HashMap<TShapeId, Colour>,
    names: HashMap<TShapeId, String>,
}

impl Document {
    /// An empty document over an empty model.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// A document over an existing model.
    #[must_use]
    pub fn over(model: Model) -> Self {
        Self {
            model,
            products: Vec::new(),
            colours: HashMap::new(),
            names: HashMap::new(),
        }
    }

    /// The model under the document.
    #[must_use]
    pub const fn model(&self) -> &Model {
        &self.model
    }

    /// The model, for construction and modification.
    pub const fn model_mut(&mut self) -> &mut Model {
        &mut self.model
    }

    /// Add a part: a named product carrying a shape.
    pub fn add_part(&mut self, name: impl Into<String>, shape: Shape) -> ProductId {
        self.push(Product {
            name: name.into(),
            colour: None,
            kind: ProductKind::Part { shape },
        })
    }

    /// Add an empty assembly.
    pub fn add_assembly(&mut self, name: impl Into<String>) -> ProductId {
        self.push(Product {
            name: name.into(),
            colour: None,
            kind: ProductKind::Assembly {
                children: Vec::new(),
            },
        })
    }

    fn push(&mut self, product: Product) -> ProductId {
        self.products.push(product);
        #[allow(clippy::cast_possible_truncation)]
        ProductId(self.products.len() as u32 - 1)
    }

    /// Place `product` inside `assembly` at `at`.
    ///
    /// The transform becomes a datum in the model's own store, so the
    /// instance's placement is structural — comparable by identity, shared by
    /// every traversal — rather than a matrix to be compared with an epsilon.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// `assembly` is not an assembly, or placing `product` there would make a
    /// product contain itself;
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if either id
    /// is not in this document.
    pub fn add_instance(
        &mut self,
        assembly: ProductId,
        product: ProductId,
        at: Transform,
        name: Option<String>,
    ) -> OgeomResult<()> {
        if self.get(product).is_none() {
            ogeom_bail!(Dangling, "the product to place is not in this document");
        }
        if self.contains_product(product, assembly) {
            ogeom_bail!(
                Construction,
                "placing this product here would make it contain itself"
            );
        }
        let location = if at == Transform::IDENTITY {
            Location::identity()
        } else {
            Location::of(self.model.add_datum(at))
        };
        let Some(entry) = self.products.get_mut(assembly.0 as usize) else {
            ogeom_bail!(Dangling, "the assembly is not in this document");
        };
        let ProductKind::Assembly { children } = &mut entry.kind else {
            ogeom_bail!(Construction, "instances go inside assemblies, not parts");
        };
        children.push(Instance {
            product,
            location,
            name,
        });
        Ok(())
    }

    /// Whether the tree under `haystack` reaches `needle`.
    fn contains_product(&self, haystack: ProductId, needle: ProductId) -> bool {
        if haystack == needle {
            return true;
        }
        match self.get(haystack).map(|p| &p.kind) {
            Some(ProductKind::Assembly { children }) => children
                .iter()
                .any(|i| self.contains_product(i.product, needle)),
            _ => false,
        }
    }

    /// The product behind an id.
    #[must_use]
    pub fn get(&self, id: ProductId) -> Option<&Product> {
        self.products.get(id.0 as usize)
    }

    /// Every product, in the order added.
    pub fn products(&self) -> impl Iterator<Item = (ProductId, &Product)> {
        self.products
            .iter()
            .enumerate()
            .map(|(i, p)| (ProductId(u32::try_from(i).unwrap_or(u32::MAX)), p))
    }

    /// The products no instance places: the top of the tree.
    #[must_use]
    pub fn roots(&self) -> Vec<ProductId> {
        let mut placed = vec![false; self.products.len()];
        for product in &self.products {
            if let ProductKind::Assembly { children } = &product.kind {
                for instance in children {
                    placed[instance.product.0 as usize] = true;
                }
            }
        }
        placed
            .iter()
            .enumerate()
            .filter(|&(_, &used)| !used)
            .map(|(i, _)| ProductId(u32::try_from(i).unwrap_or(u32::MAX)))
            .collect()
    }

    /// Every placed part under `product`, shapes moved into world space.
    ///
    /// The flattening every consumer of an assembly wants: two instances of
    /// one part come back as two occurrences whose shapes share a topology
    /// node and differ only in their location chains.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if `product`
    /// is not in this document.
    pub fn occurrences_of(&self, product: ProductId) -> OgeomResult<Vec<Occurrence>> {
        let Some(root) = self.get(product) else {
            ogeom_bail!(Dangling, "the product is not in this document");
        };
        let mut out = Vec::new();
        self.flatten(
            product,
            root,
            &Location::identity(),
            &root.name.clone(),
            &mut out,
        );
        Ok(out)
    }

    fn flatten(
        &self,
        id: ProductId,
        product: &Product,
        above: &Location,
        path: &str,
        out: &mut Vec<Occurrence>,
    ) {
        match &product.kind {
            ProductKind::Part { shape } => out.push(Occurrence {
                part: id,
                shape: shape.moved(above),
                path: path.to_string(),
            }),
            ProductKind::Assembly { children } => {
                for instance in children {
                    let Some(child) = self.get(instance.product) else {
                        continue;
                    };
                    let below = above.then(&instance.location);
                    let step = instance.name.as_deref().unwrap_or(&child.name);
                    let path = format!("{path}/{step}");
                    self.flatten(instance.product, child, &below, &path, out);
                }
            }
        }
    }

    /// Colour a product: the fallback for everything in it.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if `product`
    /// is not in this document.
    pub fn set_product_colour(&mut self, product: ProductId, colour: Colour) -> OgeomResult<()> {
        let Some(entry) = self.products.get_mut(product.0 as usize) else {
            ogeom_bail!(Dangling, "the product is not in this document");
        };
        entry.colour = Some(colour);
        Ok(())
    }

    /// Colour a shape — a whole part's shape or one sub-shape of it.
    ///
    /// Keyed by the topology node, so every occurrence of an instanced shape
    /// shows the colour: the colour belongs to the entity, not to one
    /// placement of it.
    pub fn set_colour(&mut self, shape: &Shape, colour: Colour) {
        self.colours.insert(shape.node(), colour);
    }

    /// The colour set directly on a shape's node, if any.
    #[must_use]
    pub fn colour_of(&self, shape: &Shape) -> Option<Colour> {
        self.colours.get(&shape.node()).copied()
    }

    /// The colour a sub-shape of a part actually shows.
    ///
    /// Most specific wins: the sub-shape's own colour, else the colour of the
    /// part's whole shape, else the product's, else nothing.
    #[must_use]
    pub fn resolved_colour(&self, part: ProductId, sub: &Shape) -> Option<Colour> {
        if let Some(own) = self.colour_of(sub) {
            return Some(own);
        }
        let product = self.get(part)?;
        if let ProductKind::Part { shape } = &product.kind
            && let Some(whole) = self.colour_of(shape)
        {
            return Some(whole);
        }
        product.colour
    }

    /// Name a shape's node — a face someone will want to find again.
    pub fn set_name(&mut self, shape: &Shape, name: impl Into<String>) {
        self.names.insert(shape.node(), name.into());
    }

    /// The name set on a shape's node, if any.
    #[must_use]
    pub fn name_of(&self, shape: &Shape) -> Option<&str> {
        self.names.get(&shape.node()).map(String::as_str)
    }

    /// Every node-attached colour, for a writer to carry out.
    pub fn colours(&self) -> impl Iterator<Item = (TShapeId, Colour)> + '_ {
        self.colours.iter().map(|(&node, &colour)| (node, colour))
    }

    /// Every node-attached name, for a writer to carry out.
    pub fn names(&self) -> impl Iterator<Item = (TShapeId, &str)> {
        self.names.iter().map(|(&node, name)| (node, name.as_str()))
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_core::Tolerances;
    use ogeom_math::{Frame, Point, Vector};
    use ogeom_topo::{Filter, ShapeType, explore};

    const T: Tolerances = Tolerances::millimetres();

    fn box_part(document: &mut Document, name: &str, size: f64) -> (ProductId, Shape) {
        let shape = ogeom_algo::make_box(document.model_mut(), Frame::WORLD, (size, size, size), T)
            .unwrap()
            .shape;
        (document.add_part(name, shape.clone()), shape)
    }

    #[test]
    fn two_instances_of_one_part_share_the_node_and_differ_in_placement() {
        let mut document = Document::new();
        let (bolt, shape) = box_part(&mut document, "bolt", 1.0);
        let assembly = document.add_assembly("plate");
        document
            .add_instance(assembly, bolt, Transform::IDENTITY, Some("bolt-1".into()))
            .unwrap();
        document
            .add_instance(
                assembly,
                bolt,
                Transform::translation(Vector::new(10.0, 0.0, 0.0)),
                Some("bolt-2".into()),
            )
            .unwrap();

        let occurrences = document.occurrences_of(assembly).unwrap();
        assert_eq!(occurrences.len(), 2);
        assert_eq!(occurrences[0].path, "plate/bolt-1");
        assert_eq!(occurrences[1].path, "plate/bolt-2");
        // One node, two placements: the instancing claim itself.
        assert_eq!(occurrences[0].shape.node(), shape.node());
        assert_eq!(occurrences[1].shape.node(), shape.node());
        let at = |o: &Occurrence| {
            o.shape
                .transform(document.model().datums())
                .unwrap()
                .apply(Point::new(0.0, 0.0, 0.0))
        };
        assert!(at(&occurrences[0]).is_equal(Point::new(0.0, 0.0, 0.0), T));
        assert!(at(&occurrences[1]).is_equal(Point::new(10.0, 0.0, 0.0), T));
    }

    #[test]
    fn nested_placements_compose_outer_then_inner() {
        let mut document = Document::new();
        let (part, _) = box_part(&mut document, "washer", 1.0);
        let sub = document.add_assembly("stack");
        let top = document.add_assembly("machine");
        document
            .add_instance(
                sub,
                part,
                Transform::translation(Vector::new(0.0, 5.0, 0.0)),
                None,
            )
            .unwrap();
        document
            .add_instance(
                top,
                sub,
                Transform::translation(Vector::new(100.0, 0.0, 0.0)),
                None,
            )
            .unwrap();

        let occurrences = document.occurrences_of(top).unwrap();
        assert_eq!(occurrences.len(), 1);
        assert_eq!(occurrences[0].path, "machine/stack/washer");
        let world = occurrences[0]
            .shape
            .transform(document.model().datums())
            .unwrap()
            .apply(Point::new(0.0, 0.0, 0.0));
        assert!(world.is_equal(Point::new(100.0, 5.0, 0.0), T));
    }

    #[test]
    fn a_product_cannot_contain_itself() {
        let mut document = Document::new();
        let a = document.add_assembly("a");
        let b = document.add_assembly("b");
        document
            .add_instance(a, b, Transform::IDENTITY, None)
            .unwrap();
        let err = document.add_instance(b, a, Transform::IDENTITY, None);
        assert!(err.is_err(), "a cycle must be refused");
        let direct = document.add_instance(a, a, Transform::IDENTITY, None);
        assert!(direct.is_err(), "self-containment must be refused");
    }

    #[test]
    fn instances_go_only_inside_assemblies() {
        let mut document = Document::new();
        let (part, _) = box_part(&mut document, "block", 1.0);
        let (other, _) = box_part(&mut document, "pin", 1.0);
        assert!(
            document
                .add_instance(part, other, Transform::IDENTITY, None)
                .is_err()
        );
    }

    #[test]
    fn roots_are_the_products_nothing_places() {
        let mut document = Document::new();
        let (part, _) = box_part(&mut document, "gear", 1.0);
        let assembly = document.add_assembly("gearbox");
        document
            .add_instance(assembly, part, Transform::IDENTITY, None)
            .unwrap();
        assert_eq!(document.roots(), vec![assembly]);
    }

    #[test]
    fn colour_resolution_walks_sub_shape_then_shape_then_product() {
        let mut document = Document::new();
        let (part, shape) = box_part(&mut document, "housing", 2.0);
        let face = explore(document.model(), &shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .remove(0);

        let red = Colour::rgb(1.0, 0.0, 0.0);
        let green = Colour::rgb(0.0, 1.0, 0.0);
        let blue = Colour::rgb(0.0, 0.0, 1.0);

        assert_eq!(document.resolved_colour(part, &face), None);
        document.set_product_colour(part, blue).unwrap();
        assert_eq!(document.resolved_colour(part, &face), Some(blue));
        document.set_colour(&shape, green);
        assert_eq!(document.resolved_colour(part, &face), Some(green));
        document.set_colour(&face, red);
        assert_eq!(document.resolved_colour(part, &face), Some(red));
        // The whole shape keeps its own colour under the face's override.
        assert_eq!(document.colour_of(&shape), Some(green));
    }

    #[test]
    fn names_attach_to_nodes_and_survive_occurrences() {
        let mut document = Document::new();
        let (_, shape) = box_part(&mut document, "bracket", 1.0);
        let face = explore(document.model(), &shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .remove(0);
        document.set_name(&face, "mounting-face");
        assert_eq!(document.name_of(&face), Some("mounting-face"));
        // The same node reached through a moved occurrence still answers.
        let moved = face.moved(&Location::identity());
        assert_eq!(document.name_of(&moved), Some("mounting-face"));
    }
}
