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

impl ProductId {
    /// The id's position in the document's own product order — the index a
    /// file format writes and rebinds by re-adding in order.
    #[must_use]
    pub const fn index(self) -> u32 {
        self.0
    }
}

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
    pmi: crate::pmi::Pmi,
    properties: HashMap<TShapeId, Vec<crate::attributes::Property>>,
    materials: Vec<crate::attributes::Material>,
    material_of: HashMap<TShapeId, crate::attributes::MaterialId>,
    layers: Vec<crate::attributes::Layer>,
    on_layer: HashMap<TShapeId, Vec<crate::attributes::LayerId>>,
    validation: HashMap<TShapeId, crate::attributes::ValidationProperties>,
    textures: Vec<crate::attributes::Texture>,
    texture_of: HashMap<TShapeId, crate::attributes::TextureId>,
    views: Vec<crate::view::View>,
    notes: Vec<crate::view::Note>,
    /// Document states an undo can return to, oldest first, and how far
    /// back through them the caller currently stands.
    history: Vec<State>,
    /// How many of `history`'s tail have been undone: the redo depth.
    undone: usize,
}

/// Everything a document holds *about* a model, which is everything an undo
/// can restore.
///
/// The model itself is not in here, and that is the design rather than an
/// omission. Geometry arenas are append-only — a boolean's result does not
/// erase its inputs, it stands beside them — so undoing an operation means
/// putting back what the document *said*, not unmaking what the model
/// holds. The nodes the undone operation built stay where they are,
/// unreferenced, which is what a garbage-collected arena is for.
#[derive(Debug, Clone, Default)]
struct State {
    products: Vec<Product>,
    colours: HashMap<TShapeId, Colour>,
    names: HashMap<TShapeId, String>,
    pmi: crate::pmi::Pmi,
    properties: HashMap<TShapeId, Vec<crate::attributes::Property>>,
    materials: Vec<crate::attributes::Material>,
    material_of: HashMap<TShapeId, crate::attributes::MaterialId>,
    layers: Vec<crate::attributes::Layer>,
    on_layer: HashMap<TShapeId, Vec<crate::attributes::LayerId>>,
    validation: HashMap<TShapeId, crate::attributes::ValidationProperties>,
    textures: Vec<crate::attributes::Texture>,
    texture_of: HashMap<TShapeId, crate::attributes::TextureId>,
    views: Vec<crate::view::View>,
    notes: Vec<crate::view::Note>,
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
            ..Self::default()
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

    /// Place `product` inside `assembly` at an already-resolved location.
    ///
    /// The persistence path: a file carries the instance's location chain
    /// verbatim, and re-minting a datum for it would renumber what the file
    /// preserved. Checks are as [`Document::add_instance`].
    ///
    /// # Errors
    ///
    /// As [`Document::add_instance`].
    pub fn add_instance_at(
        &mut self,
        assembly: ProductId,
        product: ProductId,
        location: Location,
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

    /// Add a texture the document can lay on shapes.
    pub fn add_texture(
        &mut self,
        texture: crate::attributes::Texture,
    ) -> crate::attributes::TextureId {
        self.textures.push(texture);
        crate::attributes::TextureId(self.textures.len() - 1)
    }

    /// Lay a texture on a shape, replacing whatever was on it.
    pub fn set_texture(&mut self, shape: &Shape, texture: crate::attributes::TextureId) {
        self.texture_of.insert(shape.node(), texture);
    }

    /// The texture on a shape, if it has one.
    #[must_use]
    pub fn texture_of(&self, shape: &Shape) -> Option<&crate::attributes::Texture> {
        let id = self.texture_of.get(&shape.node())?;
        self.textures.get(id.0)
    }

    /// Every texture the document holds, in the order they were added.
    #[must_use]
    pub fn textures(&self) -> &[crate::attributes::Texture] {
        &self.textures
    }

    /// Mark the document's current state as one an undo can return to.
    ///
    /// Call it *before* the change a caller might want back. Anything
    /// undone and not redone is dropped at the next checkpoint, which is
    /// the usual rule: a new branch replaces the abandoned one.
    pub fn checkpoint(&mut self) {
        if self.undone > 0 {
            let keep = self.history.len() - self.undone;
            self.history.truncate(keep);
            self.undone = 0;
        }
        let state = self.state();
        self.history.push(state);
    }

    /// Step back to the last checkpoint. `false` when there is none.
    pub fn undo(&mut self) -> bool {
        if self.history.len() <= self.undone {
            return false;
        }
        // The state being left is kept in its place, so redo has somewhere
        // to go.
        let index = self.history.len() - self.undone - 1;
        let current = self.state();
        let restored = core::mem::replace(&mut self.history[index], current);
        self.apply(restored);
        self.undone += 1;
        true
    }

    /// Step forward again. `false` when nothing was undone.
    pub fn redo(&mut self) -> bool {
        if self.undone == 0 {
            return false;
        }
        let index = self.history.len() - self.undone;
        let current = self.state();
        let restored = core::mem::replace(&mut self.history[index], current);
        self.apply(restored);
        self.undone -= 1;
        true
    }

    /// How many steps back are available, and how many forward.
    #[must_use]
    pub fn undo_depth(&self) -> (usize, usize) {
        (self.history.len() - self.undone, self.undone)
    }

    /// Everything the document says about its model, copied out.

    /// The position of a product in write order — how the native format
    /// refers to one across a save.
    #[must_use]
    pub fn product_index(&self, id: ProductId) -> usize {
        id.index() as usize
    }

    /// Add a saved view; its index is how STEP and the native format refer
    /// to it.
    pub fn add_view(&mut self, view: crate::view::View) -> usize {
        self.views.push(view);
        self.views.len() - 1
    }

    /// The saved views, in order.
    #[must_use]
    pub fn views(&self) -> &[crate::view::View] {
        &self.views
    }

    /// Add a note.
    pub fn add_note(&mut self, note: crate::view::Note) -> usize {
        self.notes.push(note);
        self.notes.len() - 1
    }

    /// The notes, in order.
    #[must_use]
    pub fn notes(&self) -> &[crate::view::Note] {
        &self.notes
    }

    fn state(&self) -> State {
        State {
            products: self.products.clone(),
            colours: self.colours.clone(),
            names: self.names.clone(),
            pmi: self.pmi.clone(),
            properties: self.properties.clone(),
            materials: self.materials.clone(),
            material_of: self.material_of.clone(),
            layers: self.layers.clone(),
            on_layer: self.on_layer.clone(),
            validation: self.validation.clone(),
            textures: self.textures.clone(),
            texture_of: self.texture_of.clone(),
            views: self.views.clone(),
            notes: self.notes.clone(),
        }
    }

    /// Put a state back.
    fn apply(&mut self, state: State) {
        self.products = state.products;
        self.colours = state.colours;
        self.names = state.names;
        self.pmi = state.pmi;
        self.properties = state.properties;
        self.materials = state.materials;
        self.material_of = state.material_of;
        self.layers = state.layers;
        self.on_layer = state.on_layer;
        self.validation = state.validation;
        self.textures = state.textures;
        self.texture_of = state.texture_of;
        self.views = state.views;
        self.notes = state.notes;
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

    /// Replace a part's shape — the modification step of an edit.
    ///
    /// The old shape's node-attached colours, names and PMI stay where they
    /// are: entities that survived the modification keep their annotations,
    /// entities that did not simply no longer resolve.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// `product` is not a part;
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if it is not
    /// in this document.
    pub fn replace_part_shape(&mut self, product: ProductId, shape: Shape) -> OgeomResult<()> {
        let Some(entry) = self.products.get_mut(product.0 as usize) else {
            ogeom_bail!(Dangling, "the product is not in this document");
        };
        let ProductKind::Part { shape: slot } = &mut entry.kind else {
            ogeom_bail!(Construction, "only a part carries a shape to replace");
        };
        *slot = shape;
        Ok(())
    }

    /// The document's semantic PMI.
    #[must_use]
    pub const fn pmi(&self) -> &crate::pmi::Pmi {
        &self.pmi
    }

    /// The PMI, for annotating.
    pub const fn pmi_mut(&mut self) -> &mut crate::pmi::Pmi {
        &mut self.pmi
    }

    // --- attributes ---------------------------------------------------------

    /// Pin a user-defined property to a shape. Properties accumulate; a
    /// repeated name replaces the earlier value.
    pub fn set_property(&mut self, shape: &Shape, property: crate::attributes::Property) {
        let list = self.properties.entry(shape.node()).or_default();
        if let Some(held) = list.iter_mut().find(|p| p.name == property.name) {
            *held = property;
        } else {
            list.push(property);
        }
    }

    /// The properties pinned to a shape.
    #[must_use]
    pub fn properties_of(&self, shape: &Shape) -> &[crate::attributes::Property] {
        self.properties
            .get(&shape.node())
            .map_or(&[], Vec::as_slice)
    }

    /// Every shape with properties, for persistence.
    pub fn properties(&self) -> impl Iterator<Item = (TShapeId, &[crate::attributes::Property])> {
        self.properties.iter().map(|(k, v)| (*k, v.as_slice()))
    }

    /// Add a material to the document's list.
    pub fn add_material(
        &mut self,
        material: crate::attributes::Material,
    ) -> crate::attributes::MaterialId {
        self.materials.push(material);
        crate::attributes::MaterialId(self.materials.len() - 1)
    }

    /// A material by id.
    #[must_use]
    pub fn material(
        &self,
        id: crate::attributes::MaterialId,
    ) -> Option<&crate::attributes::Material> {
        self.materials.get(id.0)
    }

    /// The materials, in id order.
    #[must_use]
    pub fn materials(&self) -> &[crate::attributes::Material] {
        &self.materials
    }

    /// The id at a list position, for rebinding persisted references.
    #[must_use]
    pub fn material_id(&self, index: usize) -> Option<crate::attributes::MaterialId> {
        (index < self.materials.len()).then_some(crate::attributes::MaterialId(index))
    }

    /// Assign a shape its material.
    pub fn assign_material(&mut self, shape: &Shape, id: crate::attributes::MaterialId) {
        self.material_of.insert(shape.node(), id);
    }

    /// The material a shape is assigned, if any.
    #[must_use]
    pub fn material_of(&self, shape: &Shape) -> Option<crate::attributes::MaterialId> {
        self.material_of.get(&shape.node()).copied()
    }

    /// Every material assignment, for persistence.
    pub fn material_assignments(
        &self,
    ) -> impl Iterator<Item = (TShapeId, crate::attributes::MaterialId)> + '_ {
        self.material_of.iter().map(|(k, v)| (*k, *v))
    }

    /// Add a layer, visible by default.
    pub fn add_layer(&mut self, name: impl Into<String>) -> crate::attributes::LayerId {
        self.layers.push(crate::attributes::Layer {
            name: name.into(),
            visible: true,
        });
        crate::attributes::LayerId(self.layers.len() - 1)
    }

    /// A layer by id.
    #[must_use]
    pub fn layer(&self, id: crate::attributes::LayerId) -> Option<&crate::attributes::Layer> {
        self.layers.get(id.0)
    }

    /// The layers, in id order.
    #[must_use]
    pub fn layers(&self) -> &[crate::attributes::Layer] {
        &self.layers
    }

    /// The id at a list position, for rebinding persisted references.
    #[must_use]
    pub fn layer_id(&self, index: usize) -> Option<crate::attributes::LayerId> {
        (index < self.layers.len()).then_some(crate::attributes::LayerId(index))
    }

    /// Show or hide a layer.
    pub fn set_layer_visible(&mut self, id: crate::attributes::LayerId, visible: bool) {
        if let Some(layer) = self.layers.get_mut(id.0) {
            layer.visible = visible;
        }
    }

    /// Put a shape on a layer. A shape may sit on several.
    pub fn place_on_layer(&mut self, shape: &Shape, layer: crate::attributes::LayerId) {
        let list = self.on_layer.entry(shape.node()).or_default();
        if !list.contains(&layer) {
            list.push(layer);
        }
    }

    /// The layers a shape sits on.
    #[must_use]
    pub fn layers_of(&self, shape: &Shape) -> &[crate::attributes::LayerId] {
        self.on_layer.get(&shape.node()).map_or(&[], Vec::as_slice)
    }

    /// Every layer membership, for persistence.
    pub fn layer_memberships(
        &self,
    ) -> impl Iterator<Item = (TShapeId, &[crate::attributes::LayerId])> {
        self.on_layer.iter().map(|(k, v)| (*k, v.as_slice()))
    }

    /// Record validation values for a shape.
    pub fn set_validation(
        &mut self,
        shape: &Shape,
        values: crate::attributes::ValidationProperties,
    ) {
        self.validation.insert(shape.node(), values);
    }

    /// The recorded validation values for a shape.
    #[must_use]
    pub fn validation_of(&self, shape: &Shape) -> Option<crate::attributes::ValidationProperties> {
        self.validation.get(&shape.node()).copied()
    }

    /// Every validation record, for persistence.
    pub fn validations(
        &self,
    ) -> impl Iterator<Item = (TShapeId, crate::attributes::ValidationProperties)> + '_ {
        self.validation.iter().map(|(k, v)| (*k, *v))
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
