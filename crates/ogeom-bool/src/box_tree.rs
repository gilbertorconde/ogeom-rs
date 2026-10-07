//! Boxes held in a tree, so the boxes meeting one box are found by
//! descending where it reaches rather than by testing every box.

use ogeom_math::Aabb;

/// Below this many boxes a node holds its boxes as they are.
const LEAF: usize = 8;

/// A tree over a set of boxes, each node holding the union of the boxes
/// beneath it.
pub(crate) struct BoxTree {
    nodes: Vec<Node>,
    /// The boxes' indices, each node's a contiguous run.
    order: Vec<usize>,
    boxes: Vec<Aabb>,
}

struct Node {
    bound: Aabb,
    /// The node's run in `order`.
    run: (usize, usize),
    /// The two children's node indices, or none for a leaf.
    children: Option<(usize, usize)>,
}

impl BoxTree {
    /// A tree over `boxes`. An empty box meets nothing and is left out.
    pub(crate) fn new(boxes: &[Aabb]) -> Self {
        let mut order: Vec<usize> = (0..boxes.len()).filter(|&i| !boxes[i].is_empty()).collect();
        let mut tree = Self {
            nodes: Vec::new(),
            order: Vec::new(),
            boxes: boxes.to_vec(),
        };
        if !order.is_empty() {
            let len = order.len();
            tree.build(&mut order, 0, len);
        }
        tree.order = order;
        tree
    }

    /// Build the node over `order[start..end]`, returning its index.
    fn build(&mut self, order: &mut [usize], start: usize, end: usize) -> usize {
        let bound = order[start..end]
            .iter()
            .fold(Aabb::EMPTY, |acc, &i| acc.union(&self.boxes[i]));
        let at = self.nodes.len();
        self.nodes.push(Node {
            bound,
            run: (start, end),
            children: None,
        });
        if end - start <= LEAF {
            return at;
        }
        // Split at the median centre along the box's longest side.
        let size = bound.size();
        let axis = if size.x >= size.y && size.x >= size.z {
            0
        } else if size.y >= size.z {
            1
        } else {
            2
        };
        let centre = |i: usize| -> f64 {
            let c = self.boxes[i].centre().unwrap_or(ogeom_math::Point::ORIGIN);
            [c.x, c.y, c.z][axis]
        };
        let middle = start + (end - start) / 2;
        order[start..end]
            .select_nth_unstable_by(middle - start, |&p, &q| centre(p).total_cmp(&centre(q)));
        let left = self.build(order, start, middle);
        let right = self.build(order, middle, end);
        self.nodes[at].children = Some((left, right));
        at
    }

    /// The indices of the boxes meeting `probe`, ascending, into `out`
    /// (cleared first).
    pub(crate) fn meeting(&self, probe: &Aabb, out: &mut Vec<usize>) {
        out.clear();
        if self.nodes.is_empty() || probe.is_empty() {
            return;
        }
        let mut stack = vec![0_usize];
        while let Some(at) = stack.pop() {
            let node = &self.nodes[at];
            if !node.bound.intersects(probe) {
                continue;
            }
            match node.children {
                Some((left, right)) => {
                    stack.push(left);
                    stack.push(right);
                }
                None => {
                    for &i in &self.order[node.run.0..node.run.1] {
                        if self.boxes[i].intersects(probe) {
                            out.push(i);
                        }
                    }
                }
            }
        }
        out.sort_unstable();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ogeom_math::Point;

    #[test]
    fn finds_what_every_pair_test_finds() {
        // A grid of unit boxes and probes of several sizes across it.
        let mut boxes = Vec::new();
        for i in 0..20 {
            for j in 0..7 {
                let low = Point::new(f64::from(i) * 1.5, f64::from(j) * 2.0, f64::from(i % 3));
                boxes.push(Aabb::of_corners(
                    low,
                    low + ogeom_math::Vector::new(1.0, 1.0, 1.0),
                ));
            }
        }
        boxes.push(Aabb::EMPTY);
        let tree = BoxTree::new(&boxes);
        let mut found = Vec::new();
        for k in 0..40 {
            let low = Point::new(f64::from(k) * 0.8 - 2.0, f64::from(k % 5) * 3.0, 0.5);
            let reach = 0.2 + f64::from(k % 4);
            let probe = Aabb::of_corners(low, low + ogeom_math::Vector::new(reach, reach, reach));
            tree.meeting(&probe, &mut found);
            let every: Vec<usize> = (0..boxes.len())
                .filter(|&i| boxes[i].intersects(&probe))
                .collect();
            assert_eq!(found, every);
        }
    }
}
