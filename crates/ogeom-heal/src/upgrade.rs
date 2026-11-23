//! The upgrade family: same-domain face unification, collinear edge
//! merging, and tolerance reduction. Defined in `ogeom-algo`, where the
//! booleans reach them too.

pub(crate) use ogeom_algo::measure_edge;
pub use ogeom_algo::{
    UnifyOptions, merge_edges, reduce_tolerances, same_domain, unify_same_domain,
    unify_same_domain_around, unify_same_domain_with,
};
