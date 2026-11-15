# Density Sparseness of a Cluster

Computes the density sparseness for a given cluster.

## Usage

``` r
.fn_density_sparseness(cls_inds, dists, d, use_igraph_mst)
```

## Arguments

- cls_inds:

  Integer vector of cluster indices.

- dists:

  Numeric matrix of distances.

- d:

  Integer, the dimensionality.

- use_igraph_mst:

  Logical flag to use MST implementation in igraph. Currently only mst
  from igraph is implemented.

## Value

A list containing the density sparseness, internal core distances, and
internal node indices.
