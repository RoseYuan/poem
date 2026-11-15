# Get Internal Objects

Computes the internal nodes and edges using Minimum Spanning Tree.

## Usage

``` r
.get_internal_objects(mutual_reach_dists, use_igraph_mst = TRUE)
```

## Arguments

- mutual_reach_dists:

  Numeric matrix representing mutual reachability distances.

- use_igraph_mst:

  Logical flag to use MST implementation in igraph. Currently only mst
  from igraph is implemented.

## Value

A list containing the indices of internal nodes and their edge weights.
