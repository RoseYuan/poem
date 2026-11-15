# Density Separation of a Pair of Clusters

Computes the density separation between two clusters.

## Usage

``` r
.fn_density_separation(
  cls_i,
  cls_j,
  dists,
  internal_core_dists_i,
  internal_core_dists_j
)
```

## Arguments

- cls_i:

  Integer, first cluster index.

- cls_j:

  Integer, second cluster index.

- dists:

  Numeric matrix of distances.

- internal_core_dists_i:

  Numeric vector of core distances for cluster i.

- internal_core_dists_j:

  Numeric vector of core distances for cluster j.

## Value

A list containing the cluster indices and their density separation.
