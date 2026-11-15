# Compute Pair to Pair Distances

Compute the pairwise distances between points in matrix `X`.

## Usage

``` r
.compute_pair_to_pair_dists(X, distance = "euclidean")
```

## Arguments

- X:

  Numeric matrix.

- distance:

  String specifying the metric to compute the distances.

## Value

Numeric matrix of pairwise distances with self-distances set to `Inf`.
