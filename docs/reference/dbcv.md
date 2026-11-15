# Calculate DBCV Metric

Compute the DBCV (Density-Based Clustering Validation) metric.

## Usage

``` r
dbcv(
  X,
  labels,
  distance = "euclidean",
  noise_id = -1,
  check_duplicates = FALSE,
  use_igraph_mst = TRUE,
  BPPARAM = BiocParallel::SerialParam(),
  ...
)
```

## Arguments

- X:

  Numeric matrix of samples.

- labels:

  Integer vector of cluster IDs.

- distance:

  String specifying the distance metric. `"sqeuclidean"`, or possible
  `method` in [`stats::dist()`](https://rdrr.io/r/stats/dist.html). By
  default `"euclidean"`.

- noise_id:

  Integer, the cluster ID in `y` for noise (default `-1`).

- check_duplicates:

  Logical flag to check for duplicate samples.

- use_igraph_mst:

  Logical flag to use `igraph`'s MST implementation. Currently only
  `mst` from `igraph` is implemented.

- BPPARAM:

  BiocParallel params for multithreading (default none)

- ...:

  Ignored

## Value

A list:

- vcs:

  Numeric vector of validity index for each cluster.

- dbcv:

  Numeric value representing the overall DBCV metric.

## Details

This implementation will not fully reproduce the results of other
existing implementations (e.g. <https://github.com/FelSiq/DBCV>) due to
the different algorithms used for computing the Minimum Spanning Tree.

## References

Davoud Moulavi, et al. 2014; 10.1137/1.9781611973440.96.

## Examples

``` r
data(noisy_moon)
data <- noisy_moon
dbcv(data[, c("x", "y")], data$kmeans_label)
#> $vcs
#> [1] -0.4383721 -0.4077112
#> 
#> $dbcv
#> [1] -0.4230416
#> 
dbcv(data[, c("x", "y")], data$hdbscan_label)
#> $vcs
#>  [1] -0.5889023  0.3726825  0.5500422  0.7884686  0.4887283  0.7682203
#>  [7]  0.7246492  0.7246492  0.9349664  0.4749650  0.5500422  0.3726825
#> 
#> $dbcv
#> [1] 0.4214685
#> 
```
