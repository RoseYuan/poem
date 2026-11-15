# getEmbeddingClassMetrics

Computes class-level, embedding-based metrics.

## Usage

``` r
getEmbeddingClassMetrics(
  x,
  labels,
  metrics = c("meanSW", "minSW", "pnSW", "dbcv"),
  distance = "euclidean",
  ...
)
```

## Arguments

- x:

  A data.frame or matrix (with features as columns and items as rows)
  from which the metrics will be computed.

- labels:

  A vector containing the labels of the predicted clusters. Must be a
  vector of characters, integers, numerics, or a factor, but not a list.

- metrics:

  The metrics to compute.

- distance:

  The distance metric to use (default euclidean).

- ...:

  Optional arguments. See details.

## Value

A data.frame of metrics for each node/element of `x`.
