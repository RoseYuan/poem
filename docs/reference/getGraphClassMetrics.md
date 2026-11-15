# getGraphClassMetrics

Computes a selection of supervised graph evaluation metrics using ground
truth class labels. The metrics are reported (as average) per class.

## Usage

``` r
getGraphClassMetrics(
  x,
  labels,
  metrics = c("SI", "NP", "AMSP", "PWC", "NCE"),
  directed = NULL,
  ...
)

# S4 method for class 'list'
getGraphClassMetrics(x, labels, metrics, directed = NULL, k = NULL, ...)

# S4 method for class 'data.frame'
getGraphClassMetrics(
  x,
  labels,
  metrics,
  directed = NULL,
  k,
  shared = FALSE,
  ...
)

# S4 method for class 'matrix'
getGraphClassMetrics(
  x,
  labels,
  metrics,
  directed = NULL,
  k,
  shared = FALSE,
  ...
)

# S4 method for class 'igraph'
getGraphClassMetrics(
  x,
  labels,
  metrics = c("SI", "NP", "AMSP", "PWC", "NCE"),
  directed = NULL,
  ...
)

# S4 method for class 'dist'
getGraphClassMetrics(
  x,
  labels,
  metrics = c("SI", "NP", "AMSP", "PWC", "NCE"),
  directed = NULL,
  ...
)
```

## Arguments

- x:

  Either an igraph object, a list of nearest neighbors (see details
  below), or a data.frame or matrix (with features as columns and items
  as rows) from which nearest neighbors will be computed.

- labels:

  Either a factor or a character vector indicating the true class label
  of each element (i.e. row or vertex) of `x`.

- metrics:

  The metrics to compute. See details.

- directed:

  Logical; whether to compute the metrics in a directed fashion. If left
  to NULL, conventional choices will be made per metric (adhesion,
  cohesion, PWC AMSP undirected, others directed).

- ...:

  Optional arguments for
  [`emb2knn()`](https://roseyuan.github.io/poem/reference/emb2knn.md) or
  [`emb2snn()`](https://roseyuan.github.io/poem/reference/emb2snn.md).

- k:

  The number of nearest neighbors to compute and/or use. Can be omitted
  if `x` is a graph or list of nearest neighbors.

- shared:

  Logical; whether to use a shared nearest neighbor network instead of a
  nearest neighbor network. Ignored if `x` is not an embedding or dist
  object.

## Value

A data.frame of metrics for each class.
