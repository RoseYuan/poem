# Compute graph-based metrics

Computes a selection of graph evaluation metrics using class labels.

## Usage

``` r
getGraphMetrics(
  x,
  labels,
  metrics = NULL,
  directed = NULL,
  k = 10,
  shared = FALSE,
  level = "class",
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

- k:

  The number of nearest neighbors to compute and/or use. Can be omitted
  if `x` is a graph or list of nearest neighbors.

- shared:

  Logical; whether to use a shared nearest neighbor network instead of a
  nearest neighbor network. Ignored if `x` is not an embedding or dist
  object.

- level:

  The level to calculate the metrics. Options include `"element"`,
  `"class"` and `"dataset"`.

- ...:

  Optional arguments for
  [`emb2knn()`](https://roseyuan.github.io/poem/reference/emb2knn.md) or
  [`emb2snn()`](https://roseyuan.github.io/poem/reference/emb2snn.md).

## Value

A data.frame of metrics.

## Details

The allowed values for `metrics` depend on the value of `level`:

- If `level = "element"`, the allowed `metrics` are:
  `"SI"`,`"ISI"`,`"NP"`, `"NCE"` (see below for details).

- If `level = "class"`, the allowed `metrics` are:

  - `"SI"`: Simpson’s Index.

  - `"ISI"`: Inverse Simpson’s Index

  - `"NP"`: Neighborhood Purity

  - `"AMSP"`: Adjusted Mean Shortest Path

  - `"PWC"`: Proportion of Weakly Connected

  - `"NCE"`: Neighborhood Class Enrichment

  - `"adhesion"`: adhesion of a graph, is the minumum number of nodes
    that must be removed to split a graph.

  - `"cohesion"`: cohesion of a graph, is the minumum number of edges
    that must be removed to split a graph.

- If `level = "dataset"`, the allowed `metrics` are: `"SI"`,`"ISI"`,
  `"NP"`,`"AMSP"`,`"PWC"`,`"NCE"`, `"adhesion"`,`"cohesion"`.

## Examples

``` r
d1 <- mockData()
getGraphMetrics(d1[,seq_len(2)], labels=d1$class, level="class")
#>   class     SI    NP      AMSP PWC       NCE
#> A     A 0.7192 0.748 0.1136000 0.2 0.1645283
#> B     B 0.5600 0.520 0.1726984 0.6 0.3089594
```
