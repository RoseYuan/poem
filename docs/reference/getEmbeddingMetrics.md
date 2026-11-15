# Compute embedding-based metrics

Computes embedding-based metrics for the specified level.

## Usage

``` r
getEmbeddingMetrics(
  x,
  labels,
  metrics = NULL,
  distance = "euclidean",
  level = "class",
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

  The metrics to compute. See details.

- distance:

  The distance metric to use (default euclidean).

- level:

  The level to calculate the metrics. Options include `"element"`,
  `"class"` and `"dataset"`.

- ...:

  Optional arguments. See details.

## Value

A data.frame of metrics.

## Details

The allowed values for `metrics` depend on the value of `level`:

- If `level = "element"`, the allowed `metrics` are: `"SW"`.

- If `level = "class"`, the allowed `metrics` are: `"meanSW"`,
  `"minSW"`, `"pnSW"`, `"dbcv"`.

- If `level = "dataset"`, the allowed `metrics` are: `"meanSW"`,
  `"meanClassSW"`, `"pnSW"`, `"minClassSW"`, `"cdbw"`, `"cohesion"`,
  `"compactness"`, `"sep"`, `"dbcv"`.

The function(s) that the optional arguments `...` passed to depend on
the value of `level`:

- If `level = "element"`, optional arguments are passed to
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html).

- If `level = "class"`, optional arguments are passed to
  [`dbcv()`](https://roseyuan.github.io/poem/reference/dbcv.md).

- If `level = "dataset"`, optional arguments are passed to
  [`dbcv()`](https://roseyuan.github.io/poem/reference/dbcv.md) or
  [`CDbw()`](https://roseyuan.github.io/poem/reference/CDbw.md).

## Examples

``` r
d1 <- mockData()
getEmbeddingMetrics(d1[,seq_len(2)], labels=d1$class, 
metrics=c("meanSW", "minSW", "pnSW", "dbcv"), level="class")
#>   class    meanSW      minSW pnSW       dbcv
#> 1     A 0.3334155 -0.1926139 0.12 -0.1250971
#> 2     B 0.1589520 -0.4874509 0.20 -0.4673145
```
