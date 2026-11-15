# Compute dataset-level external evaluation metrics for spatially-resolved data

Computes a selection of external clustering evaluation metrics for
spatial data at the dataset level. Options include a series of fuzzy
pair-counting metrics and set matching-based accuracy.

## Usage

``` r
getSpatialGlobalExternalMetrics(
  true,
  pred,
  location,
  k = 6,
  alpha = 0.5,
  metrics = c("nsRI", "nsARI", "nsWH", "nsAWH", "nsWC", "nsAWC", "nsAccuracy",
    "SpatialRI", "SpatialARI"),
  fuzzy_true = TRUE,
  fuzzy_pred = FALSE,
  lowMemory = NULL,
  ...
)
```

## Arguments

- true:

  A vector containing the labels of the true classes. Must be a vector
  of characters, integers, numerics, or a factor, but not a list.

- pred:

  A vector containing the labels of the predicted clusters. Must be a
  vector of characters, integers, numerics, or a factor, but not a list.

- location:

  A matrix or data.frame of coordinates

- k:

  The number of neighbors used when calculating the fuzzy class
  memberships for fuzzy metrics, or when calculating the weighted
  accuracy.

- alpha:

  The parameter to control to what extend the spot itself contribute to
  the class composition calculation. `"equal"` means it is weighted the
  same as other neighbors. A numeric value between `0` and `1` means the
  weight of the frequency contribution for the spot itself, and the
  frequency contribution for its knn is then `1-alpha`. By default
  `0.5`.

- metrics:

  a vector of metric names to compute.

- fuzzy_true:

  Logical; whether to compute fuzzy class memberships for `true`.

- fuzzy_pred:

  Logical; whether to compute fuzzy class memberships for `pred`.

- lowMemory:

  Logical, whether to use a low memory mode. This is only useful
  when`hardTrue` and `fuzzyPred` is used. If `TRUE`, the function will
  compute the metrics in a low memory mode, which is slower but uses
  less memory. If `FALSE`, the function will compute the metrics in a
  high memory mode, which is faster but uses more memory. By default it
  is set automatically based on the size of the input data. See
  [`fuzzyHardMetrics`](https://roseyuan.github.io/poem/reference/fuzzyHardMetrics.md).

- ...:

  Optional params for
  [`fuzzyPartitionMetrics`](https://roseyuan.github.io/poem/reference/fuzzyPartitionMetrics.md)
  or
  [`findSpatialKNN`](https://roseyuan.github.io/poem/reference/findSpatialKNN.md).

## Value

A data.frame of metrics.
