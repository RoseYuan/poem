# getSpatialElementExternalMetrics

Computes a selection of external clustering evaluation metrics for
spatial data at the element level.

## Usage

``` r
getSpatialElementExternalMetrics(
  true,
  pred,
  location,
  k = 6,
  alpha = 0.5,
  metrics = c("nsSPC", "NPC", "SpatialSPC"),
  fuzzy_true = TRUE,
  fuzzy_pred = FALSE,
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

- ...:

  Optional params for
  [`getFuzzyPartitionElementMetrics()`](https://roseyuan.github.io/poem/reference/getFuzzyPartitionElementMetrics.md)
  or
  [`findSpatialKNN()`](https://roseyuan.github.io/poem/reference/findSpatialKNN.md).

## Value

A data.frame of metrics.
