# nnWeightedAccuracy

Computes an accuracy score which weighs elements/spots that are
misclassified by the proportion of their (spatial) neighborhood that is
not of the element/spot's predicted class. This reduces the weight of
misclassifications happening at the boundary between domains.

## Usage

``` r
nnWeightedAccuracy(true, pred, location, k = 5, ...)
```

## Arguments

- true:

  True class labels (vector coercible to factor)

- pred:

  Predicted labels (vector coercible to factor)

- location:

  The spatial coordinates to compute the nearest neighbors.

- k:

  Number of nearest neighbors

- ...:

  Optional params passed to
  [`findSpatialKNN()`](https://roseyuan.github.io/poem/reference/findSpatialKNN.md).

## Value

A scalar representing the weighted accuracy.
