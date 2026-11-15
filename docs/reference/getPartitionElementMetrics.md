# getPartitionElementMetrics

Computes a selection of external evaluation metrics for partition. The
metrics are reported per element.

## Usage

``` r
getPartitionElementMetrics(
  true,
  pred,
  metrics = c("SPC"),
  usePairs = TRUE,
  useNegatives = TRUE
)
```

## Arguments

- true:

  A vector of true class labels

- pred:

  A vector of predicted clusters

- metrics:

  The metrics to compute.

- usePairs:

  Logical; whether to compute over pairs instead of elements Recommended
  and TRUE by default.

- useNegatives:

  Logical; whether to include the consistency of negative pairs in the
  score (default FALSE).

## Value

A dataframe of metrics.
