# getPartitionClassMetrics

Computes a selection of external evaluation metrics for partition. The
metrics are reported per class.

## Usage

``` r
getPartitionClassMetrics(
  true,
  pred,
  metrics = c("WC", "WH", "AWC", "AWH", "FM")
)
```

## Arguments

- true:

  A vector containing the labels of the true classes. Must be a vector
  of characters, integers, numerics, or a factor, but not a list.

- pred:

  A vector containing the labels of the predicted clusters. Must be a
  vector of characters, integers, numerics, or a factor, but not a list.

- metrics:

  The metrics to compute. If omitted, main metrics will be computed.

## Value

A dataframe of metrics.
