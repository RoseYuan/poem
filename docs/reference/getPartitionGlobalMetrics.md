# getPartitionGlobalMetrics

Computes a selection of external evaluation metrics for partition. The
metrics are reported per dataset.

## Usage

``` r
getPartitionGlobalMetrics(
  true,
  pred,
  metrics = c("RI", "WC", "WH", "ARI", "NCR", "AWC", "AWH", "MI", "AMI", "VI", "EH",
    "EC", "VM", "FM"),
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

- metrics:

  The metrics to compute. If omitted, main metrics will be computed. See
  below for more details.

- ...:

  Optional arguments for MI, VI, or VM. See
  [`clevr::mutual_info()`](https://rdrr.io/pkg/clevr/man/mutual_info.html),
  [`clevr::variation_info()`](https://rdrr.io/pkg/clevr/man/variation_info.html)
  and
  [`clevr::v_measure()`](https://rdrr.io/pkg/clevr/man/v_measure.html)
  for more details.

## Value

A dataframe of metric results. Possible metrics are:

- RI:

  Rand Index

- WC:

  Wallace Completeness

- WH:

  Wallace Homogeneity

- ARI:

  Adjusted Rand Index

- AWC:

  Adjusted Wallace Completeness

- AWH:

  Adjusted Wallace Homogeneity

- NCR:

  Normalized class size Rand index

- MI:

  Mutual Information

- AMI:

  Adjusted Mutual Information

- VI:

  Variation of Information

- EH:

  (Entropy-based) Homogeneity

- EC:

  (Entropy-based) Completeness

- VM:

  V-measure

- FM:

  F-measure/weighted average F1 score

- VDM:

  Van Dongen Measure

- MHM:

  Meila-Heckerman Measure

- MMM:

  Maximum-Match Measure

- Mirkin:

  Mirkin Metric

- Accuracy:

  Set Matching Accuracy
