# Compute partition-based metrics

Computes a selection of external evaluation metrics for partition.

## Usage

``` r
getPartitionMetrics(true, pred, metrics = NULL, level = "class", ...)
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
  details.

- level:

  The level to calculate the metrics. Options include "element",
  `"class"` and `"dataset"`.

- ...:

  Optional arguments for MI, VI, or VM. See
  [`clevr::mutual_info()`](https://rdrr.io/pkg/clevr/man/mutual_info.html),
  [`clevr::variation_info()`](https://rdrr.io/pkg/clevr/man/variation_info.html)
  and
  [`clevr::v_measure()`](https://rdrr.io/pkg/clevr/man/v_measure.html)
  for more details.

## Value

A data.frame of metrics.

## Details

The allowed values for `metrics` depend on the value of `level`:

- If `level = "element"`, the allowed `metrics` are:

  - `"SPC"`: Spot-wise Pair Concordance.

  - `"ASPC"`: Adjusted Spot-wise Pair Concordance.

- If `level = "class"`, the allowed `metrics` are:
  `"WC"`,`"WH"`,`"AWC"`, `"AWH"`,`"FM"` (see below for details).

- If `level = "dataset"`, the allowed `metrics` are:

  - `"RI"`: Rand Index

  - `"WC"`: Wallace Completeness

  - `"WH"`: Wallace Homogeneity

  - `"ARI"`: Adjusted Rand Index

  - `"AWC"`: Adjusted Wallace Completeness

  - `"AWH"`: Adjusted Wallace Homogeneity

  - `"NCR"`: Normalized class size Rand index

  - `"MI"`: Mutual Information

  - `"AMI"`: Adjusted Mutual Information

  - `"VI"`: Variation of Information

  - `"EH"`: (Entropy-based) Homogeneity

  - `"EC"`: (Entropy-based) Completeness

  - `"VM"`: V-measure

  - `"FM"`: F-measure/weighted average F1 score

  - `"VDM"`: Van Dongen Measure

  - `"MHM"`: Meila-Heckerman Measure

  - `"MMM"`: Maximum-Match Measure

  - `"Mirkin"`: Mirkin Metric

  - `"Accuracy"`: Set Matching Accuracy

## Examples

``` r
true <- rep(LETTERS[seq_len(3)], each=10)
pred <- c(rep("A", 8), rep("B", 9), rep("C", 3), rep("D", 10))
getPartitionMetrics(true, pred, level="class")
#> Warning: `sortPairs()` was deprecated in aricode 1.1.0.
#> ℹ Please use `sort_pairs()` instead.
#> ℹ The deprecated feature was likely used in the poem package.
#>   Please report the issue at <https://github.com/RoseYuan/poem/issues>.
#>          WC       AWC        FM class        WH       AWH cluster
#> 1 0.6444444 0.5211558 0.8888889     A        NA        NA    <NA>
#> 2 0.5333333 0.3715170 0.7368421     B        NA        NA    <NA>
#> 3 1.0000000 1.0000000 1.0000000     C        NA        NA    <NA>
#> 4        NA        NA        NA  <NA> 1.0000000 1.0000000       A
#> 5        NA        NA        NA  <NA> 0.6111111 0.4361111       B
#> 6        NA        NA        NA  <NA> 1.0000000 1.0000000       C
#> 7        NA        NA        NA  <NA> 1.0000000 1.0000000       D
getPartitionMetrics(true, pred, level="dataset")
#>          RI        WC    WH       ARI       NCR      AWC     AWH        MI
#> 1 0.8827586 0.7259259 0.875 0.7126481 0.8396296 0.630891 0.81875 0.9397004
#>         AMI        VI        EH        EC        VM        FM
#> 1 0.6887211 0.5293341 0.8553522 0.7172614 0.7802439 0.8752437
```
