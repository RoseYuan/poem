# Compute fuzzy-hard versions of pair-sorting partition metrics

Computes fuzzy-hard versions of pair-sorting partition metrics to
compare a hard clustering with both a fuzzy and hard truth. This was
especially designed for cases where the fuzzy truth represents an
uncertainty of a hard truth. Briefly put, the maximum of the pair
concordance between the clustering and either the hard or the fuzzy
truth is used, and the hard truth is used to compute completeness. See
[`fuzzyPartitionMetrics`](https://roseyuan.github.io/poem/reference/fuzzyPartitionMetrics.md)
for the more standard implementation of the metrics.

## Usage

``` r
fuzzyHardMetrics(
  hardTrue,
  fuzzyTrue,
  hardPred,
  nperms = NULL,
  returnElementPairAccuracy = FALSE,
  lowMemory = NULL,
  verbose = TRUE,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- hardTrue:

  An atomic vector coercible to a factor or integer vector containing
  the true hard labels. Must have the same length as `hardPred`.

- fuzzyTrue:

  A object coercible to a numeric matrix with membership probability of
  elements (rows) in clusters (columns). Must have the same number of
  rows as the length of `hardTrue`. Also note that the columns of
  `fuzzyTrue` should be in the order of the levels (or integer values)
  of `hardTrue`.

- hardPred:

  An atomic vector coercible to a factor or integer vector containing
  the predicted hard labels.

- nperms:

  The number of permutations (for correction for chance). If NULL
  (default), a first set of 10 permutations will be run to estimate
  whether the variation across permutations is above 0.0025, in which
  case more (max 1000) permutations will be run.

- returnElementPairAccuracy:

  Logical. If TRUE, returns the per-element pair accuracy instead of the
  various parition-level and dataset-level metrics. Default FALSE.

- lowMemory:

  Logical; whether to use the slower, low-memory algorithm. By default
  this is enabled if the projected memory usage is higher than ~2GB.

- verbose:

  Logical; whether to print info and warnings, including the standard
  error of the mean across permutations (giving an idea of the precision
  of the adjusted metrics).

- BPPARAM:

  BiocParallel params for multithreading (default none)

## Value

A list of metrics:

- NDC:

  Hullermeier's NDC (fuzzy rand index)

- ACI:

  Ambrosio's Adjusted Concordance Index (ACI), i.e. a permutation-based
  fuzzy version of the adjusted Rand index.

- fuzzyWH:

  Fuzzy Wallace Homogeneity index

- fuzzyWC:

  Fuzzy Wallace Completeness index

- fuzzyAWH:

  Adjusted fuzzy Wallace Homogeneity index

- fuzzyAWC:

  Adjusted fuzzy Wallace Completeness index

## References

Hullermeier et al. 2012; 10.1109/TFUZZ.2011.2179303;

D'Ambrosio et al. 2021; 10.1007/s00357-020-09367-0

## See also

[`fuzzyPartitionMetrics()`](https://roseyuan.github.io/poem/reference/fuzzyPartitionMetrics.md).

## Author

Pierre-Luc Germain

## Examples

``` r
# generate a fuzzy truth:
fuzzyTrue <- matrix(c(
  0.95, 0.025, 0.025, 
  0.98, 0.01, 0.01, 
  0.96, 0.02, 0.02, 
  0.95, 0.04, 0.01, 
  0.95, 0.01, 0.04, 
  0.99, 0.005, 0.005, 
  0.025, 0.95, 0.025, 
  0.97, 0.02, 0.01, 
  0.025, 0.025, 0.95), 
  ncol = 3, byrow=TRUE)
# a hard truth:
hardTrue <- apply(fuzzyTrue,1,FUN=which.max)
# some predicted labels:
hardPred <- c(1,1,1,1,1,1,2,2,2)
fuzzyHardMetrics(hardTrue, fuzzyTrue, hardPred, nperms=3)
#> Standard error of the mean NDC across permutations:0.0352
#> You might want to increase the number of permutations to increase the robustness of the adjusted metrics.
#> $NDC
#> [1] 0.7581944
#> 
#> $ACI
#> [1] 0.5394992
#> 
#> $fuzzyWH
#> $fuzzyWH$global
#> [1] 0.8436111
#> 
#> $fuzzyWH$perPartition
#>          1          2 
#> 1.00000000 0.06166667 
#> 
#> 
#> $fuzzyWC
#> $fuzzyWC$global
#> [1] 0.7322727
#> 
#> $fuzzyWC$perPartition
#>         1         2         3 
#> 0.7195238 1.0000000 1.0000000 
#> 
#> 
#> $fuzzyAWH
#> $fuzzyAWH$global
#> [1] 0.6403322
#> 
#> $fuzzyAWH$perPartition
#>         1         2 
#>  1.000000 -3.456464 
#> 
#> 
#> $fuzzyAWC
#> $fuzzyAWC$global
#> [1] 0.4682516
#> 
#> $fuzzyAWC$perPartition
#>         1         2         3 
#> 0.4682516       NaN       NaN 
#> 
#> 
```
