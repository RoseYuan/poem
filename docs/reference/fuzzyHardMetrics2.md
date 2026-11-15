# Compute fuzzy-hard metrics with lower memory requirement

This is a slightly slower, but low-memory version of
[`fuzzyHardMetrics`](https://roseyuan.github.io/poem/reference/fuzzyHardMetrics.md).

## Usage

``` r
fuzzyHardMetrics2(
  hardTrue,
  fuzzyTrue,
  hardPred,
  nperms = 10,
  returnElementPairAccuracy = FALSE,
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

[`fuzzyHardMetrics()`](https://roseyuan.github.io/poem/reference/fuzzyHardMetrics.md)

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
poem:::fuzzyHardMetrics2(hardTrue, fuzzyTrue, hardPred, nperms=3)
#> Standard error of the mean NDC across permutations:0.00016
#> $NDC
#> [1] 0.7581944
#> 
#> $ACI
#> [1] 0.4675841
#> 
#> $fuzzyWH
#> $fuzzyWH$global
#> [1] 0.7581944
#> 
#> $fuzzyWH$perPartition
#>         1         2 
#> 0.8772917 0.5200000 
#> 
#> 
#> $fuzzyWC
#> $fuzzyWC$global
#> [1] 0.7647297
#> 
#> $fuzzyWC$perPartition
#>         1         2         3 
#> 0.7558929 0.7922222 0.7922222 
#> 
#> 
#> $fuzzyAWH
#> $fuzzyAWH$global
#> [1] 0.4675841
#> 
#> $fuzzyAWH$perPartition
#>         1         2 
#> 0.6921871 0.1507556 
#> 
#> 
#> $fuzzyAWC
#> $fuzzyAWC$global
#> [1] 0.4675841
#> 
#> $fuzzyAWC$perPartition
#>         1         2         3 
#> 0.4763455 0.3407756 0.5024390 
#> 
#> 
```
