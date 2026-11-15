# Compute fuzzy-fuzzy versions of pair-sorting partition metrics

Computes fuzzy versions of pair-sorting partition metrics. This is
largely based on the permutation-based implementation by Antonio
D'Ambrosio from the `ConsRankClass` package, modified to also compute
the fuzzy versions of the adjusted Wallace indices, implement
multithreading, and adjust the number of permutations according to their
variability.

## Usage

``` r
fuzzyPartitionMetrics(
  P,
  Q,
  computeWallace = TRUE,
  nperms = NULL,
  verbose = TRUE,
  returnElementPairAccuracy = FALSE,
  BPPARAM = BiocParallel::SerialParam(),
  tnorm = c("product", "min", "lukasiewicz")
)
```

## Arguments

- P:

  A object coercible to a numeric matrix with membership probability of
  elements (rows) in ground-truth classes (columns).

- Q:

  A object coercible to a numeric matrix with membership probability of
  elements (rows) in predicted clusters (columns). Must have the same
  number of rows as `P`.

- computeWallace:

  Logical; whether to compute the individual fuzzy versions of the
  Wallace indices (increases running time).

- nperms:

  The number of permutations (for correction for chance). If NULL
  (default), a first set of 10 permutations will be run to estimate
  whether the variation across permutations is above 0.0025, in which
  case more (max 1000) permutations will be run.

- verbose:

  Logical; whether to print info and warnings, including the standard
  error of the mean across permutations (giving an idea of the precision
  of the adjusted metrics).

- returnElementPairAccuracy:

  Logical. If TRUE, returns the per-element pair accuracy instead of the
  various parition-level and dataset-level metrics. Default FALSE.

- BPPARAM:

  BiocParallel params for multithreading (default none)

- tnorm:

  Which type of t-norm operation to use for class membership of pairs
  (either product, min, or lukasiewicz) when calculating the Wallace
  indices. Does not influence the NDC/ACI metrics.

## Value

When `returnElementPairAccuracy` is `FALSE`, return a list of metrics:

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

## Author

Pierre-Luc Germain

## Examples

``` r
# generate fuzzy partitions:
```
