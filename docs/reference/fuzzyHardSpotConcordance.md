# Per-element maximal concordance between a hard and a fuzzy partition

Per-element maximal concordance between a hard clustering and hard and
fuzzy ground truth labels.

## Usage

``` r
fuzzyHardSpotConcordance(
  hardTrue,
  fuzzyTrue,
  hardPred,
  useNegatives = TRUE,
  verbose = TRUE
)
```

## Arguments

- hardTrue:

  A vector of true cluster labels

- fuzzyTrue:

  A object coercible to a numeric matrix with membership probability of
  elements (rows) in clusters (columns). Must have the same number of
  rows as the length of `hardTrue`.

- hardPred:

  A vector of predicted cluster labels

- useNegatives:

  Logical; whether to include negative pairs in the concordance score
  (tends to result in a larger overall concordance and lower dynamic
  range of the score). Default TRUE.

- verbose:

  Logical; whether to print expected memory usage for large datasets.

## Value

A numeric vector of concordance scores for each element of `hardPred`

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
fuzzyHardSpotConcordance(hardTrue, fuzzyTrue, hardPred)
#> [1] 0.87750 0.87625 0.87625 0.87750 0.87875 0.87750 0.76625 0.02750 0.76625
```
