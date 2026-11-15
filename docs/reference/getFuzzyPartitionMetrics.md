# Compute external metrics for fuzzy clusterings

Computes a selection of external fuzzy clustering evaluation metrics.

## Usage

``` r
getFuzzyPartitionMetrics(
  hardTrue = NULL,
  fuzzyTrue = NULL,
  hardPred = NULL,
  fuzzyPred = NULL,
  metrics = c("fuzzyWH", "fuzzyAWH", "fuzzyWC", "fuzzyAWC"),
  level = "class",
  nperms = NULL,
  verbose = TRUE,
  BPPARAM = BiocParallel::SerialParam(),
  useNegatives = TRUE,
  usePairs = NULL,
  lowMemory = NULL,
  ...
)
```

## Arguments

- hardTrue:

  An atomic vector coercible to a factor or integer vector containing
  the true hard labels.

- fuzzyTrue:

  A object coercible to a numeric matrix with membership probability of
  elements (rows) in clusters (columns).

- hardPred:

  An atomic vector coercible to a factor or integer vector containing
  the predicted hard labels.

- fuzzyPred:

  A object coercible to a numeric matrix with membership probability of
  elements (rows) in clusters (columns).

- metrics:

  The metrics to compute. See details.

- level:

  The level to calculate the metrics. Options include `"element"`,
  `"class"` and `"dataset"`.

- nperms:

  The number of permutations (for correction for chance). If NULL
  (default), a first set of 10 permutations will be run to estimate
  whether the variation across permutations is above 0.0025, in which
  case more (max 1000) permutations will be run.

- verbose:

  Logical; whether to print info and warnings, including the standard
  error of the mean across permutations (giving an idea of the precision
  of the adjusted metrics).

- BPPARAM:

  BiocParallel params for multithreading (default none)

- useNegatives:

  Logical; whether to include negative pairs in the concordance score
  (tends to result in a larger overall concordance and lower dynamic
  range of the score). Default TRUE.

- usePairs:

  Logical; whether to compute over pairs instead of elements Recommended
  and TRUE by default.

- lowMemory:

  Logical, whether to use a low memory mode. This is only useful
  when`hardTrue` and `fuzzyPred` is used. If `TRUE`, the function will
  compute the metrics in a low memory mode, which is slower but uses
  less memory. If `FALSE`, the function will compute the metrics in a
  high memory mode, which is faster but uses more memory. By default it
  is set automatically based on the size of the input data. See
  [`fuzzyHardMetrics`](https://roseyuan.github.io/poem/reference/fuzzyHardMetrics.md).

- ...:

  Optional arguments for
  [`fuzzyPartitionMetrics`](https://roseyuan.github.io/poem/reference/fuzzyPartitionMetrics.md):
  `tnorm`. Only useful when `fuzzyTrue` and `fuzzyPred` is used.

## Value

A dataframe of metric results.

## Details

The allowed values for `metrics` depend on the value of `level`:

- If `level = "element"`, the allowed `metrics` are: `"fuzzySPC"`.

- If `level = "class"`, the allowed `metrics` are: `"fuzzyWH"`,
  `"fuzzyAWH"`, `"fuzzyWC"`, `"fuzzyAWC"`.

- If `level = "dataset"`, the allowed `metrics` are: `"fuzzyRI"`,
  `"fuzzyARI"`, `"fuzzyWH"`, `"fuzzyAWH"`, `"fuzzyWC"`, `"fuzzyAWC"`.

## Examples

``` r
# generate fuzzy partitions:
m1 <- matrix(c(0.95, 0.025, 0.025, 
               0.98, 0.01, 0.01, 
               0.96, 0.02, 0.02, 
               0.95, 0.04, 0.01, 
               0.95, 0.01, 0.04, 
               0.99, 0.005, 0.005, 
               0.025, 0.95, 0.025, 
               0.97, 0.02, 0.01, 
               0.025, 0.025, 0.95), 
               ncol = 3, byrow=TRUE)
m2 <- matrix(c(0.95, 0.025, 0.025,  
               0.98, 0.01, 0.01, 
               0.96, 0.02, 0.02, 
               0.025, 0.95, 0.025, 
               0.02, 0.96, 0.02, 
               0.01, 0.98, 0.01, 
               0.05, 0.05, 0.95, 
               0.02, 0.02, 0.96, 
               0.01, 0.01, 0.98), 
               ncol = 3, byrow=TRUE)
colnames(m1) <- colnames(m2) <- LETTERS[seq_len(3)]
getFuzzyPartitionMetrics(fuzzyTrue=m1,fuzzyPred=m2, level="class")
#> Comparing between a fuzzy truth and a fuzzy prediction...
#> Running 100 extra permutations.
#> Standard error of the mean NDC across permutations:0.00215
#>     fuzzyWC    fuzzyAWC class   fuzzyWH   fuzzyAWH cluster
#> 1 0.3445840  0.05010597     1        NA         NA      NA
#> 2 0.7242508 -0.02250224     2        NA         NA      NA
#> 3 0.7520319  0.06388584     3        NA         NA      NA
#> 4        NA          NA    NA 0.9359492  0.8286227       1
#> 5        NA          NA    NA 0.9214151  0.8274947       2
#> 6        NA          NA    NA 0.1588990 -1.1160180       3

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
getFuzzyPartitionMetrics(hardPred=hardPred, hardTrue=hardTrue, 
fuzzyTrue=fuzzyTrue, nperms=3, level="class")
#> Comparing between a fuzzy truth and a hard prediction...
#> Standard error of the mean NDC across permutations:0.0947
#> You might want to increase the number of permutations to increase the robustness of the adjusted metrics.
#>     fuzzyWC  fuzzyAWC class    fuzzyWH   fuzzyAWH cluster
#> 1 0.7195238 0.3542847     1         NA         NA      NA
#> 2 1.0000000       NaN     2         NA         NA      NA
#> 3 1.0000000       NaN     3         NA         NA      NA
#> 4        NA        NA    NA 1.00000000  1.0000000       1
#> 5        NA        NA    NA 0.06166667 -0.8006397       2
getFuzzyPartitionMetrics(hardTrue=hardPred, hardPred=hardTrue, 
fuzzyPred=fuzzyTrue, nperms=3, level="class")
#> Comparing between a hard truth and a fuzzy prediction...
#> Standard error of the mean NDC across permutations:0.106
#> You might want to increase the number of permutations to increase the robustness of the adjusted metrics.
#>      fuzzyWC  fuzzyAWC class   fuzzyWH  fuzzyAWH cluster
#> 1 1.00000000  1.000000     1        NA        NA      NA
#> 2 0.06166667 -1.978836     2        NA        NA      NA
#> 3         NA        NA    NA 0.7195238 0.3967224       1
#> 4         NA        NA    NA 1.0000000       NaN       2
#> 5         NA        NA    NA 1.0000000       NaN       3
```
