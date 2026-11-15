# Per-element concordance between two fuzzy partitions

Per-element concordance between two fuzzy partitions

## Usage

``` r
fuzzySpotConcordance(P, Q)
```

## Arguments

- P:

  A object coercible to a numeric matrix with membership probability of
  elements (rows) in clusters (columns)

- Q:

  A object coercible to a numeric matrix with membership probability of
  elements (rows) in clusters (columns). Must have the same number of
  rows as `P`

## Value

A numeric vector of concordance scores for each row of `P`.

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
fuzzySpotConcordance(m1,m2)
#> Error in fuzzySpotConcordance(m1, m2): could not find function "fuzzySpotConcordance"
```
