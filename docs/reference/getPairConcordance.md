# Per-element pair concordance score

Per-element pair concordance between a clustering and a ground truth.
Note that by default, negative pairs (i.e. that are split in both the
predicted and true groupings) are not counted. To count it (as in the
standard Rand Index), use `useNegatives=TRUE`.

## Usage

``` r
getPairConcordance(
  true,
  pred,
  usePairs = TRUE,
  useNegatives = FALSE,
  adjust = FALSE
)
```

## Arguments

- true:

  A vector of true class labels

- pred:

  A vector of predicted clusters

- usePairs:

  Logical; whether to compute over pairs instead of elements Recommended
  and TRUE by default.

- useNegatives:

  Logical; whether to include the consistency of negative pairs in the
  score (default FALSE).

- adjust:

  Logical; whether to adjust for chance. Only implemented for
  `useNegatives=FALSE` (doesn't make sense on a element-level
  otherwise).

## Value

A vector of concordance scores
