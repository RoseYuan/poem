# Per-element agreement score

Per-element agreement between a clustering and a ground truth

## Usage

``` r
getAgreement(true, pred, usePairs = TRUE, useNegatives = FALSE, adjust = FALSE)
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

A vector of agreement scores
