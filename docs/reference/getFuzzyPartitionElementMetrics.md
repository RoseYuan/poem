# getFuzzyPartitionElementMetrics

Computes a selection of external fuzzy clustering evaluation metrics at
the element level.

## Usage

``` r
getFuzzyPartitionElementMetrics(
  hardTrue = NULL,
  fuzzyTrue = NULL,
  hardPred = NULL,
  fuzzyPred = NULL,
  fuzzy_true = TRUE,
  fuzzy_pred = FALSE,
  metrics = c("fuzzySPC"),
  useNegatives = TRUE,
  verbose = TRUE,
  usePairs = TRUE
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

- fuzzyPred:

  A object coercible to a numeric matrix with membership probability of
  elements (rows) in clusters (columns).

- fuzzy_true:

  Logical; whether the truth is fuzzy.

- fuzzy_pred:

  Logical; whether the prediction is fuzzy.

- metrics:

  The metrics to compute. Currently only `"fuzzySPC"` is included at the
  element level.

- useNegatives:

  Logical; whether to include negative pairs in the concordance score
  (tends to result in a larger overall concordance and lower dynamic
  range of the score). Default TRUE.

- verbose:

  Logical; whether to print expected memory usage for large datasets.

- usePairs:

  Logical; whether to compute over pairs instead of elements. Only
  useful when `fuzzy_true=TRUE` and `fuzzy_pred=FALSE`.

## Value

A dataframe of metric values.
