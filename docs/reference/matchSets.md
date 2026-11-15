# Match two partitions using Hungarian algorithm

Match sets from a partitions to a reference partition using the
Hungarian algorithm to optimize F1 scores.

## Usage

``` r
matchSets(pred, true, forceMatch = TRUE, returnIndices = is.integer(true))
```

## Arguments

- pred:

  An integer or factor of cluster labels

- true:

  An integer or factor of reference labels

- forceMatch:

  Logical; whether to enforce a match for every set of `pred`

- returnIndices:

  Logical; whether to return indices rather than levels

## Value

A vector of matching sets (i.e. level) from `true` for every set (i.e.
level) of `pred`.
