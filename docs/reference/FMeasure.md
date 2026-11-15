# Calculate F measure

Compute the F measure between two clustering results. This is directly
copied from the package `FlowSOM`.

## Usage

``` r
FMeasure(true, pred, silent = TRUE)
```

## Arguments

- true:

  Array containing real cluster labels for each sample

- pred:

  Array containing predicted cluster labels for each sample

- silent:

  Logical, if FALSE, print some information about precision and recall

## Value

F measure score
