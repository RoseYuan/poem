# silhouetteWidths

Computes the silhouette widths. If the dataset is sufficiently small for
the
[`cluster::silhouette`](https://rdrr.io/pkg/cluster/man/silhouette.html)
implementation to work, this will be used. Otherwise a slower chunked
implementation is used.

## Usage

``` r
silhouetteWidths(x, labels)
```

## Arguments

- x:

  A numeric matrix or data.frame with observations as rows.

- labels:

  An integer/factor vector of clustering labels, or length equal to the
  number of rows in `x`.

## Value

A numeric vector of silhouette widths for each element of `x`.

## Examples

``` r
# generate dummy data
m <- matrix(rnorm(100*3),ncol=3)
labels <- sample.int(3,100,replace=TRUE)
# calculate SWs:
sw <- silhouetteWidths(m, labels)
```
