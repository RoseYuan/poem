# Compute class-level internal evaluation metrics for spatially-resolved data

Computes a selection of internal clustering evaluation metrics for
spatial data for each class.

## Usage

``` r
getSpatialClassInternalMetrics(
  labels,
  location,
  k = 6,
  metrics = c("CHAOS", "PAS", "ELSA"),
  ...
)
```

## Arguments

- labels:

  A vector containing the labels to be evaluated.

- location:

  A numerical matrix containing the location information, with rows as
  samples and columns as location dimensions.

- k:

  The size of the spatial neighborhood to look at for each spot. This is
  used for calculating PAS and ELSA scores.

- metrics:

  Possible metrics: "CHAOS", "PAS" and "ELSA".

- ...:

  Optional params for
  [`PAS()`](https://roseyuan.github.io/poem/reference/PAS.md).

## Value

A dataframe of metric values.
