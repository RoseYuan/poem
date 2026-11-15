# Compute spot-level internal evaluation metrics for spatially-resolved data

Computes a selection of internal clustering evaluation metrics for
spatial data at each spot level.

## Usage

``` r
getSpatialElementInternalMetrics(
  labels,
  location,
  k = 6,
  metrics = c("PAS", "ELSA"),
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

  Possible metrics: "PAS" and "ELSA".

- ...:

  Optional params for
  [`PAS()`](https://roseyuan.github.io/poem/reference/PAS.md).

## Value

A dataframe containing the metric values for all samples in the dataset.
If PAS is calculated, the value is a Boolean about the abnormality of a
spot. If ELSA is calculated, Ea, Ec and ELSA for all spots will be
returned.
