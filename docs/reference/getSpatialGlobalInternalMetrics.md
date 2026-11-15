# Compute dataset-level internal evaluation metrics for spatially-resolved data

Computes a selection of internal clustering evaluation metrics for
spatial data at the dataset level. MPC, PC and PE are internal metrics
for fuzzy clustering, and their implementations in package `fclust` are
used.

## Usage

``` r
getSpatialGlobalInternalMetrics(
  labels,
  location,
  k = 6,
  metrics = c("PAS", "ELSA", "CHAOS"),
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

  The metrics to compute. See below for more details.

- ...:

  Optional arguments for
  [`PAS()`](https://roseyuan.github.io/poem/reference/PAS.md).

## Value

A named vector containing metric values. Possible metrics are:

- PAS:

  Proportion of abnormal spots (PAS score).

- ELSA:

  Entropy-based Local indicator of Spatial Association (ELSA score).

- CHAOS:

  Spatial Chaos Score.

- MPC:

  Modified partition coefficient

- PC:

  Partition coefficient

- PE:

  Partition entropy

## References

Yuan, Zhiyuan, et al., 2024; 10.1038/s41592-024-02215-8

Naimi, Babak, et al., 2019; 10.1016/j.spasta.2018.10.001

Wang, et al., 2022; 10.1016/j.ins.2022.11.010
