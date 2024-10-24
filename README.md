# poem: POpulation-based Evaluation Metrics
<!-- badges: start -->
[![R-CMD-check](https://github.com/RoseYuan/poem/actions/workflows/check.yaml/badge.svg)](https://github.com/RoseYuan/poem/actions/workflows/check.yaml)
<!-- badges: end -->

## Overview
The `poem` package provides multiple metrics for evaluating subpopulation structure identification in a dataset. These include:
- Metrics for comparing two partitions[^1] of the same dataset, or metrics evaluating the alignment between a dataset’s embedding or graph representations with its partition.
- Metrics for comparing two fuzzy partitions[^2], or for comparing a hard partition with a fuzzy partition. This allows the evaluation of fuzzy partition results by assessing its agreement to a fuzzy or a hard ground-truth partition.
- Metrics tailored for domain detection in spatially-resolved -omics data. These include especially external evaluation metrics (i.e. based on a comparison to ground truth labels), but also internal metrics. 

[^1]: A partition is a way of organizing the data points of a dataset into distinct, non-overlapping, and non-empty subsets. For example, a clustering is a partition. 

[^2]: In 'hard' partitions, each data point belongs to one and only one subset. However, clustering can also generate fuzzy partitions, in which data points can belong to multiple subsets with varying degrees (or probability) of membership. 

For a detailed introduction of the package, see the vignettes.



## Installation

You can install the development version of poem from
[GitHub](https://github.com/RoseYuan/poem) with:

``` r
# install.packages("devtools")
devtools::install_github("RoseYuan/poem")
```

## Contact

In case you have any questions or suggestions to `poem`, please consider opening an
issue to the [GitHub repository](https://github.com/RoseYuan/poem).
