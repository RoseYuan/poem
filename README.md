# poem: POpulation-based Evaluation Metrics
<!-- badges: start -->
[![R-CMD-check](https://github.com/RoseYuan/poem/actions/workflows/check.yaml/badge.svg)](https://github.com/RoseYuan/poem/actions/workflows/check.yaml)
<!-- badges: end -->

## Overview
The `poem` package provides multiple metrics for evaluating subpopulation structure identification in a dataset. These include:

* Metrics for comparing two partitions[^1] of the same dataset, or metrics evaluating the alignment between a dataset’s embedding or graph representations with its partition.

* Metrics for comparing two fuzzy partitions[^2], or for comparing a hard partition with a fuzzy partition. This allows the evaluation of fuzzy partition results by assessing its agreement to a fuzzy or a hard ground-truth partition.

* Metrics tailored for domain detection in spatially-resolved omics data. These include especially external evaluation metrics (i.e. based on a comparison to ground truth labels), but also internal metrics. 

[^1]: A partition is a way of organizing the data points of a dataset into distinct, non-overlapping, and non-empty subsets. For example, a clustering is a partition. 

[^2]: In 'hard' partitions, each data point belongs to one and only one subset. However, clustering can also generate fuzzy partitions, in which data points can belong to multiple subsets with varying degrees (or probability) of membership. 

For a detailed introduction of the package, see the [online docs](https://roseyuan.github.io/poem/).



## Installation

`poem` is available through Bioconductor, and can be installed using the following commands:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("poem")
```

## Contact

In case you have any questions or suggestions to `poem`, please consider opening an
issue to the [GitHub repository](https://github.com/RoseYuan/poem).

## Citation

If you use *poem* or find it useful, please cite:

Luo S, Germain P, von Meyenn F, Robinson M (2025). “On metrics for subpopulation detection in single-cell and spatial omics data.” Nucleic Acids Research, 53(18). doi:10.1093/nar/gkaf921, https://doi.org/10.1093/nar/gkaf921.

### Bioconductor package
poem: POpulation-based Evaluation Metrics (v1.3.0). Bioconductor DOI: 10.18129/B9.bioc.poem, https://bioconductor.org/packages/poem/.


### BibTeX

```bibtex
@article{Luo2025Poem,
  title   = {On metrics for subpopulation detection in single-cell and spatial omics data},
  author  = {Luo, Siyuan and Germain, Pierre-Luc and von Meyenn, Ferdinand and Robinson, Mark D.},
  journal = {Nucleic Acids Research},
  volume  = {53},
  number  = {18},
  pages   = {gkaf921},
  year    = {2025},
  doi     = {10.1093/nar/gkaf921}
}

@manual{poemBioc,
  title        = {{poem}: POpulation-based Evaluation Metrics},
  author       = {Luo, Siyuan and Germain, Pierre-Luc},
  year         = {2025},
  note         = {R package version 1.3.0, Bioconductor},
  doi          = {10.18129/B9.bioc.poem},
  url          = {https://bioconductor.org/packages/poem/}
}
