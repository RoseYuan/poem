# Package index

## Main functions

All metrics can be retrieved via the following six main wrapper
functions.

- [`getEmbeddingMetrics()`](https://roseyuan.github.io/poem/reference/getEmbeddingMetrics.md)
  : Compute embedding-based metrics
- [`getGraphMetrics()`](https://roseyuan.github.io/poem/reference/getGraphMetrics.md)
  : Compute graph-based metrics
- [`getPartitionMetrics()`](https://roseyuan.github.io/poem/reference/getPartitionMetrics.md)
  : Compute partition-based metrics
- [`getFuzzyPartitionMetrics()`](https://roseyuan.github.io/poem/reference/getFuzzyPartitionMetrics.md)
  : Compute external metrics for fuzzy clusterings
- [`getSpatialExternalMetrics()`](https://roseyuan.github.io/poem/reference/getSpatialExternalMetrics.md)
  : Calculate Spatial External Metrics
- [`getSpatialInternalMetrics()`](https://roseyuan.github.io/poem/reference/getSpatialInternalMetrics.md)
  : Compute internal metrics for spatial data

## Functions for individual metrics

These are some useful functions for calculating individual metrics.

### Embedding metrics

- [`CDbw()`](https://roseyuan.github.io/poem/reference/CDbw.md) :
  Calculate CDbw index
- [`dbcv()`](https://roseyuan.github.io/poem/reference/dbcv.md) :
  Calculate DBCV Metric
- [`silhouetteWidths()`](https://roseyuan.github.io/poem/reference/silhouetteWidths.md)
  : silhouetteWidths

### Fuzzy metrics

- [`fuzzyHardMetrics()`](https://roseyuan.github.io/poem/reference/fuzzyHardMetrics.md)
  : Compute fuzzy-hard versions of pair-sorting partition metrics
- [`fuzzyPartitionMetrics()`](https://roseyuan.github.io/poem/reference/fuzzyPartitionMetrics.md)
  : Compute fuzzy-fuzzy versions of pair-sorting partition metrics

### Internal metrics for spatial domain detection

- [`PAS()`](https://roseyuan.github.io/poem/reference/PAS.md) :
  Calculate PAS score
- [`CHAOS()`](https://roseyuan.github.io/poem/reference/CHAOS.md) :
  Calculate CHAOS score
- [`ELSA()`](https://roseyuan.github.io/poem/reference/ELSA.md) :
  Calculate ELSA scores

### External metrics for spatial domain detection

- [`spatialARI()`](https://roseyuan.github.io/poem/reference/spatialARI.md)
  : Spatially aware ARI from Yan, Yinqiao, et al. (2025).

### Element-wise external metrics for spatial domain detection, particularly useful for visualizations.

- [`fuzzyHardSpotConcordance()`](https://roseyuan.github.io/poem/reference/fuzzyHardSpotConcordance.md)
  : Per-element maximal concordance between a hard and a fuzzy partition
- [`getPairConcordance()`](https://roseyuan.github.io/poem/reference/getPairConcordance.md)
  : Per-element pair concordance score
- [`getNeighboringPairConcordance()`](https://roseyuan.github.io/poem/reference/getNeighboringPairConcordance.md)
  : Per-element local concordance between a clustering and a ground
  truth

## Datasets

The package provides several datasets for demonstration purposes.

- [`noisy_moon`](https://roseyuan.github.io/poem/reference/noisy_moon.md)
  : The noisy moon dataset
- [`sp_toys`](https://roseyuan.github.io/poem/reference/sp_toys.md) :
  Toy examples of spatial data
- [`toyExamples`](https://roseyuan.github.io/poem/reference/toyExamples.md)
  : Toy embedding examples
- [`metric_info`](https://roseyuan.github.io/poem/reference/metric_info.md)
  : Metrics Information

## Other utility functions

- [`mockData()`](https://roseyuan.github.io/poem/reference/mockData.md)
  : Generate mock multidimensional data
- [`matchSets()`](https://roseyuan.github.io/poem/reference/matchSets.md)
  : Match two partitions using Hungarian algorithm
- [`knnComposition()`](https://roseyuan.github.io/poem/reference/knnComposition.md)
  : Compute neighborhood composition
- [`getFuzzyLabel()`](https://roseyuan.github.io/poem/reference/getFuzzyLabel.md)
  : Get fuzzy representation of labels
- [`findSpatialKNN()`](https://roseyuan.github.io/poem/reference/findSpatialKNN.md)
  : Find the k nearest spatial neighbors
- [`emb2knn()`](https://roseyuan.github.io/poem/reference/emb2knn.md) :
  Computes k nearest neighbors from embedding
- [`emb2snn()`](https://roseyuan.github.io/poem/reference/emb2snn.md) :
  Computes shared nearest neighbors from embedding
