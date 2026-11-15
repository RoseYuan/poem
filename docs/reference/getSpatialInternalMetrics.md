# Compute internal metrics for spatial data

A generic function to compute a selection of internal clustering
evaluation metrics for spatial data. It can be applied to raw components
(`labels`, `location`) or directly to a `SpatialExperiment` object.

## Usage

``` r
getSpatialInternalMetrics(
  object = NULL,
  labels,
  location = NULL,
  k = 6,
  level = "class",
  metrics = c("CHAOS", "PAS", "ELSA"),
  ...
)

# S4 method for class 'missing'
getSpatialInternalMetrics(
  object = NULL,
  labels,
  location = NULL,
  k = 6,
  level = "class",
  metrics = c("CHAOS", "PAS", "ELSA"),
  ...
)

# S4 method for class 'ANY'
getSpatialInternalMetrics(
  object = NULL,
  labels,
  location = NULL,
  k = 6,
  level = "class",
  metrics = c("CHAOS", "PAS", "ELSA"),
  ...
)
```

## Arguments

- object:

  The main input. Can be a `SpatialExperiment` object or missing (when
  using `labels`, and `location` directly).

- labels:

  When `object` is missing: a vector containing the labels of the
  predicted clusters. Must be a vector of characters, integers,
  numerics, or a factor, but not a list. When `object` is a
  `SpatialExperiment` object: the column name in `colData(object)`
  containing the labels.

- location:

  A numerical matrix containing the location information, with rows as
  samples and columns as location dimensions.

- k:

  The size of the spatial neighborhood to look at for each spot. This is
  used for calculating PAS and ELSA scores.

- level:

  The level to calculate the metrics. Options include `"element"`,
  `"class"` and `"dataset"`.

- metrics:

  The metrics to compute. See details.

- ...:

  Optional params for
  [`PAS()`](https://roseyuan.github.io/poem/reference/PAS.md).

## Value

A data.frame of metrics.

## Details

The allowed values for `metrics` depend on the value of `level`:

- If `level = "element"`, the allowed `metrics` are: `"PAS"`, `"ELSA"`.

- If `level = "class"`, the allowed `metrics` are: `"CHAOS"`, `"PAS"`,
  `"ELSA"`.

- If `level = "dataset"`, the allowed `metrics` are:

  - `"PAS"`: Proportion of abnormal spots (PAS score)

  - `"ELSA"`: Entropy-based Local indicator of Spatial Association (ELSA
    score)

  - `"CHAOS"`: Spatial Chaos Score.

  - `"MPC"`: Modified partition coefficient

  - `"PC"`: Partition coefficient

  - `"PE"`: Partition entropy

## Examples

``` r
# Example with individual components
data(sp_toys)
data <- sp_toys
getSpatialInternalMetrics(labels=data$label, location=data[,c("x", "y")], 
                          k=6, level="class")
#> the specified variable is considered as categorical...
#>   class         PAS    ELSA.Ea   ELSA.Ec       ELSA     CHAOS
#> 1     1 0.008333333 0.07106927 0.1940084 0.06040663 0.2305866
#> 2     2 0.008333333 0.06918805 0.1944579 0.05840334 0.2305866

# Example with SpatialExperiment object
se_object <- SpatialExperiment::SpatialExperiment(assays=matrix(NA, 
                                             ncol = nrow(data[,c("x", "y")]), 
                                             nrow = ncol(data[,c("x", "y")])), 
                               spatialCoords=as.matrix(data[,c("x", "y")]))
SummarizedExperiment::colData(se_object) <- 
cbind(SummarizedExperiment::colData(se_object), 
      data.frame(label=data$label))
getSpatialInternalMetrics(object=se_object, labels="label", k=6, 
                          level="class")
#> the specified variable is considered as categorical...
#>   class         PAS    ELSA.Ea   ELSA.Ec       ELSA     CHAOS
#> 1     1 0.008333333 0.07106927 0.1940084 0.06040663 0.2305866
#> 2     2 0.008333333 0.06918805 0.1944579 0.05840334 0.2305866
```
