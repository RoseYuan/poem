# Calculate Spatial External Metrics

A generic function to calculate spatial external metrics. It can be
applied to raw components (`true`, `pred`, `location`) or directly to a
`SpatialExperiment` object.

## Usage

``` r
getSpatialExternalMetrics(
  object = NULL,
  true,
  pred,
  location = NULL,
  k = 6,
  alpha = 0.5,
  level = "class",
  metrics = NULL,
  fuzzy_true = TRUE,
  fuzzy_pred = FALSE,
  ...
)

# S4 method for class 'missing'
getSpatialExternalMetrics(
  object = NULL,
  true,
  pred,
  location = NULL,
  k = 6,
  alpha = 0.5,
  level = "class",
  metrics = NULL,
  fuzzy_true = TRUE,
  fuzzy_pred = FALSE,
  ...
)

# S4 method for class 'ANY'
getSpatialExternalMetrics(
  object = NULL,
  true,
  pred,
  location = NULL,
  k = 6,
  alpha = 0.5,
  level = "class",
  metrics = NULL,
  fuzzy_true = TRUE,
  fuzzy_pred = FALSE,
  ...
)
```

## Arguments

- object:

  The main input. Can be a `SpatialExperiment` object or missing (when
  using `true`, `pred`, and `location` directly).

- true:

  When `object` is missing: a vector containing the labels of the true
  classes. Must be a vector of characters, integers, numerics, or a
  factor, but not a list. When `object` is a `SpatialExperiment` object:
  the column name in `colData(object)` containing the true labels.

- pred:

  When `object` is missing: a vector containing the labels of the
  predicted clusters. Must be a vector of characters, integers,
  numerics, or a factor, but not a list. When `object` is a
  `SpatialExperiment` object: the column name in `colData(object)`
  containing the predicted labels.

- location:

  A matrix or data.frame of coordinates

- k:

  The number of neighbors used when calculating the fuzzy class
  memberships for fuzzy metrics, or when calculating the weighted
  accuracy.

- alpha:

  The parameter to control to what extend the spot itself contribute to
  the class composition calculation. `"equal"` means it is weighted the
  same as other neighbors. A numeric value between `0` and `1` means the
  weight of the frequency contribution for the spot itself, and the
  frequency contribution for its knn is then `1-alpha`. By default
  `0.5`.

- level:

  The level to calculate the metrics. Options include `"element"`,
  `"class"` and `"dataset"`.

- metrics:

  The metrics to compute. See details.

- fuzzy_true:

  Logical; whether to compute fuzzy class memberships for `true`.

- fuzzy_pred:

  Logical; whether to compute fuzzy class memberships for `pred`.

- ...:

  Additional arguments passed to specific methods.

## Value

A data.frame of metrics based on the specified input.

## Details

The allowed values for `metrics` depend on the value of `level`:

- If `level = "element"`, the allowed `metrics` are: `"nsSPC"`,
  `"NPC"`,`"SpatialSPC"`.

- If `level = "class"`, the allowed `metrics` are: `"nsWH"`, `"nsAWH"`,
  `"nsWC"`,`"nsAWC"`.

- If `level = "dataset"`, the allowed `metrics` are: `"nsRI"`,
  `"nsARI"`,`"nsWH"`,`"nsAWH"`, `"nsWC"`,`"nsAWC"`,
  `"nsAccuracy"`,`"SpatialRI"`,`"SpatialARI"`.

## Examples

``` r
# Example with individual components
data(sp_toys)
data <- sp_toys
getSpatialExternalMetrics(true=data$label, pred=data$p1, 
location=data[,c("x", "y")], k=6, level="class")
#> Comparing between a fuzzy truth and a hard prediction...
#> Standard error of the mean NDC across permutations:0.00118
#>        nsWH     nsAWH      nsWC     nsAWC class cluster
#> 1        NA        NA 0.8078698 0.5929239     1      NA
#> 2        NA        NA 1.0000000 1.0000000     2      NA
#> 3 1.0000000 1.0000000        NA        NA    NA       1
#> 4 0.8323893 0.6513298        NA        NA    NA       2

# Example with SpatialExperiment object
se_object <- SpatialExperiment::SpatialExperiment(assays=matrix(NA, 
                                             ncol = nrow(data[,c("x", "y")]), 
                                             nrow = ncol(data[,c("x", "y")])), 
                               spatialCoords=as.matrix(data[,c("x", "y")]))
SummarizedExperiment::colData(se_object) <- 
                      cbind(SummarizedExperiment::colData(se_object), 
                            data.frame(true=data$label, pred=data$p1))
getSpatialExternalMetrics(object=se_object, true="true", pred="pred", k=6, 
                          level="class")
#> Comparing between a fuzzy truth and a hard prediction...
#> Standard error of the mean NDC across permutations:0.000813
#>        nsWH     nsAWH      nsWC     nsAWC class cluster
#> 1        NA        NA 0.8078698 0.5921384     1      NA
#> 2        NA        NA 1.0000000 1.0000000     2      NA
#> 3 1.0000000 1.0000000        NA        NA    NA       1
#> 4 0.8323893 0.6500703        NA        NA    NA       2
```
