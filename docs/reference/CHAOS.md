# Calculate CHAOS score

CHAOS score measures the clustering performance by calculating the mean
length of the graph edges in the 1-nearest neighbor (1NN) graph for each
cluster, averaged across clusters. Lower CHAOS score indicates better
spatial domain clustering performance.

## Usage

``` r
CHAOS(labels, location, BNPARAM = NULL)
```

## Arguments

- labels:

  Cluster labels.

- location:

  A numeric data matrix containing location information, where rows are
  points and columns are location dimensions.

- BNPARAM:

  BNPARAM object passed to
  [`findKNN`](https://rdrr.io/pkg/BiocNeighbors/man/findKNN.html)
  specifying the kNN approximation method to use. Defaults to exact for
  small datasets, and Annoy for larger ones.

## Value

A numeric value for CHAOS score.

## Examples

``` r
data(sp_toys)
data <- sp_toys
CHAOS(data$label, data[,c("x", "y")])
#> $CHAOS
#> [1] 0.2305866
#> 
#> $CHAOS_class
#>         2         1 
#> 0.2305866 0.2305866 
#> 
CHAOS(data$p1, data[,c("x", "y")])
#> $CHAOS
#> [1] 0.2305866
#> 
#> $CHAOS_class
#>         2         1 
#> 0.2305866 0.2305866 
#> 
CHAOS(data$p2, data[,c("x", "y")])
#> $CHAOS
#> [1] 0.2305866
#> 
#> $CHAOS_class
#>         2         1 
#> 0.2305866 0.2305866 
#> 
```
