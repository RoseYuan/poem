# Calculate CDbw index

Computes the CDbw-index (Halkidi and Vazirgiannis 2008; Halkidi,
Vazirgiannis and Hennig, 2015). This function is directly copied from
the `fpc` CRAN package and was written by Christian Hennig. It is
included here to reduce the package dependencies (since `fpc` has a few
not-so-light dependencies that aren't required here).

## Usage

``` r
CDbw(
  x,
  labels,
  r = 10,
  s = seq(0.1, 0.8, by = 0.1),
  clusterstdev = TRUE,
  trace = FALSE
)
```

## Arguments

- x:

  Something that can be coerced into a numerical matrix, with elements
  as rows.

- labels:

  A vector of integers with length `=nrow(x)` indicating the cluster for
  each observation.

- r:

  Number of cluster border representatives.

- s:

  Vector of shrinking factors.

- clusterstdev:

  Logical. If `TRUE`, the neighborhood radius for intra-cluster density
  is the within-cluster estimated squared distance from the mean of the
  cluster; otherwise it is the average of these over all clusters.

- trace:

  Logical; whether to print processing info.

## Value

A vector with the following values (see refs for details):

- cdbw:

  value of CDbw index (the higher the better).

- cohesion:

  cohesion.

- compactness:

  compactness.

- sep:

  separation.

## References

Halkidi, M. and Vazirgiannis, M. (2008) A density-based cluster validity
approach using multi-representatives. Pattern Recognition Letters 29,
773-786.

Halkidi, M., Vazirgiannis, M. and Hennig, C. (2015) Method-independent
indices for cluster validation. In C. Hennig, M. Meila, F. Murtagh, R.
Rocci (eds.) Handbook of Cluster Analysis, CRC Press/Taylor `&` Francis,
Boca Raton.

## Author

Christian Hennig

## Examples

``` r
d1 <- mockData()
CDbw(d1[,seq_len(2)], d1[,3])
#>        cdbw    cohesion compactness         sep 
#>  0.03805049  0.23717624  0.24202068  0.66288255 
```
