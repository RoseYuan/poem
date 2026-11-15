# Generate mock multidimensional data

Generates mock multidimensional data of a given number of classes of
points, for testing.

## Usage

``` r
mockData(
  Ns = c(25, 15),
  classDiff = 2,
  Sds = 1,
  ndims = 2,
  spread = c(1, 2),
  rndFn = rnorm
)
```

## Arguments

- Ns:

  A vector of more than one positive integers specifying the number of
  elements of each class.

- classDiff:

  The distances between the classes. If there are more than 2 classes,
  this can be a `dist` object or a symmetric matrix of
  `length(Ns)-1`columns/rows where the lower triangle indicates the
  desired distances between classes.

- Sds:

  The standard deviation. Can either be a fixed value, a value per
  class, or a matrix of values for each class (rows) and dimension
  (column).

- ndims:

  The number of dimensions to generate (default 2).

- spread:

  The spread of the points. Can either be a fixed value, a value per
  class, or a matrix of values for each class (rows) and dimension
  (col).

- rndFn:

  The random function, by default `rnorm`, but should also work for
  `rlnorm` and similar.

## Value

A data.frame with coordinates and a `class` column.

## Examples

``` r
d <- mockData()
```
