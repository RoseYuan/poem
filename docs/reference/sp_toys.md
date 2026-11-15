# Toy examples of spatial data

Toy examples of spatial data. The code to generate the dataset is at
system.file('inst/scripts/', 'sp_toys.R', package='poem').

## Usage

``` r
data(sp_toys)
```

## Format

### `sp_toys`

A data frame with 240 rows and 11 columns, representing a 16 x 15 array
of spots:

- x, y:

  Coordinates of each spots.

- row, col:

  The row and column index of each spots.

- label:

  Ground truth labels. Either 1 or 2.

- p1-p6:

  Hypothetical predicted spatial clustering labels.
