# The noisy moon dataset

A simple toy dataset consists of two interleaving half circles. The code
to generate the dataset is at system.file('inst/scripts/',
'noisy_moon.R', package='poem').

## Usage

``` r
data(noisy_moon)
```

## Format

### `noisy_moon`

A data frame with 100 rows and 5 columns:

- x, y:

  Coordinates of each observations.

- label:

  Ground truth labels. Either 1 or 2.

- kmeans_label:

  Predicted clustering labels using kmeans with 2 centers.

- hdbscan_label:

  Predicted clustering labels using hdbscan with `minPts = 5`.
