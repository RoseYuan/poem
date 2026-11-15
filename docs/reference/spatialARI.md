# Spatially aware ARI from Yan, Yinqiao, et al. (2025).

Computes the spatial Rand Index and spatial ARI (Yan, Feng and Luo,
2025). Note that by default, the decay functions are different from
those of the original publication (see details for more information),
but the latter can be replicated with `original=TRUE`.

## Usage

``` r
spatialARI(
  true,
  pred,
  location,
  normCoords = TRUE,
  lambda = 0.8,
  fbeta = 4,
  hbeta = 1,
  spotWise = FALSE,
  nChunks = NULL,
  original = FALSE,
  f = function(x) {
lambda * exp(-x * fbeta)
 },
  h = function(x) {
lambda * (1 - exp(-x * hbeta))
 }
)
```

## Arguments

- true:

  A vector of true class labels

- pred:

  A vector of predicted clusters

- location:

  A matrix of spatial coordinates, with dimensions as columns

- normCoords:

  Logical; whether to normalize the coordinates to 0-1.

- lambda:

  The `alpha` used in the `f` and `h` functions (default 0.8) in Yan,
  Feng and Luo, 2025.

- fbeta, hbeta:

  Additional factors used in the exponential decay functions (see
  details). A higher value means a faster decay. These are ignored if
  `original=TRUE`.

- spotWise:

  Logical; whether to return the spot-wise spatial concordance (not
  adjusted for chance).

- nChunks:

  The number of processing chunks. If NULL, this will be determined
  automatically based on the size of the dataset, so as to remain below
  2GB RAM usage.

- original:

  Logical; whether to use the original h/f functions from Yan, Feng and
  Luo (default FALSE). If set to TRUE, the arguments `fbeta`, `hbeta`,
  `f` and `h` are ignored.

- f:

  The f function, which determines the positive contribution of pairs
  that are in different partitions in the reference, but grouped
  together in the clustering, based on the distance between mates.

- h:

  The h function, which determines the positive contribution of pairs
  that are in the same partition in the reference, but different ones in
  the clustering, based on the distance between mates.

## Value

A vector containing the spatial Rand Index (spRI) and spatial adjusted
Rand Index (spARI). Alternatively, if `spotWise=TRUE`, a vector of
spatial pair concordances for each spot.

## Details

This is a reimplementation of the method from the `spARI` package, made
more scalable (i.e. a bit slower but more memory-efficient) through
chunk-based processing, extensible to more than 2 dimensions, and with
some additional options. Note that by default, this will not produce the
same results as the original method: to do so, set `original=TRUE`. In
our exploration of the method and its behavior, we found the decay to be
too slow, and we therefore 1) do not square the distances, and 2)
introduced a beta parameter in each function which allows to scale it (a
higher beta parameter means a faster decay).

By default, chunking to keep RAM usage roughly below 2GB. Higher speed
can be achieved (at higher memory costs) for larger datasets by limiting
the number of chunks. The memory usage if done in a single chunk should
be roughly `4e-5*nrow(location)^2` Mb, and this scales down linearly
with the number of chunks.

## References

Yan, Feng and Luo, biorxiv 2025,
https://doi.org/10.1101/2025.03.25.645156

## Author

Pierre-Luc Germain

## Examples

``` r
data(sp_toys)
spatialARI(true=sp_toys$label, pred=sp_toys$p2, location = sp_toys[,1:2])
#>      spRI     spARI 
#> 0.8988816 0.7369664 
```
