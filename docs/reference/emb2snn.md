# Computes shared nearest neighbors from embedding

computes shared nearest neighbors from embedding.

## Usage

``` r
emb2snn(x, k, type = "rank", BNPARAM = NULL)
```

## Arguments

- x:

  A numeric matrix (with features as columns and items as rows) from
  which nearest neighbors will be computed.

- k:

  The number of nearest neighbors.

- type:

  A string specifying the type of weighting scheme to use for shared
  neighbors. Possible choices include "rank", "number", and "jaccard".
  See `type` in
  [`bluster::neighborsToSNNGraph()`](https://rdrr.io/pkg/bluster/man/makeSNNGraph.html)
  for details.

- BNPARAM:

  A BiocNeighbors parameter object to compute kNNs. Ignored unless the
  input is a matrix or data.frame. If omitted, the Annoy approximation
  will be used if there are more than 500 elements.

## Value

An igraph object.

## Examples

``` r
d1 <- mockData()
emb2snn(as.matrix(d1[,seq_len(2)]),k=5)
#> IGRAPH 13beedc U-W- 40 257 -- 
#> + attr: weight (e/n)
#> + edges from 13beedc:
#>  [1]  1-- 3  3-- 4  2-- 4  1-- 5  3-- 6  1-- 6  4-- 6  1-- 7  5-- 7  6-- 7
#> [11]  3-- 7  4-- 7  4-- 8  3-- 8  7-- 8  5-- 9  7-- 9  8-- 9  2--10  3--10
#> [21]  6--10  8--11  9--11  7--11  3--11  4--11  1--12  3--12  6--12  7--12
#> [31]  4--12  8--12 11--12  2--13 10--13  4--13  4--14  8--14 11--14  9--14
#> [41]  2--15 10--15 13--15  4--15  3--15  8--15  2--16  4--16 13--16 15--16
#> [51] 10--16  3--16  8--16  3--17  6--17 10--17 12--17  4--17  1--17  7--17
#> [61] 17--18  1--18  6--18 12--18  3--18  7--18  7--19  9--19 11--19  8--19
#> [71]  5--19  3--19  4--19 12--19 14--19  1--19  2--20 10--20 13--20 15--20
#> + ... omitted several edges
```
