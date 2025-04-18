% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbcv.R
\name{.get_internal_objects}
\alias{.get_internal_objects}
\title{Get Internal Objects}
\usage{
.get_internal_objects(mutual_reach_dists, use_igraph_mst = TRUE)
}
\arguments{
\item{mutual_reach_dists}{Numeric matrix representing mutual reachability
distances.}

\item{use_igraph_mst}{Logical flag to use MST implementation
in igraph. Currently only mst from igraph is implemented.}
}
\value{
A list containing the indices of internal nodes and their edge
weights.
}
\description{
Computes the internal nodes and edges using Minimum Spanning
Tree.
}
\keyword{internal}
