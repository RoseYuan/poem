#' Compute graph-based metrics
#' 
#' Computes a selection of graph evaluation metrics using class labels.
#'
#' @param x Either an igraph object, a list of nearest neighbors (see details
#'   below), or a data.frame or matrix (with features as columns and items as 
#'   rows) from which nearest neighbors will be computed.
#' @param labels Either a factor or a character vector indicating the true class
#'   label of each element (i.e. row or vertex) of `x`.
#' @param metrics The metrics to compute. See details.
#' @param directed Logical; whether to compute the metrics in a directed 
#' fashion. If left to NULL, conventional choices will be made per metric 
#' (adhesion, cohesion, PWC AMSP undirected, others directed).
#' @param level The level to calculate the metrics. Options include 
#' `"element"`, `"class"` and `"dataset"`.
#' @param k The number of nearest neighbors to compute and/or use. Can be 
#'   omitted if `x` is a graph or list of nearest neighbors.
#' @param shared Logical; whether to use a shared nearest neighbor network 
#'   instead of a nearest neighbor network. Ignored if `x` is not an embedding 
#'   or dist object.
#' @param ... Optional arguments for [emb2knn()] or [emb2snn()].
#' @return A data.frame of metrics.
#' @details
#' The allowed values for `metrics` depend on the value of `level`:
#'   - If `level = "element"`, the allowed `metrics` are: `"SI"`,`"ISI"`,`"NP"`,
#'   `"NCE"` (see below for details).
#'   - If `level = "class"`, the allowed `metrics` are: 
#'      - `"SI"`: Simpson’s Index.
#'      - `"ISI"`: Inverse Simpson’s Index
#'      - `"NP"`: Neighborhood Purity
#'      - `"AMSP"`: Adjusted Mean Shortest Path
#'      - `"PWC"`: Proportion of Weakly Connected 
#'      - `"NCE"`: Neighborhood Class Enrichment
#'      - `"adhesion"`: adhesion of a graph, is the minumum number of nodes 
#'      that must be removed to split a graph.
#'      - `"cohesion"`: cohesion of a graph, is the minumum number of edges 
#'      that must be removed to split a graph.
#'   - If `level = "dataset"`, the allowed `metrics` are: `"SI"`,`"ISI"`,
#'   `"NP"`,`"AMSP"`,`"PWC"`,`"NCE"`, `"adhesion"`,`"cohesion"`.
#' @export
#' @examples
#' d1 <- mockData()
#' getGraphMetrics(d1[,seq_len(2)], labels=d1$class, level="class")
getGraphMetrics <-function(x, labels, metrics=NULL, 
                           directed=NULL, k=10, shared=FALSE, 
                           level="class", ...){
  # Map level to the corresponding function
  level <- match.arg(level, c("dataset","class","element"))
  level_functions <- list(
    "element" = getGraphElementMetrics,
    "class" = getGraphClassMetrics,
    "dataset" = getGraphGlobalMetrics
  )
  if(is.null(metrics))
    metrics <- switch(level,
                      "dataset"=c("SI","NP","AMSP","PWC","NCE"),
                      "class"=c("SI","NP","AMSP","PWC","NCE"),
                      "element"=c("SI","NP"),
                      stop("Unknown `level` specified.")
    )
  
  .checkMetricsLevel(metrics, level, level_functions, use_default=FALSE, 
                     use_attribute=TRUE, attr_name="allowed_metrics")
  # Collect all arguments into a list
  args <- list(x = x, labels = labels, metrics = metrics, directed=directed, 
               k=k, shared=shared, ...)
  do.call(level_functions[[level]], args)
}

getGraphGlobalMetrics <- function(x, labels, 
                                  metrics=c("SI","NP","AMSP","PWC","NCE"),
                                  directed=NULL, k=10, shared=FALSE,...){
  .class2global(getGraphClassMetrics(x=x, labels=labels, metrics=metrics,
                                     directed=directed, k=k, 
                                     shared=shared, ...))
}

attr(getGraphGlobalMetrics, "allowed_metrics") <- c("SI", "ISI", "NP", "NCE", 
                                                    "AMSP", "PWC", "adhesion", 
                                                    "cohesion")
