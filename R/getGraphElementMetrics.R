#' getGraphElementMetrics
#' 
#' Computes a selection of supervised graph evaluation metrics using ground
#' truth class labels. The metrics are reported (as average) per node/element.
#'
#' @inheritParams getGraphMetrics
#' @param k The number of nearest neighbors to compute and/or use. Can be 
#'   omitted if `x` is a graph or list of nearest neighbors.
#' @param shared Logical; whether to use a shared nearest neighbor network 
#'   instead of a nearest neighbor network. Ignored if `x` is not an embedding 
#'   or dist object.
#' @param ... Optional arguments for [emb2knn()] or [emb2snn()].
#'
#' @return A data.frame of metrics for each node/element of `x`.
#' 
#' @rdname getGraphElementMetrics
setGeneric("getGraphElementMetrics", signature="x",
           def=function(x, labels, metrics=c("SI","NP","NCE"), 
                        directed=NULL, ...)
             standardGeneric("getGraphElementMetrics"))

attr(getGraphElementMetrics, "allowed_metrics") <- c("SI", "ISI", "NP","NCE")

#' @importFrom bluster neighborsToKNNGraph
#' @importFrom igraph adhesion cohesion set_vertex_attr
.getGraphElementMetricsFromKnn <- function(x, labels, metrics, directed=NULL, ...){
  .checkInputs(x,labels,checkNNcl=FALSE)
  x$nncl <- matrix(labels[x$index], nrow=nrow(x$index))
  labels <- as.factor(labels)
  if(isFALSE(directed)){
    # if the metrics are all to be computed in an undirected fashion, it's 
    # faster to start from an igraph:
    g <- bluster::neighborsToKNNGraph(x$index, directed=TRUE)
    return(.getGraphElementMetricsFromKnn(x, labels, metrics, directed=FALSE,...))
  }
  res <- as.data.frame(lapply(setNames(metrics,metrics), FUN=function(m){
    switch(m,
           SI=.simpsonIndex(x, labels, directed=directed),
           ISI=1/.simpsonIndex(x, labels, directed=directed),
           NP=.nPurity(x, labels, directed=directed),
           NCE=.nlog2Enrichment(x, labels, directed=directed),
           stop("The metric '",m,"' is either unknown or not available at the",
                " node-level.")
           )
  }))
  cbind(class=labels, res)
}

#' @rdname getGraphElementMetrics
setMethod("getGraphElementMetrics", signature="list",
          definition=function(x, labels, metrics, directed=NULL, k=NULL, ...){
            .checkInputs(x,labels,checkNNcl=FALSE)
            if(!is.null(k)){
              if(k>ncol(x$index))
                stop("The requested `k` is greater than the number of ",
                     "computed neighbors.")
              x <- lapply(x, FUN=function(x) x[,seq_len(k)])
            }
            .getGraphElementMetricsFromKnn(x, labels=labels, metrics=metrics,
                                    directed=directed, ...)
          })


.getGraphElementMetricsFromEmbedding <- function(x, labels, metrics, directed=NULL,
                                               k, shared=FALSE, ...){
  stopifnot(is.character(labels) || is.factor(labels) || is.integer(labels))
  labels <- as.factor(labels)
  stopifnot(length(labels)==nrow(x))
  if(is.data.frame(x)){
    stopifnot(all(vapply(x, FUN.VALUE=logical(1), FUN=is.numeric)))
    x <- as.matrix(x)
  }
  if(shared){
    g <- emb2snn(x, k=k, ...)
    res <- .getGraphElementMetricsFromGraph(g, labels=labels, metrics=metrics,
                                     directed=directed)
  }else{
    g <- emb2knn(x, k=k, ...)
    res <- .getGraphElementMetricsFromKnn(g, labels=labels, metrics=metrics,
                                   directed=directed)
  }
  res
}

#' @rdname getGraphElementMetrics
setMethod("getGraphElementMetrics", signature="data.frame",
          definition=.getGraphElementMetricsFromEmbedding)

#' @rdname getGraphElementMetrics
setMethod("getGraphElementMetrics", signature="matrix",
          definition=.getGraphElementMetricsFromEmbedding)


#' @importFrom igraph gorder set_vertex_attr
.getGraphElementMetricsFromGraph <- function(x, labels, 
                                           metrics=c("AMSP","PWC"),
                                           directed=NULL, ...){
  stopifnot(is.character(labels) || is.factor(labels))
  stopifnot(length(labels)==gorder(x))
  labels <- as.factor(labels)
  x <- set_vertex_attr(x, "class", value=labels)
  res <- as.data.frame(lapply(setNames(metrics,metrics), FUN=function(m){
    switch(m,
           SI=.simpsonIndex(x, labels, directed=directed),
           ISI=1/.simpsonIndex(x, labels, directed=directed),
           NP=.nPurity(x, labels, directed=directed),
           NCE=.nlog2Enrichment(x, labels, directed=directed),
           stop("The metric '",m,"' is either unknown or not available at the",
                " node-level.")
    )
  }))
  cbind(class=labels, res)
}

#' @rdname getGraphElementMetrics
setMethod("getGraphElementMetrics", signature="igraph",
          definition=function(x, labels, ...){
            stopifnot(is(x,"igraph"))
            .getGraphElementMetricsFromGraph(x, labels=labels, ...)
          })


.getGraphElementMetricsFromDist <- function(x, labels, shared=FALSE, k=10, ...){
  stopifnot(is.character(labels) || is.factor(labels))
  stopifnot(length(labels)==nrow(x))
  if(is.data.frame(x)){
    stopifnot(all(!vapply(x, FUN.VALUE=logical(1), FUN=is.numeric)))
    x <- as.matrix(x)
  }
  if(shared){
    g <- .dist2snn(x, k=k)
  }else{
    g <- .dist2knn(x, k=k)
  }
  getGraphElementMetrics(g, labels=labels, ...)
}

#' @rdname getGraphElementMetrics
setMethod("getGraphElementMetrics", signature="dist",
          definition=function(x, labels, ...){
            stopifnot(is(x,"dist"))
            .getGraphElementMetricsFromDist(x, labels=labels, ...)
          })