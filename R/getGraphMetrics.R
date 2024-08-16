#' getGraphMetrics
#' 
#' Computes a selection of supervised graph evaluation metrics using ground
#' truth class labels. The metrics are reported (as average) per node/element.
#'
#' @param x Either an igraph object, a list of nearest neighbors (see details
#'   below), or a data.frame or matrix (with features as columns and items as 
#'   rows) from which nearest neighbors will be computed.
#' @param labels Either a factor or a character vector indicating the true class
#'   label of each element (i.e. row or vertex) of `x`.
#' @param metrics The metrics to compute. If omitted, main metrics will be 
#'   computed.
#' @param directed Logical; whether to compute the metrics in a directed fashion.
#'   If left to NULL, conventional choices will be made per metric (adhesion, 
#'   cohesion, PWC AMSP undirected, others directed).
#' @param k The number of nearest neighbors to compute and/or use. Can be 
#'   omitted if `x` is a graph or list of nearest neighbors.
#' @param BNPARAM A BiocNeighbors parameter object to compute kNNs. Ignored 
#'   unless the input is a matrix or data.frame. If omitted, the Annoy 
#'   approximation will be used if there are more than 500 elements.
#' @param shared Logical; whether to use a shared nearest neighbor network 
#'   instead of a nearest neighbor network. Ignored if `x` is not an embedding 
#'   or dist object.
#' @param ... Passed to other methods
#'
#' @return A data.frame of metrics for each node/element of `x`.
#' 
#' @details
#' Additional details...
#' 
#' 
#' @export
#' @examples
#' # generate random data and labels:
#' dat <- matrix(rnorm(300), nrow=30)
#' dat[1:10,1:3] <- dat[1:10,1:3]+2
#' labels <- rep(LETTERS[1:3], each=10)
#' # compute graph metrics:
#' head(getGraphMetrics(dat, labels, k=4))
setGeneric("getGraphMetrics", signature="x",
           def=function(x, labels, metrics=c("SI","NP","NCE"), 
                        directed=NULL, ...)
             standardGeneric("getGraphMetrics"))


#' @importFrom bluster neighborsToKNNGraph
#' @importFrom igraph adhesion cohesion set_vertex_attr
.getGraphMetricsFromKnn <- function(x, labels, metrics, directed=NULL, ...){
  .checkInputs(x,labels,checkNNcl=FALSE)
  x$nncl <- matrix(labels[x$index], nrow=nrow(x$index))
  labels <- as.factor(labels)
  if(isFALSE(directed)){
    # if the metrics are all to be computed in an undirected fashion, it's 
    # faster to start from an igraph:
    g <- bluster::neighborsToKNNGraph(x$index, directed=TRUE)
    return(.getGraphMetricsFromKnn(x, labels, metrics, directed=FALSE,...))
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

setMethod("getGraphMetrics", signature="list",
          definition=function(x, labels, metrics, directed=NULL, k=NULL, ...){
            .checkInputs(x,labels,checkNNcl=FALSE)
            if(!is.null(k)){
              if(k>ncol(x$index))
                stop("The requested `k` is greater than the number of ",
                     "computed neighbors.")
              x <- lapply(x, FUN=function(x) x[,seq_len(k)])
            }
            .getGraphMetricsFromKnn(x, labels=labels, metrics=metrics,
                                    directed=directed, ...)
          })


.getGraphMetricsFromEmbedding <- function(x, labels, metrics, directed=NULL,
                                               k, shared=FALSE, ...){
  stopifnot(is.character(labels) || is.factor(labels))
  stopifnot(length(labels)==nrow(x))
  if(is.data.frame(x)){
    stopifnot(all(vapply(x, FUN.VALUE=logical(1), FUN=is.numeric)))
    x <- as.matrix(x)
  }
  if(shared){
    g <- .emb2snn(x, k=k, ...)
  }else{
    g <- .emb2knn(x, k=k, ...)
  }
  .getGraphMetricsFromKnn(g, labels=labels, metrics=metrics, directed=directed)
}

setMethod("getGraphMetrics", signature="data.frame",
          definition=.getGraphMetricsFromEmbedding)
setMethod("getGraphMetrics", signature="matrix",
          definition=.getGraphMetricsFromEmbedding)


#' @importFrom igraph gorder set_vertex_attr
.getGraphMetricsFromGraph <- function(x, labels, 
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

setMethod("getGraphMetrics", signature="igraph",
          definition=function(x, labels, ...){
            stopifnot(is(x,"igraph"))
            .getGraphMetricsFromGraph(x, labels=labels, ...)
          })


.getGraphMetricsFromDist <- function(x, labels, shared=FALSE, k=10, ...){
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
  getGraphMetrics(g, labels=labels, ...)
}

setMethod("getGraphMetrics", signature="dist",
          definition=function(x, labels, ...){
            stopifnot(is(x,"dist"))
            .getGraphMetricsFromDist(x, labels=labels, ...)
          })