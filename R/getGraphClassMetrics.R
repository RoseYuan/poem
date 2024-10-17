setOldClass("igraph")

#' getGraphClassMetrics
#' 
#' Computes a selection of supervised graph evaluation metrics using ground
#' truth class labels. The metrics are reported (as average) per class.
#'
#' @inheritParams getGraphMetrics
#' @param k The number of nearest neighbors to compute and/or use. Can be 
#'   omitted if `x` is a graph or list of nearest neighbors.
#' @param shared Logical; whether to use a shared nearest neighbor network 
#'   instead of a nearest neighbor network. Ignored if `x` is not an embedding 
#'   or dist object.
#' @param ... Optional arguments for [emb2knn()] or [emb2snn()].
#'
#' @return A data.frame of metrics for each class.
#' 
#' 
#' @rdname getGraphClassMetrics
setGeneric("getGraphClassMetrics", signature="x",
           def=function(x, labels, metrics=c("SI","NP","AMSP","PWC","NCE"), 
                        directed=NULL, ...)
             standardGeneric("getGraphClassMetrics"))

attr(getGraphClassMetrics, "allowed_metrics") <- c("SI", "ISI", "NP", "NCE", "AMSP", "PWC", "adhesion", "cohesion")

#' @importFrom bluster neighborsToKNNGraph
#' @importFrom igraph adhesion cohesion set_vertex_attr
.getGraphClassMetricsFromKnn <- function(x, labels, metrics, directed=NULL, ...){
  .checkInputs(x,labels,checkNNcl=FALSE)
  x$nncl <- matrix(labels[x$index], nrow=nrow(x$index))
  labels <- as.factor(labels)
  if(isFALSE(directed) || any(c("adhesion","cohesion","AMSP") %in% metrics)){
    # convert to graph
    g <- bluster::neighborsToKNNGraph(x$index, directed=FALSE)
    # if the metrics are all to be computed in an undirected fashion, it's 
    # faster to start from an igraph:
    if(isFALSE(directed))
      return(.getGraphClassMetricsFromGraph(g, labels, metrics,
                                            directed=FALSE, ...))
    g <- set_vertex_attr(g, "class", value=labels)
  }
  tt <- as.integer(table(labels))
  res <- as.data.frame(lapply(setNames(metrics,metrics), FUN=function(m){
    switch(m,
           SI=rowsum(.simpsonIndex(x, labels, directed=directed), labels)[,1]/tt,
           ISI=rowsum(1/.simpsonIndex(x, labels, directed=directed), labels)[,1]/tt,
           NP=rowsum(.nPurity(x, labels, directed=directed), labels)[,1]/tt,
           NCE=rowsum(.nlog2Enrichment(x, labels, directed=directed), labels)[,1]/tt,
           adhesion=.igraphFunPerClass(g, FUN=igraph::adhesion, directed=directed),
           cohesion=.igraphFunPerClass(g, FUN=igraph::cohesion, directed=directed),
           AMSP=.igraphFunPerClass(g, FUN=.adjMeanShortestPath, directed=directed),
           PWC=as.numeric(rowsum(as.integer(.nPurity(x,labels, directed=directed)<=0.5),
                                 labels)[,1]/tt),
           stop("Unknown metric ", m)
           )
  }))
  row.names(res) <- levels(labels)
  cbind(class=levels(labels),res)
}


#' @rdname getGraphClassMetrics
setMethod("getGraphClassMetrics", signature="list",
          definition=function(x, labels, metrics, directed=NULL, k=NULL, ...){
            .checkInputs(x,labels,checkNNcl=FALSE)
            if(!is.null(k)){
              if(k>ncol(x$index))
                stop("The requested `k` is greater than the number of ",
                     "computed neighbors.")
              x <- lapply(x, FUN=function(x) x[,seq_len(k)])
            }
            .getGraphClassMetricsFromKnn(x, labels=labels, metrics=metrics,
                                         directed=directed, ...)
          })


.getGraphClassMetricsFromEmbedding <- function(x, labels, metrics, directed=NULL,
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
    res <- .getGraphClassMetricsFromGraph(g, labels=labels, metrics=metrics,
                                          directed=directed)
  }else{
    g <- emb2knn(x, k=k, ...)
    res <- .getGraphClassMetricsFromKnn(g, labels=labels, metrics=metrics,
                                        directed=directed)
  }
  res
}

#' @rdname getGraphClassMetrics
setMethod("getGraphClassMetrics", signature="data.frame",
          definition=.getGraphClassMetricsFromEmbedding)

#' @rdname getGraphClassMetrics
setMethod("getGraphClassMetrics", signature="matrix",
          definition=.getGraphClassMetricsFromEmbedding)


#' @importFrom igraph gorder set_vertex_attr
.getGraphClassMetricsFromGraph <- function(x, labels, 
                                           metrics=c("AMSP","PWC"),
                                           directed=NULL, ...){
  stopifnot(is.character(labels) || is.factor(labels))
  stopifnot(length(labels)==gorder(x))
  labels <- as.factor(labels)
  x <- set_vertex_attr(x, "class", value=labels)
  tt <- as.integer(table(labels))
  res <- as.data.frame(lapply(setNames(metrics,metrics), FUN=function(m){
    switch(m,
           SI=rowsum(.simpsonIndex(x, labels, directed=directed), labels)[,1]/tt,
           ISI=rowsum(1/.simpsonIndex(x, labels, directed=directed), labels)[,1]/tt,
           NP=rowsum(.nPurity(x, labels, directed=directed), labels)[,1]/tt,
           NCE=rowsum(.nlog2Enrichment(x, labels, directed=directed), labels)[,1]/tt,
           adhesion=.igraphFunPerClass(x, FUN=igraph::adhesion, directed=directed),
           cohesion=.igraphFunPerClass(x, FUN=igraph::cohesion, directed=directed),
           AMSP=.igraphFunPerClass(x, FUN=.adjMeanShortestPath, directed=directed),
           PWC=as.numeric(rowsum(as.integer(.nPurity(x,labels, directed=directed)<=0.5),
                                 labels)[,1]/tt),
           stop("Unknown metric ", m)
    )
  }))
  row.names(res) <- levels(labels)
  cbind(class=levels(labels), res)
}

#' @rdname getGraphClassMetrics
setMethod("getGraphClassMetrics", signature="igraph",
          definition=function(x, labels, ...){
            stopifnot(is(x,"igraph"))
            .getGraphClassMetricsFromGraph(x, labels=labels, ...)
          })


.getGraphClassMetricsFromDist <- function(x, labels, shared=FALSE, k=10, ...){
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
  getGraphClassMetrics(g, labels=labels, ...)
}

#' @rdname getGraphClassMetrics
setMethod("getGraphClassMetrics", signature="dist",
          definition=function(x, labels, ...){
            stopifnot(is(x,"dist"))
            .getGraphClassMetricsFromDist(x, labels=labels, ...)
          })
