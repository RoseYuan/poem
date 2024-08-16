#' getEmbeddingMetrics
#' 
#' Computes element-level, embedding-based metrics.
#'
#' @param x A data.frame or matrix (with features as columns and items as rows) 
#'  from which the metrics will be computed.
#' @param labels A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param metrics The metrics to compute. Currently, only the silhouette width
#'   is supported at the node-level.
#' @param distance The distance metric to use (default euclidean).
#'   
#' @return A data.frame of metrics for each node/element of `x`.
#' @details
#' Additional details...
#' 
#' @importFrom cluster silhouette
#' @export
#' @examples
#' d1 <- mockData()
#' head(getEmbeddingMetrics(d1[,1:2], labels=d1$class))
getEmbeddingMetrics <-function(x, labels, metrics=c("SW"), distance="euclidean",
                                ...){
  stopifnot(is.atomic(labels) && (is.factor(labels) | is.integer(labels)))
  x <- as.matrix(x)
  stopifnot(mode(x) %in% c("numeric","integer"))
  stopifnot(length(labels)==nrow(x))
  stopifnot(length(match.arg(metrics))>0)
  d <- dist(x, method=distance, ...)
  data.frame(row.names=row.names(x), class=labels,
             SW=cluster::silhouette(as.integer(labels), dist=d)[,3])
}

#' getClassEmbeddingMetrics
#' 
#' Computes class-level, embedding-based metrics.
#'
#' @param x A data.frame or matrix (with features as columns and items as rows) 
#'  from which the metrics will be computed.
#' @param labels A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param metrics The metrics to compute.
#' @param distance The distance metric to use (default euclidean).
#'   
#' @return A data.frame of metrics for each node/element of `x`.
#' @details
#' Additional details...
#' 
#' @export
#' @examples
#' d1 <- mockData()
#' getClassEmbeddingMetrics(d1[,1:2], labels=d1$class)
getClassEmbeddingMetrics <-function(x, labels,
                                    metrics=c("meanSW", "minSW", "pnSW"),
                                    distance="euclidean", ...){
  stopifnot(is.atomic(labels) && (is.factor(labels) | is.integer(labels)))
  metrics <- match.arg(metrics, several.ok = TRUE)
  x <- as.matrix(x)
  stopifnot(mode(x) %in% c("numeric","integer"))
  stopifnot(length(labels)==nrow(x))
  stopifnot(length(match.arg(metrics))>0)
  res <- data.frame(class=sort(unique(labels)))
  nsw <- getEmbeddingMetrics(x, labels, distance=distance, metrics="SW")
  if("meanSW" %in% metrics)
    res$meanSW <- aggregate(nsw$SW, by=list(nsw$class), FUN=mean)[,2]
  if("minSW" %in% metrics)
    res$minSW <- aggregate(nsw$SW, by=list(nsw$class), FUN=min)[,2]
  if("pnSW" %in% metrics)
    res$pnSW <- aggregate(nsw$SW, by=list(nsw$class), FUN=function(x){
      sum(x<0)/length(x)
  })[,2]
  res
}

#' getGlobalEmbeddingMetrics
#' 
#' Computes global, embedding-based metrics.
#'
#' @param x A data.frame or matrix (with features as columns and items as rows) 
#'  from which the metrics will be computed.
#' @param labels A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param metrics The metrics to compute.
#' @param distance The distance metric to use (default euclidean).
#'   
#' @return A data.frame (with 1 row) of metrics.
#' @details
#' Additional details...
#' 
#' @export
#' @examples
#' d1 <- mockData()
#' getGlobalEmbeddingMetrics(d1[,1:2], labels=d1$class)
getGlobalEmbeddingMetrics <-function(x, labels,
                                     metrics=c("meanSW", "meanClassSW", "pnSW",
                                               "minClassSW", "cdbw", "cohesion",
                                               "compactness", "sep"),
                                     distance="euclidean", ...){
  stopifnot(is.atomic(labels) && (is.factor(labels) | is.integer(labels)))
  metrics <- match.arg(metrics, several.ok = TRUE)
  x <- as.matrix(x)
  stopifnot(mode(x) %in% c("numeric","integer"))
  stopifnot(length(labels)==nrow(x))
  stopifnot(length(match.arg(metrics))>0)
  res <- data.frame(class=sort(unique(labels)))
  nsw <- getEmbeddingMetrics(x, labels, distance=distance, metrics="SW")
  ag <- aggregate(nsw$SW, by=list(nsw$class), FUN=mean)[,2]
  res <- c(CDbw(x, as.integer(labels), ...),
    meanSW=mean(nsw$SW), pnSW=sum(nsw$SW<0)/nrow(nsw),
    meanClassSW=mean(ag), minClassSW=min(ag))[metrics]
  as.data.frame(t(res))
}