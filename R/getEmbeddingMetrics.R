#' getEmbeddingMetrics
#' 
#' Computes embedding-based metrics for the specified level.
#'
#' @param x A data.frame or matrix (with features as columns and items as rows) 
#'  from which the metrics will be computed.
#' @param labels A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param metrics The metrics to compute. See details. 
#' @param distance The distance metric to use (default euclidean).
#' @param level The level to calculate the metrics. Options include `"element"`, 
#' `"class"` and `"dataset"`.
#' @param ... Optional arguments. See details.
#'   
#' @return A data.frame of metrics.
#' @details
#' The allowed values for `metrics` depend on the value of `level`:
#'   - If `level = "element"`, the allowed `metrics` are: `"SW"`.
#'   - If `level = "class"`, the allowed `metrics` are: `"meanSW"`, `"minSW"`, `"pnSW"`, `"dbcv"`.
#'   - If `level = "dataset"`, the allowed `metrics` are: `"meanSW"`, 
#'   `"meanClassSW"`, `"pnSW"`, `"minClassSW"`, `"cdbw"`, `"cohesion"`, `"compactness"`, `"sep"`, `"dbcv"`.
#'   
#' The function(s) that the optional arguments `...` passed to depend on the 
#' value of `level`:
#'   - If `level = "element"`, optional arguments are passed to [stats::dist()].
#'   - If `level = "class"`, optional arguments are passed to [dbcv()].
#'   - If `level = "dataset"`, optional arguments are passed to [dbcv()] or [CDbw()].
#' @export
#' @examples
#' d1 <- mockData()
#' getEmbeddingMetrics(d1[,1:2], labels=d1$class, 
#' metrics=c("meanSW", "minSW", "pnSW", "dbcv"), level="class")
getEmbeddingMetrics <-function(x, labels, metrics=c("meanSW", "minSW", "pnSW", "dbcv"), 
                               distance="euclidean", level="class", ...){
  # Map level to the corresponding function
  level_functions <- list(
    "element" = getEmbeddingElementMetrics,
    "class" = getEmbeddingClassMetrics,
    "dataset" = getEmbeddingGlobalMetrics
  )
  .checkMetricsLevel(metrics, level, level_functions, use_default=TRUE, 
                     use_attribute=FALSE)
  # Collect all arguments into a list
  args <- list(x = x, labels = labels, metrics = metrics, distance = distance, ...)
  do.call(level_functions[[level]], args)
}

#' getEmbeddingElementMetrics
#' 
#' Computes element-level, embedding-based metrics.
#' @param metrics The metrics to compute. Currently, only the silhouette width
#'   is supported at the node-level.
#' @inheritParams getEmbeddingMetrics
#'   
#' @return A data.frame of metrics for each node/element of `x`.
#' 
#' @importFrom cluster silhouette
#' @importFrom stats dist
#' @examples
#' d1 <- mockData()
#' head(getEmbeddingElementMetrics(d1[,1:2], labels=d1$class))
getEmbeddingElementMetrics <-function(x, labels, metrics=c("SW"), 
                                      distance="euclidean", ...){
  stopifnot(is.atomic(labels) && (is.factor(labels) | is.integer(labels)))
  x <- as.matrix(x)
  stopifnot(mode(x) %in% c("numeric","integer"))
  stopifnot(length(labels)==nrow(x))
  stopifnot(length(match.arg(metrics))>0)
  .checkInvalidArgs(metrics, 
                  .get_allowed_args(getEmbeddingElementMetrics, "metrics"), 
                  "metrics")
  d <- dist(x, method=distance, ...)
  data.frame(row.names=row.names(x), class=labels,
             SW=cluster::silhouette(as.integer(labels), dist=d)[,3])
}

#' getEmbeddingClassMetrics
#' 
#' Computes class-level, embedding-based metrics.
#'
#' @inheritParams getEmbeddingMetrics
#' @param metrics The metrics to compute.
#'   
#' @return A data.frame of metrics for each node/element of `x`.
#' 
#' @importFrom stats aggregate
#' @examples
#' d1 <- mockData()
#' getEmbeddingClassMetrics(d1[,1:2], labels=d1$class)
getEmbeddingClassMetrics <-function(x, labels,
                                    metrics=c("meanSW", "minSW", "pnSW", "dbcv"),
                                    distance="euclidean", ...){
  stopifnot(is.atomic(labels) && (is.factor(labels) | is.integer(labels)))
  
  x <- as.matrix(x)
  stopifnot(mode(x) %in% c("numeric","integer"))
  stopifnot(length(labels)==nrow(x))
  .checkInvalidArgs(metrics, 
                    .get_allowed_args(getEmbeddingClassMetrics, "metrics"), "metrics")
  metrics <- match.arg(metrics, several.ok = TRUE)
  stopifnot(length(metrics)>0)
  res <- data.frame(class=sort(unique(labels)))
  nsw <- getEmbeddingElementMetrics(x, labels, distance=distance, metrics="SW")
  if("meanSW" %in% metrics)
    res$meanSW <- aggregate(nsw$SW, by=list(nsw$class), FUN=mean)[,2]
  if("minSW" %in% metrics)
    res$minSW <- aggregate(nsw$SW, by=list(nsw$class), FUN=min)[,2]
  if("pnSW" %in% metrics)
    res$pnSW <- aggregate(nsw$SW, by=list(nsw$class), FUN=function(x){
      sum(x<0)/length(x)
  })[,2]
  if("dbcv" %in% metrics){
    y <- as.integer(factor(labels, levels = sort(unique(labels))))
    res$dbcv <- dbcv(x, y, distance = distance, ...)$vcs}
  res
}

#' getEmbeddingGlobalMetrics
#' 
#' Computes dataset-level, embedding-based metrics.
#' 
#' @param metrics The metrics to compute.
#' @inheritParams getEmbeddingMetrics
#'   
#' @return A data.frame (with 1 row) of metrics.
#' @examples
#' d1 <- mockData()
#' getEmbeddingGlobalMetrics(d1[,1:2], labels=d1$class)
getEmbeddingGlobalMetrics <-function(x, labels,
                                     metrics=c("meanSW", "meanClassSW", "pnSW",
                                               "minClassSW", "cdbw", "cohesion",
                                               "compactness", "sep", "dbcv"),
                                     distance="euclidean", ...){
  stopifnot(is.atomic(labels) && (is.factor(labels) | is.integer(labels)))
  x <- as.matrix(x)
  stopifnot(mode(x) %in% c("numeric","integer"))
  stopifnot(length(labels)==nrow(x))
  .checkInvalidArgs(metrics, 
                    .get_allowed_args(getEmbeddingGlobalMetrics, "metrics"), 
                    "metrics")
  metrics <- match.arg(metrics, several.ok = TRUE)
  stopifnot(length(metrics)>0)
  res <- data.frame(class=sort(unique(labels)))
  nsw <- getEmbeddingElementMetrics(x, labels, distance=distance, metrics="SW")
  ag <- aggregate(nsw$SW, by=list(nsw$class), FUN=mean)[,2]
  res <- c(CDbw(x, as.integer(labels), ...),
    meanSW=mean(nsw$SW), pnSW=sum(nsw$SW<0)/nrow(nsw),
    meanClassSW=mean(ag), minClassSW=min(ag))[metrics]
  if("dbcv" %in% metrics){ # lazy computation of dbcv, because it's relatively slow
  y <- as.integer(factor(labels, levels = sort(unique(labels))))
  dbcv <- dbcv(x, y, distance = distance, ...)$dbcv
  res <- c(CDbw(x, as.integer(labels), ...),
           meanSW=mean(nsw$SW), pnSW=sum(nsw$SW<0)/nrow(nsw),
           meanClassSW=mean(ag), minClassSW=min(ag), dbcv=dbcv)[metrics]
  }else{
    res <- c(CDbw(x, as.integer(labels), ...),
             meanSW=mean(nsw$SW), pnSW=sum(nsw$SW<0)/nrow(nsw),
             meanClassSW=mean(ag), minClassSW=min(ag))[metrics]
  }
  as.data.frame(t(res))
}
