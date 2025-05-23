#' Calculate Spatial External Metrics
#'
#' A generic function to calculate spatial external metrics. It can be applied 
#' to raw components (`true`, `pred`, `location`) or directly to a 
#' `SpatialExperiment` object.
#' @param object The main input. Can be a `SpatialExperiment` object or missing
#'   (when using `true`, `pred`, and `location` directly).
#' @param true When `object` is missing: a vector containing the labels of the 
#' true classes. Must be a vector of characters, integers, numerics, or a factor, 
#' but not a list. When `object` is a `SpatialExperiment` object: the column name
#' in `colData(object)` containing the true labels. 
#' @param pred When `object` is missing: a vector containing the labels of the
#' predicted clusters. Must be a vector of characters, integers, numerics, or a
#' factor, but not a list. When `object` is a `SpatialExperiment` object: the
#' column name in `colData(object)` containing the predicted labels.
#' @param metrics The metrics to compute. See details.
#' @param level The level to calculate the metrics. Options include `"element"`,
#' `"class"` and `"dataset"`.
#' @details
#' The allowed values for `metrics` depend on the value of `level`:
#'   - If `level = "element"`, the allowed `metrics` are: `"nsSPC"`, 
#'   `"NPC"`,`"SpatialSPC"`.
#'   - If `level = "class"`, the allowed `metrics` are: `"nsWH"`,
#'   `"nsAWH"`, `"nsWC"`,`"nsAWC"`.
#'   - If `level = "dataset"`, the allowed `metrics` are: `"nsRI"`,
#'   `"nsARI"`,`"nsWH"`,`"nsAWH"`, `"nsWC"`,`"nsAWC"`,
#'   `"nsAccuracy"`,`"SpatialRI"`,`"SpatialARI"`. 
#' @inheritParams getSpatialElementExternalMetrics
#' @param ... Additional arguments passed to specific methods.
#' @return A data.frame of metrics based on the specified input.
#' @importFrom SpatialExperiment spatialCoords SpatialExperiment
#' @importFrom SummarizedExperiment colData "colData<-" 
#' @examples
#' # Example with individual components
#' data(sp_toys)
#' data <- sp_toys
#' getSpatialExternalMetrics(true=data$label, pred=data$p1, 
#' location=data[,c("x", "y")], k=6, level="class")
#' 
#' # Example with SpatialExperiment object
#' se_object <- SpatialExperiment::SpatialExperiment(assays=matrix(NA, 
#'                                              ncol = nrow(data[,c("x", "y")]), 
#'                                              nrow = ncol(data[,c("x", "y")])), 
#'                                spatialCoords=as.matrix(data[,c("x", "y")]))
#' SummarizedExperiment::colData(se_object) <- 
#'                       cbind(SummarizedExperiment::colData(se_object), 
#'                             data.frame(true=data$label, pred=data$p1))
#' getSpatialExternalMetrics(object=se_object, true="true", pred="pred", k=6, 
#'                           level="class")
#'
#' @export
setGeneric("getSpatialExternalMetrics", signature="object",
           def=function(object=NULL, true, pred, location=NULL, 
           k=6, alpha=0.5, level="class",
           metrics=NULL, fuzzy_true=TRUE, fuzzy_pred=FALSE, ...) {
  level <- match.arg(level, c("dataset","class","element"))
  if(is.null(metrics))
    metrics <- switch(level,
      "dataset"=c("nsWH", "nsAWH", "nsWC", "nsAWC", "SpatialRI", "SpatialARI"),
      "class"=c("nsWH", "nsAWH", "nsWC", "nsAWC"),
      "element"=c("NPC", "nsSPC", "SpatialSPC"),
      stop("Unknown `level` specified.")
    )
  standardGeneric("getSpatialExternalMetrics")
})

#' @rdname getSpatialExternalMetrics
setMethod("getSpatialExternalMetrics", signature(object="missing"), 
          function(object, true, pred, location, k, 
                   alpha, level, metrics, fuzzy_true, fuzzy_pred, ...) {
                    # input validation
                    if (anyNA(true) | anyNA(pred)) stop("NA are not supported.")
                    if (is.character(true)) true <- as.factor(true)
                    if (is.character(pred)) pred <- as.factor(pred)
                    if (!is.atomic(true) || (!is.factor(true) && !is.integer(true)) ||
                        !is.atomic(pred) || (!is.factor(pred) && !is.integer(pred)) )
                      stop("true and pred must be vectors or factors but not lists.")
                    if(length(true) != length(pred)){
                      stop("The two input vectors should have the same length.")
                    }
                    # Extract the level functions
                    level_functions <- list(
                      "element" = getSpatialElementExternalMetrics,
                      "class" = getSpatialClassExternalMetrics,
                      "dataset" = getSpatialGlobalExternalMetrics
                    )
                    .checkMetricsLevel(metrics, level, level_functions, 
                                      use_default=TRUE, use_attribute=FALSE)
                    # Collect all arguments into a list
                    args <- list(true=true, pred=pred, location=location, 
                                k=k, alpha=alpha,
                                metrics=metrics, fuzzy_true=fuzzy_true, 
                                fuzzy_pred=fuzzy_pred, ...)
                    do.call(level_functions[[level]], args)
          })

#' @rdname getSpatialExternalMetrics
setMethod("getSpatialExternalMetrics", signature(object="SpatialExperiment"), 
          function(object, true, pred, k, alpha, level, metrics, 
          fuzzy_true, fuzzy_pred, ...) {
            if (!true %in% colnames(colData(object))) {
              stop("The column", true, "is not present.")
            }
            if (!pred %in% colnames(colData(object))) {
              stop("The column", pred, "is not present.")
            }
            # Extract true, pred, and location from the SpatialExperiment object
            true <- colData(object)[, true]
            pred <- colData(object)[, pred]
            location <- data.frame(spatialCoords(object))
            # Call the main function
            getSpatialExternalMetrics(true=true, pred=pred, location=location, 
                                       k=k, alpha=alpha, level=level, 
                                       metrics=metrics, fuzzy_true=fuzzy_true, 
                                       fuzzy_pred=fuzzy_pred, ...)
          })


#' Compute dataset-level external evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data at the dataset level. Options include a series of fuzzy pair-counting 
#' metrics and set matching-based accuracy.
#' @inheritParams getFuzzyPartitionMetrics
#' @inheritParams getFuzzyLabel
#' @keywords internal
#' @param true A vector containing the labels of the true classes. Must be a 
#'  vector of characters, integers, numerics, or a factor, but not a list.
#' @param pred A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param k The number of neighbors used when calculating the fuzzy 
#' class memberships for fuzzy metrics, or when calculating the weighted
#' accuracy.
#' @param metrics a vector of metric names to compute. 
#' @param fuzzy_true Logical; whether to compute fuzzy class memberships 
#' for `true`.
#' @param fuzzy_pred Logical; whether to compute fuzzy class memberships 
#' for `pred`.
#' @param ... Optional params for \code{\link{fuzzyPartitionMetrics}} or 
#'   \code{\link{findSpatialKNN}}.
#' @return A data.frame of metrics.
getSpatialGlobalExternalMetrics <- function(true, pred, location, 
                                            k=6, alpha=0.5,
                                            metrics=c("nsRI","nsARI",
                                                      "nsWH","nsAWH", 
                                                      "nsWC","nsAWC",
                                                      "nsAccuracy",
                                                      "SpatialRI","SpatialARI"), 
                                            fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                            lowMemory=NULL,
                                            ...){
  argfindSpatialKNN <- .checkEllipsisArgs(fnList=list(findSpatialKNN, 
                                              fuzzyPartitionMetrics,
                                              spatialARI), ...)[[1]] 
  argfuzzyPartitionMetrics <- .checkEllipsisArgs(fnList=list(findSpatialKNN, 
                                              fuzzyPartitionMetrics,
                                              spatialARI), ...)[[2]] 
  argspatialARI <- .checkEllipsisArgs(fnList=list(findSpatialKNN, 
                                              fuzzyPartitionMetrics,
                                              spatialARI), ...)[[3]]
  if(length(intersect(metrics, c("nsRI","nsARI","nsWH",
                                 "nsAWH", "nsWC","nsAWC")))>0){
    hardTrue <- true
    hardPred <- pred
    fuzzyTrue <- do.call(getFuzzyLabel,
                         c(argfindSpatialKNN, 
                           list(labels=hardTrue, location=location, k=k, 
                                alpha=alpha)))
    fuzzyPred <- do.call(getFuzzyLabel,
                         c(argfindSpatialKNN, 
                           list(labels=hardPred, location=location, k=k,
                                alpha=alpha)))
    res <- do.call(getFuzzyPartitionGlobalMetrics,
                   c(argfuzzyPartitionMetrics,
                     list(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue, 
                          hardPred=hardPred, fuzzyPred=fuzzyPred, 
                          fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                          lowMemory=lowMemory,
                          metrics=c("fuzzyRI", "fuzzyARI", "fuzzyWH", 
                                    "fuzzyAWH", "fuzzyWC", "fuzzyAWC"))))
  }else{res <- data.frame(matrix(nrow = 1, ncol = 0))}  
  if("nsAccuracy" %in% metrics){
    res$nsAccuracy <- do.call(nnWeightedAccuracy, 
                                   c(list(true=true, pred=pred, 
                                          location=location, k=k), 
                                     argfindSpatialKNN))
  }
  if("SpatialARI" %in% metrics | "SpatialRI" %in% metrics){
    tmp <- do.call(spatialARI, c(list(true=true, pred=pred, 
                                          location=location), 
                                          argspatialARI))
    res$SpatialRI <- tmp[1]
    res$SpatialARI <- tmp[2]
  }
  colnames(res) <- sub("fuzzy", "ns", colnames(res))
  return(res[,metrics, drop=FALSE])
}

#' Compute class-level external evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data at the class/cluster level. 
#' @inheritParams getSpatialGlobalExternalMetrics
#' @keywords internal
#' @param k The number of neighbors used when calculating the fuzzy 
#' class memberships for fuzzy metrics.
#' @param ... Optional params for \code{\link{fuzzyPartitionMetrics}} or 
#'   \code{\link{findSpatialKNN}}.
#' @return A data.frame of metrics.
getSpatialClassExternalMetrics <- function(true, pred, location, k=6, alpha=0.5,
                                           metrics=c("nsWH","nsAWH", 
                                                     "nsWC","nsAWC"), 
                                           fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                           lowMemory=NULL,
                                           ...){
  
  argfindSpatialKNN <- .checkEllipsisArgs(fnList=list(findSpatialKNN, 
                                              fuzzyPartitionMetrics), ...)[[1]] 
  argfuzzyPartitionMetrics <- .checkEllipsisArgs(fnList=list(findSpatialKNN, 
                                              fuzzyPartitionMetrics), ...)[[2]]
  hardTrue <- true
  hardPred <- pred
  fuzzyTrue <- do.call(getFuzzyLabel,
                       c(argfindSpatialKNN, 
                         list(labels=hardTrue, location=location, k=k,
                              alpha=alpha)))
  fuzzyPred <- do.call(getFuzzyLabel,
                       c(argfindSpatialKNN,
                         list(labels=hardPred, location=location, k=k,
                              alpha=alpha)))
  res <- do.call(getFuzzyPartitionClassMetrics,
                 c(argfuzzyPartitionMetrics,
                   list(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue,
                        hardPred=hardPred, fuzzyPred=fuzzyPred, 
                        fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                        lowMemory=lowMemory,
                      metrics=c("fuzzyWH","fuzzyAWH", "fuzzyWC", "fuzzyAWC"))))
  colnames(res) <- sub("fuzzy", "ns",colnames(res))
  return(res[,c(metrics, "class","cluster")])
}

#' getSpatialElementExternalMetrics
#'
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data at the element level.
#' @inheritParams getSpatialGlobalExternalMetrics
#' @param ... Optional params for [getFuzzyPartitionElementMetrics()] or 
#'   [findSpatialKNN()].
#' @keywords internal
#' @return A data.frame of metrics.
getSpatialElementExternalMetrics <- function(true, pred, location, 
                                             k=6, alpha=0.5,
                                         metrics=c("nsSPC", "NPC", "SpatialSPC"),
                                             fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                             ...){
  parsedArgs <- .checkEllipsisArgs(fnList=list(findSpatialKNN,
                                               getFuzzyPartitionElementMetrics,
                                               getNeighboringPairConcordance,
                                               spatialARI),
                                   ...)
  argfindSpatialKNN <- parsedArgs[[1]] 
  arggetFuzzyPartitionElementMetrics <- parsedArgs[[2]]
  arggetNeighboringPairConcordance <- parsedArgs[[3]]
  argspatialARI <- parsedArgs[[4]]
  if("nsSPC" %in% metrics){
  hardTrue <- true
  hardPred <- pred
  fuzzyTrue <- do.call(getFuzzyLabel,
                       c(argfindSpatialKNN, 
                         list(labels=hardTrue, location=location, k=k, 
                              alpha=alpha)))
  fuzzyPred <- do.call(getFuzzyLabel,
                       c(argfindSpatialKNN, 
                         list(labels=hardPred, location=location, k=k,
                              alpha=alpha)))
  nsSPC <- do.call(getFuzzyPartitionElementMetrics,
                        c(list(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue, 
                                hardPred=hardPred, fuzzyPred=fuzzyPred, 
                                fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                                metrics=c("nsSPC")),
                          arggetFuzzyPartitionElementMetrics))$fuzzySPC  
  }

  if("NPC" %in% metrics){
    NPC <- do.call(getNeighboringPairConcordance,
                          c(list(true=true, pred=pred, location=location, k=k),
                            arggetNeighboringPairConcordance))
  }
  if("SpatialSPC" %in% metrics){
    SpatialSPC <- do.call(spatialARI,
                          c(list(true=true, pred=pred, location=location, 
                          spotWise=TRUE), argspatialARI))
  }
  res <- as.data.frame(lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           nsSPC = nsSPC,
           NPC = NPC,
           SpatialSPC = SpatialSPC,
           stop("Unknown metric.")
           )})
    )
  return(res)
}
