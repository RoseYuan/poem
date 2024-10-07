#' getSpatialExternalMetrics
#' 
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data.
#' @inheritParams getSpatialElementExternalMetrics
#' @param metrics The metrics to compute. See details.
#' @param level The level to calculate the metrics. Options include `"element"`,
#' `"class"` and `"global"`.
#' @return A data.frame of metrics.
#' @export
#' @details
#' The allowed values for `metrics` depend on the value of `level`:
#'   - If `level = "element"`, the allowed `metrics` are: `"spotAgreement"`.
#'   - If `level = "class"`, the allowed `metrics` are: `"SpatialWH"`,`"SpatialAWH"`, `"SpatialWC"`,`"SpatialAWC"`.
#'   - If `level = "global"`, the allowed `metrics` are: `"SpatialRI"`,`"SpatialARI"`,`"SpatialWH"`,`"SpatialAWH"`, `"SpatialWC"`,`"SpatialAWC"`,`"SpatialAccuracy"`. 
#' @examples
#' data <- sp_toys
#' getSpatialExternalMetrics(data$label, data$p1, data[,c("x", "y")], k=6, level="class")
getSpatialExternalMetrics <- function(true, pred, location, k=6, alpha=0.5, level="class",
                                      metrics=c("SpatialWH","SpatialAWH", 
                                                "SpatialWC","SpatialAWC"),
                                      fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                      ...){
  # Map level to the corresponding function
  level_functions <- list(
    "element" = getSpatialElementExternalMetrics,
    "class" = getSpatialClassExternalMetrics,
    "global" = getSpatialGlobalExternalMetrics
  )
  .checkMetricsLevel(metrics, level, level_functions, use_default=TRUE, 
                     use_attribute=FALSE)
  # Collect all arguments into a list
  args <- list(true=true, pred=pred, location=location, k=k, alpha=alpha,
               metrics=metrics, fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
               ...)
  do.call(level_functions[[level]], args)
}




#' Compute global-level external evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data at the global level. Options include a series of fuzzy pair-counting 
#' metrics and set matching-based accuracy.
#' @inheritParams getFuzzyPartitionMetrics
#' @inheritParams getFuzzyLabel
#' @param true A vector containing the labels of the true classes. Must be a 
#'  vector of characters, integers, numerics, or a factor, but not a list.
#' @param pred A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param k The number of neighbors used when calculating the fuzzy 
#' class memberships for fuzzy metrics, or when calculating the weighted
#' accuracy.
#' @param metrics a vector of metric names to compute. 
#' @param ... Optional params for [FuzzyPartitionMetrics()] or [findSpatialKNN()].
#' @examples
#' data <- sp_toys
#' getSpatialGlobalExternalMetrics(data$label, data$p1, data[,c("x", "y")], k=6)
#' getSpatialGlobalExternalMetrics(data$label, data$p2, data[,c("x", "y")], k=6)
getSpatialGlobalExternalMetrics <- function(true, pred, location, k=6, alpha=0.5,
                                            metrics=c("SpatialRI","SpatialARI",
                                                      "SpatialWH","SpatialAWH", 
                                                      "SpatialWC","SpatialAWC",
                                                      "SpatialAccuracy"), 
                                            fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                            ...){
  argfindSpatialKNN <- .checkEllipsisArgs(fnList=list(findSpatialKNN, fuzzyPartitionMetrics), ...)[[1]] 
  argfuzzyPartitionMetrics <- .checkEllipsisArgs(fnList=list(findSpatialKNN, fuzzyPartitionMetrics), ...)[[2]] 
  if(length(intersect(metrics, c("SpatialRI","SpatialARI","SpatialWH",
                                 "SpatialAWH", "SpatialWC","SpatialAWC")))>0){
    hardTrue <- true
    hardPred <- pred
    fuzzyTrue <- do.call(getFuzzyLabel, c(argfindSpatialKNN, 
                                          list(labels=hardTrue, location=location, k=k, alpha=alpha)))
    fuzzyPred <- do.call(getFuzzyLabel, c(argfindSpatialKNN, 
                                          list(labels=hardPred, location=location, k=k, alpha=alpha)))
    res <- do.call(getFuzzyPartitionGlobalMetrics, c(argfuzzyPartitionMetrics,
                                                     list(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue,
                                                          hardPred=hardPred, fuzzyPred=fuzzyPred, 
                                                          fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                                                          metrics=c("fuzzyRI", "fuzzyARI", "fuzzyWH", 
                                                                    "fuzzyAWH", "fuzzyWC", "fuzzyAWC"))))
  }else{res <- data.frame(matrix(nrow = 1, ncol = 0))}  
  if("SpatialAccuracy" %in% metrics){
    res$SpatialAccuracy <- do.call(nnWeightedAccuracy, 
                                   c(list(true=true, pred=pred, 
                                          location=location, k=k), 
                                     argfindSpatialKNN))
  }
  if("setMatchingAccuracy" %in% metrics){
    res$setMatchingAccuracy<- .setMatchingAccuracy(true, pred)
  }
  colnames(res) <- sub("fuzzy", "Spatial",colnames(res))
  return(res[,metrics])
}

#' Compute class-level external evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data at the class/cluster level. 
#' @inheritParams getSpatialGlobalExternalMetrics
#' @param k The number of neighbors used when calculating the fuzzy 
#' class memberships for fuzzy metrics.
#' @param ... Optional params for [FuzzyPartitionMetrics()] or [findSpatialKNN()].
#' @examples
#' data <- sp_toys
#' getSpatialClassExternalMetrics(data$label, data$p1, data[,c("x", "y")], k=6)
#' getSpatialClassExternalMetrics(data$label, data$p2, data[,c("x", "y")], k=6)
getSpatialClassExternalMetrics <- function(true, pred, location, k=6, alpha=0.5,
                                           metrics=c("SpatialWH","SpatialAWH", 
                                                     "SpatialWC","SpatialAWC"), 
                                           fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                           ...){
  
  argfindSpatialKNN <- .checkEllipsisArgs(fnList=list(findSpatialKNN, fuzzyPartitionMetrics), ...)[[1]] 
  argfuzzyPartitionMetrics <- .checkEllipsisArgs(fnList=list(findSpatialKNN, fuzzyPartitionMetrics), ...)[[2]]
  hardTrue <- true
  hardPred <- pred
  fuzzyTrue <- do.call(getFuzzyLabel, c(argfindSpatialKNN, 
                                        list(labels=hardTrue, location=location, k=k, alpha=alpha)))
  fuzzyPred <- do.call(getFuzzyLabel, c(argfindSpatialKNN, 
                                        list(labels=hardPred, location=location, k=k, alpha=alpha)))
  res <- do.call(getFuzzyPartitionClassMetrics, c(argfuzzyPartitionMetrics,
                                                   list(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue,
                                                        hardPred=hardPred, fuzzyPred=fuzzyPred, 
                                                        fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                                                        metrics=c("fuzzyWH","fuzzyAWH", "fuzzyWC", "fuzzyAWC"))))
  colnames(res) <- sub("fuzzy", "Spatial",colnames(res))
  return(res[,metrics])
}

#' getSpatialElementExternalMetrics
#'
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data at the element level.
#' @inheritParams getSpatialGlobalExternalMetrics
#' @param ... Optional params for [getFuzzyPartitionElementMetrics()] or [findSpatialKNN()].
#'
#' @examples
#' data <- sp_toys
#' getSpatialElementExternalMetrics(data$label, data$p1, data[,c("x", "y")], k=6)
getSpatialElementExternalMetrics <- function(true, pred, location, k=6, alpha=0.5,
                                             metrics=c("spotAgreement"),
                                             fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                             ...){
  argfindSpatialKNN <- .checkEllipsisArgs(fnList=list(findSpatialKNN, getFuzzyPartitionElementMetrics), ...)[[1]] 
  arggetFuzzyPartitionElementMetrics <- .checkEllipsisArgs(fnList=list(findSpatialKNN, getFuzzyPartitionElementMetrics), ...)[[2]]
  hardTrue <- true
  hardPred <- pred
  fuzzyTrue <- do.call(getFuzzyLabel, c(argfindSpatialKNN, 
                                        list(labels=hardTrue, location=location, k=k, alpha=alpha)))
  fuzzyPred <- do.call(getFuzzyLabel, c(argfindSpatialKNN, 
                                        list(labels=hardPred, location=location, k=k, alpha=alpha)))
  do.call(getFuzzyPartitionElementMetrics, c(list(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue, 
                                                  hardPred=hardPred, fuzzyPred=fuzzyPred, 
                                                  fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                                                  metrics=metrics), 
                                             arggetFuzzyPartitionElementMetrics))  
}
