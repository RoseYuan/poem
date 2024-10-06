#' Compute global-level external evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data at the global level. Options include a series of fuzzy pair-counting 
#' metrics and set matching-based accuracy.
#' @inheritParams getFuzzyPartitionMetrics
#' @param k The number of neighbors used when calculating the fuzzy 
#' class memberships for fuzzy metrics, or when calculating the weighted
#' accuracy.
#' @param metrics a vector of metric names to compute. 
#' @param ... Optional params for [getFuzzyPartitionMetrics()] or [nnWeightedAccuracy].
#' @examples
#' data <- sp_toys
#' getSpatialGlobalExternalMetrics(data$label, data$p1, data[,c("x", "y")], k=6)
#' getSpatialGlobalExternalMetrics(data$label, data$p2, data[,c("x", "y")], k=6)
getSpatialGlobalExternalMetrics <- function(true, pred, location, k=6, 
                                            metrics=c("SpatialRI","SpatialARI",
                                                      "SpatialWH","SpatialAWH", 
                                                      "SpatialWC","SpatialAWC",
                                                      "SpatialAccuracy"), 
                                            ...){
  if("SpatialAccuracy" %in% metrics){
    SpatialAccuracy <- nnWeightedAccuracy(true, pred, location, k=k, ...)
  }
  if("setMatchingAccuracy" %in% metrics){
    setMatchingAccuracy<- .setMatchingAccuracy(true, pred)
  }
  if(length(intersect(metrics, c("SpatialRI","SpatialARI","SpatialWH",
                                 "SpatialAWH", "SpatialWC","SpatialAWC")))>0){
    fuzzyMetrics <- getFuzzyPartitionMetrics(true, pred, location, k=k, ...)
    SpatialRI <- fuzzyMetrics$NDC
    SpatialARI <- fuzzyMetrics$ACI
    SpatialWH <- fuzzyMetrics$fuzzyWH$global
    SpatialAWH <- fuzzyMetrics$fuzzyAWH$global
    SpatialWC <- fuzzyMetrics$fuzzyWC$global
    SpatialAWC <- fuzzyMetrics$fuzzyAWC$global
  }
  res <- lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           SpatialAccuracy = SpatialAccuracy,
           SpatialRI = SpatialRI,
           SpatialARI = SpatialARI,
           SpatialWH = SpatialWH,
           SpatialAWH = SpatialAWH,
           SpatialWC = SpatialWC,
           SpatialAWC = SpatialAWC,
           setMatchingAccuracy = setMatchingAccuracy,
           stop("Unknown metric.")
    )
  })
  return(unlist(res))
}

#' Compute class-level external evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data at the class/cluster level. 
#' @inheritParams getSpatialGlobalExternalMetrics
#' @param k The number of neighbors used when calculating the fuzzy 
#' class memberships for fuzzy metrics.
#' @param ... Optional params for [getFuzzyPartitionMetrics()].
#' @examples
#' data <- sp_toys
#' getSpatialClassExternalMetrics(data$label, data$p1, data[,c("x", "y")], k=6)
#' getSpatialClassExternalMetrics(data$label, data$p2, data[,c("x", "y")], k=6)
getSpatialClassExternalMetrics <- function(true, pred, location, k=6, 
                                           metrics=c("SpatialWH","SpatialAWH", 
                                                     "SpatialWC","SpatialAWC"), 
                                           ...){
  fuzzyMetrics <- getFuzzyPartitionMetrics(true, pred, location, k=k, ...)
  res <- lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           SpatialWH =  fuzzyMetrics$fuzzyWH$perPartition,
           SpatialAWH = fuzzyMetrics$fuzzyAWH$perPartition,
           SpatialWC = fuzzyMetrics$fuzzyWC$perPartition,
           SpatialAWC = fuzzyMetrics$fuzzyAWC$perPartition,
           stop("Unknown metric.")
    )
  })
  return(res)
}