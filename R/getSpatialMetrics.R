#' getSpatialGlobalInternalMetrics
#' 
#' Computes a selection of internal clustering evaluation metrics for spatial 
#' data at the global level. MPC, PC and PE are internal metrics for fuzzy 
#' clustering. Use the implementations in package `fclust`.
#' @references Wang, et al. "A survey of fuzzy clustering validity evaluation 
#' methods." Information Sciences 618 (2022): 270-297.
#' @return A vector with the following values (see refs for details):
#'   \item{PAS}{}
#'   \item{ELSA}{}
#'   \item{CHAOS}{}
#'   \item{MPC}{Modified partition coefficient} 
#'   \item{PC}{Partition coefficient} 
#'   \item{PE}{Partition entropy} 
getSpatialGlobalInternalMetrics <- function(label, location, k=6,
                                            metrics=c("PAS", "ELSA", "CHAOS"),
                                            ...){
  if(length(intersect(metrics, c("MPC", "PC", "PE")))>0){require(fclust)}
  res <- lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           PAS = PAS(label, location, k=k, ...)$PAS,
           ELSA = colMeans(ELSA(label, location, k=k), na.rm = TRUE),
           CHAOS = CHAOS(label, location, BNPARAM=NULL),
           MPC = fclust::MPC(getFuzzyLabel(label, location)),
           PC = fclust::PC(getFuzzyLabel(label, location)),
           PE = fclust::PE(getFuzzyLabel(label, location)),
           stop("Unknown metric.")
    )}
    )
  res <- unlist(res)
  names(res)[names(res) == "ELSA.ELSA"] <- "ELSA"
  return(res)
}

#' getSpatialInternalMetrics
#' 
#' Computes a selection of internal clustering evaluation metrics for spatial 
#' data at the spot level.
getSpatialInternalMetrics <- function(label, location, k=6, 
                                      metrics=c("PAS", "ELSA"), ...){
  res <- as.data.frame(lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           PAS.abnormal = PAS(label, location, k=k, ...)$abnormalty,
           ELSA = ELSA(label, location, k=k),
           stop("Unknown metric.")
           )})
    )
  colnames(res)[colnames(res) == "ELSA.ELSA"] <- "ELSA"
  return(res)
}

#' getSpatiaGloballExternalMetrics
#' 
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data at the global level.
getSpatialGlobalExternalMetrics <- function(true, pred, location, k=6, 
                                            metrics=c("SpatialRI","SpatialARI",
                                                      "SpatialWH","SpatialAWH", 
                                                      "SpatialWC","SpatialAWC",
                                                      "SpatialAccuracy"), 
                                            ...){
  if("SpatialAccuracy" %in% metrics){
    SpatialAccuracy <- nnWeightedAccuracy(true, pred, location, k=k, ...)
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
           stop("Unknown metric.")
    )
  })
  return(unlist(res))
}

#' getSpatialClassExternalMetrics
#' 
#' Computes a selection of external clustering evaluation metrics for spatial 
#' data at the class/cluster level.
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