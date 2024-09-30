#' Compute global-level internal evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of internal clustering evaluation metrics for spatial 
#' data at the global level. MPC, PC and PE are internal metrics for fuzzy 
#' clustering, and their implementations in package `fclust` are used.
#' 
#' @param label A vector containing the labels to be evaluated.
#' @param location A numerical matrix containing the location information, with
#' rows as samples and columns as location dimensions.
#' @param k The size of the spatial neighborhood to look at for each spot. 
#' This is used for calculating PAS and ELSA scores.
#' @param metric The metrics to compute. See below for more details.
#' 
#' @importFrom fclust MPC PC PE
#' 
#' @references Yuan, Zhiyuan, et al., 2024; 10.1038/s41592-024-02215-8
#' @references Naimi, Babak, et al., 2019; 10.1016/j.spasta.2018.10.001
#' @references Wang, et al., 2022; 10.1016/j.ins.2022.11.010
#' 
#' @return A named vector containing metric values. Possible metrics are:
#'   \item{PAS}{Proportion of abnormal spots (PAS score).}
#'   \item{ELSA}{Entropy-based Local indicator of Spatial Association (ELSA score).}
#'   \item{CHAOS}{Spatial Chaos Score.}
#'   \item{MPC}{Modified partition coefficient} 
#'   \item{PC}{Partition coefficient} 
#'   \item{PE}{Partition entropy} 
#' @export
#' @examples
getSpatialGlobalInternalMetrics <- function(label, location, k=6,
                                            metrics=c("PAS", "ELSA", "CHAOS"),
                                            ...){
  # if(length(intersect(metrics, c("MPC", "PC", "PE")))>0){require(fclust)}
  res <- lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           PAS = PAS(label, location, k=k, ...)$PAS,
           ELSA = colMeans(ELSA(label, location, k=k), na.rm = TRUE),
           CHAOS = CHAOS(label, location, BNPARAM=NULL),
           MPC = MPC(getFuzzyLabel(label, location)),
           PC = PC(getFuzzyLabel(label, location)),
           PE = PE(getFuzzyLabel(label, location)),
           stop("Unknown metric.")
    )}
    )
  res <- unlist(res)
  names(res)[names(res) == "ELSA.ELSA"] <- "ELSA"
  return(res)
}

#' Compute spot-level internal evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of internal clustering evaluation metrics for spatial 
#' data at each spot level.
#'
#' @inheritParams getSpatialGlobalInternalMetrics
#' @param metrics Possible metrics: "PAS" and "ELSA".
#' @param ... Optional params for [PAS()].
#' @return A dataframe containing the metric values for all samples in the dataset.
#' If PAS is calculated, the value is a Boolean about the abnormality of a spot.
#' If ELSA is calculated, Ea, Ec and ELSA for all spots will be returned.
#' @examples 
getSpatialInternalMetrics <- function(label, location, k=6, 
                                      metrics=c("PAS", "ELSA"), ...){
  res <- as.data.frame(lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           PAS = PAS(label, location, k=k, ...)$abnormalty,
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