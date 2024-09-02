#' getPartitionMetrics
#' 
#' Computes a selection of external evaluation metrics for partition. The 
#' metrics are reported per dataset.
#'
#' @param true A vector containing the labels of the true classes. Must be a 
#'  vector of characters, integers, numerics, or a factor, but not a list.
#' @param pred A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param metrics The metrics to compute. If omitted, main metrics will be 
#'   computed. RI: Rand Index; WC: Wallace Completeness; WH: Wallace 
#'   Homogeneity; ARI:Adjusted Rand Index; AWC: Adjusted Wallace Completeness; 
#'   AWH: Adjusted Wallace Homogeneity; MI: Mutual Information; AMI: Adjusted 
#'   Mutual Information; VI: Variation of Information; EH: (Entropy-based)
#'   Homogeneity; EC: (Entropy-based) Completeness; VM: V-measure; 
#'   FM: F-measure/weighted average F1 score; VDM: Van Dongen Measure; 
#'   MHM: Meila-Heckerman Measure; MMM: Maximum-Match Measure; 
#'   Mirkin: Mirkin Metric.
#' @return A list of metric results.
#' @details
#' Additional details...
#' 
#' @importFrom aricode sortPairs AMI
#' @importFrom clevr mutual_info variation_info homogeneity completeness v_measure
#' @importFrom mclustcomp mclustcomp
#' @importFrom FlowSOM FMeasure
#' @export
#' @examples
#' true <- rep(LETTERS[1:3], each=10)
#' pred <- c(rep("A", 8), rep("B", 9), rep("C", 3), rep("D", 10))
#' getPartitionMetrics(true, pred)
getPartitionMetrics <-function(true, pred, 
                                metrics=c("RI","WC","WH","ARI","AWC","AWH",
                                          "MI","AMI","VI","EH","EC","VM",
                                          "FM","VDM","Mirkin","MHM","MMM"),
                                ...){
  if (anyNA(true) | anyNA(pred)) stop("NA are not supported.")
  if (is.character(true)) true <- as.factor(true)
  if (is.character(pred)) pred <- as.factor(pred)
  if (!is.atomic(true) || (!is.factor(true) && !is.integer(true)) ||
      !is.atomic(pred) || (!is.factor(pred) && !is.integer(pred)) )
    stop("true and pred must be vectors of factors but not lists.")
  if(length(true) != length(pred)){
    stop("The two input vectors should have the same length.")
  }
  res <- sortPairs(true, pred)
  n <- length(true)

  spairs <- n*(n-1)/2 # N
  stot <- sum(choose(res$nij, 2), na.rm=TRUE) # T
  srow <- sum(choose(res$ni., 2), na.rm=TRUE) # P
  scol <- sum(choose(res$n.j, 2), na.rm=TRUE) # Q
  a <- stot
  b <- srow-stot
  c <- scol-stot
  d <- spairs+stot-srow-scol
  
  res <- lapply(setNames(metrics,metrics), FUN=function(m){
    switch(m,
           WC = a/(a+b),
           AWC = (a*d-b*c)/((a+b)*(b+d)),
           WH = a/(a+c),
           AWH = (a*d-b*c)/((a+c)*(c+d)),
           RI = (a+d)/(a+b+c+d),
           ARI = 2*(a*d-b*c)/((a+b)*(b+d)+(a+c)*(c+d)),
           MI = mutual_info(true, pred, ...),
           AMI = AMI(true, pred, ...),
           VI = variation_info(true, pred, ...),
           EH = homogeneity(true, pred, ...),
           EC = completeness(true, pred, ...),
           VM = v_measure(true, pred, ...),
           VDM = mclustcomp(as.vector(true), as.vector(pred), types = "vdm")$scores,
           Mirkin = mclustcomp(as.vector(true), as.vector(pred), types = "mirkin")$scores,
           MHM = mclustcomp(as.vector(true), as.vector(pred), types = "mhm")$scores,
           MMM = mclustcomp(as.vector(true), as.vector(pred), types = "mmm")$scores,
           FM = FMeasure(.aux.conversion(true), 
                         .aux.conversion(pred), silent = TRUE),
           stop("Unknown metric.")
    )
  })
  return(res)
}