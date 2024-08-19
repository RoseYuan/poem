#' getPartitionClassMetrics
#' 
#' Computes a selection of external evaluation metrics for partition. The 
#' metrics are reported per class.
#'
#' @param true A vector containing the labels of the true classes. Must be a 
#'  vector of characters, integers, numerics, or a factor, but not a list.
#' @param pred A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param metrics The metrics to compute. If omitted, main metrics will be 
#'   computed.
#' @return A list of metrics.
#' @details
#' Additional details...
#' 
#' @importFrom aricode sortPairs
#' @export
#' @examples
#' true <- rep(LETTERS[1:3], each=10)
#' pred <- c(rep("A", 8), rep("B", 9), rep("C", 3), rep("D", 10))
#' getPartitionClassMetrics(true, pred)
getPartitionClassMetrics <-function(true, pred, metrics=c("WC","WH","AWC","AWH",
                                                          "FM"), ...){
  if (anyNA(true) | anyNA(pred))
    stop("NA are not supported.")
  if (is.character(true)) true <- as.factor(true)
  if (is.character(pred)) pred <- as.factor(pred)
  if (!is.atomic(true) || (!is.factor(true) && !is.integer(true)) ||
      !is.atomic(pred) || (!is.factor(pred) && !is.integer(pred)) )
    stop("true and pred must be vectors or factors but not lists.")
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

  .calWC <- function(){
    wi <- list()
    for (i in sort(unique(res$pair_c1))){
      idx <- which(res$pair_c1 == i)
      term1 <- sum(choose(res$nij[idx], 2)) 
      term3 <- choose(sum(res$nij[idx]), 2) 
      wi[i+1] <- term1 / term3
    }
    names(wi) <- res$levels$c1
    return(wi)
  }
  
  .calAWC <- function(){
    awi <- list()
    for (i in sort(unique(res$pair_c1))){
      idx <- which(res$pair_c1 == i)
      term1 <- spairs * sum(choose(res$nij[idx], 2)) 
      term2 <- choose(sum(res$nij[idx]), 2) * scol
      term3 <- choose(sum(res$nij[idx]), 2) * (spairs - scol)
      awi[i+1] <- (term1 - term2) / term3
    }
    names(awi) <- res$levels$c1
    return(awi)
  }

  .calWH <- function(){
    vj <- list()
    for (j in sort(unique(res$pair_c2))){
      idx <- which(res$pair_c2 == j)
      term1 <- sum(choose(res$nij[idx], 2)) 
      term3 <- choose(sum(res$nij[idx]), 2)
      vj[j+1] <- term1 / term3
    }
    names(vj) <- res$levels$c2
    return(vj)
  }
  
  .calAWH <- function(){
    avj <- list()
    for (j in sort(unique(res$pair_c2))){
      idx <- which(res$pair_c2 == j)
      term1 <- spairs * sum(choose(res$nij[idx], 2)) 
      term2 <- choose(sum(res$nij[idx]), 2) * srow
      term3 <- choose(sum(res$nij[idx]), 2) * (spairs - srow)
      avj[j+1] <- (term1 - term2) / term3
    }
    names(avj) <- res$levels$c2
    return(avj)
  }
  
  .calF1 <- function(){
    realClusters <- aux.conversion(true)
    predictedClusters <- aux.conversion(pred)
    if (sum(predictedClusters) == 0)
      return(0)
    a <- table(realClusters, predictedClusters)
    p <- t(apply(a, 1, function(x) x / colSums(a)))
    r <- apply(a, 2, function(x) x / rowSums(a))
    f <- 2 * r * p / (r + p)
    f[is.na(f)] <- 0
    return(as.list(apply(f, 1, max)))
  }

  res <- lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           WC = .calWC(),
           AWC = .calAWC(),
           WH = .calWH(),
           AWH = .calAWH(),
           FM = .calF1(),
           stop("Unknown metric.")
    )
  })
  return(res)
}
