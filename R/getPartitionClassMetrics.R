#' getPartitionClassMetrics
#' 
#' Short description...
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
#' 
#' @export
#' @examples
#' true <- rep(LETTERS[1:3], each=10)
#' pred <- c(rep("A", 8), rep("B", 9), rep("C", 3), rep("D", 10))
#' getPartitionClassMetrics(true, pred)
setGeneric("getPartitionClassMetrics", 
           function(true, pred, metrics, ...) standardGeneric("getPartitionClassMetrics"))

setMethod("getPartitionClassMetrics", signature="ANY",
          definition=function(true, pred, ...){
            if (anyNA(true) | anyNA(pred))
              stop("NA are not supported.")
            if (((!is.vector(true) & !is.factor(true)) | is.list(true)) | ((!is.vector(pred) & !is.factor(pred)) | is.list(pred)))
              stop("true and pred must be vectors or factors but not lists.")
            if(length(true) != length(pred)){
              stop("The two input vectors should have the same length.")
            }
            .getPartitionClassMetrics(true, pred, ...)
          })

#' @importFrom aricode sortPairs
.getPartitionClassMetrics <-function(true, pred, metrics=c("AW", "AV"), ...){
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
  
  .calAW <- function(){
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

  .calAV <- function(){
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

  res <- lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           AW = .calAW(),
           AV = .calAV(),
           stop("Unknown metric.")
    )
  })
  return(res)
}
