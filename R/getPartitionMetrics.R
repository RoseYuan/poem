#' getPartitionMetrics
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
#' getPartitionMetrics(true, pred)
setGeneric("getPartitionMetrics", 
           function(true, pred, metrics, ...) standardGeneric("getPartitionMetrics"))

setMethod("getPartitionMetrics", signature="ANY",
          definition=function(true, pred, metrics, ...){
            if (anyNA(true) | anyNA(pred))
              stop("NA are not supported.")
            if (((!is.vector(true) & !is.factor(true)) | is.list(true)) | ((!is.vector(pred) & !is.factor(pred)) | is.list(pred)))
              stop("true and pred must be vectors or factors but not lists.")
            if(length(true) != length(pred)){
              stop("The two input vectors should have the same length.")
            }
            .getPartitionMetrics(true, pred, metrics, ...)
          })

#' @importFrom aricode sortPairs
.getPartitionMetrics <-function(true, pred, metrics=c("RI","WC","WH","ARI","AWC","AWH"), ...){
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
           stop("Unknown metric.")
    )
  })
  return(res)
}
