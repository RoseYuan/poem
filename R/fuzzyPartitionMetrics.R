#' Compute fuzzy-fuzzy versions of pair-sorting partition metrics
#' 
#' Computes fuzzy versions of pair-sorting partition metrics. This is largely 
#' based on the permutation-based implementation by Antonio D'Ambrosio from the 
#' `ConsRankClass` package, modified to also compute the fuzzy versions of the 
#' adjusted Wallace indices, implement multithreading, and adjust the number of
#' permutations according to their variability.
#'
#' @param P A object coercible to a numeric matrix with membership probability 
#'   of elements (rows) in ground-truth classes (columns).
#' @param Q A object coercible to a numeric matrix with membership probability 
#'   of elements (rows) in predicted clusters (columns). Must have the same 
#' number of rows as `P`.
#' @param nperms The number of permutations (for correction for chance). If 
#'   NULL (default), a first set of 10 permutations will be run to estimate 
#'   whether the variation across permutations is above 0.0025, in which case 
#'   more (max 1000) permutations will be run.
#' @param computeWallace Logical; whether to compute the individual fuzzy 
#'   versions of the Wallace indices (increases running time).
#' @param verbose Logical; whether to print info and warnings, including the 
#'   standard error of the mean across permutations (giving an idea of the 
#'   precision of the adjusted metrics).
#' @param returnElementPairAccuracy Logical. If TRUE, returns the per-element
#'   pair accuracy instead of the various parition-level and dataset-level 
#'   metrics. Default FALSE.
#' @param BPPARAM BiocParallel params for multithreading (default none)
#' @param tnorm Which type of t-norm operation to use for class membership of
#'   pairs (either product, min, or lukasiewicz) when calculating the Wallace 
#'   indices. Does not influence the NDC/ACI metrics.
#' 
#' @references Hullermeier et al. 2012; 10.1109/TFUZZ.2011.2179303;
#' @references D'Ambrosio et al. 2021; 10.1007/s00357-020-09367-0
#' 
#' @author Pierre-Luc Germain
#' 
#' @return When `returnElementPairAccuracy` is `FALSE`, return a list of 
#' metrics:
#'   \item{NDC}{Hullermeier's NDC (fuzzy rand index)}
#'   \item{ACI}{Ambrosio's Adjusted Concordance Index (ACI), i.e. a 
#'     permutation-based fuzzy version of the adjusted Rand index.}
#'   \item{fuzzyWH}{Fuzzy Wallace Homogeneity index}
#'   \item{fuzzyWC}{Fuzzy Wallace Completeness index}
#'   \item{fuzzyAWH}{Adjusted fuzzy Wallace Homogeneity index}
#'   \item{fuzzyAWC}{Adjusted fuzzy Wallace Completeness index}
#'   
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom stats sd
#' @export
#' @examples
#' # generate fuzzy partitions:
# m1 <- matrix(c(0.95, 0.025, 0.025,
#                0.98, 0.01, 0.01,
#                0.96, 0.02, 0.02,
#                0.95, 0.04, 0.01,
#                0.95, 0.01, 0.04,
#                0.99, 0.005, 0.005,
#                0.025, 0.95, 0.025,
#                0.97, 0.02, 0.01,
#                0.025, 0.025, 0.95),
#                ncol = 3, byrow=TRUE)
# m2 <- matrix(c(0.95, 0.025, 0.025,
#                0.98, 0.01, 0.01,
#                0.96, 0.02, 0.02,
#                0.025, 0.95, 0.025,
#                0.02, 0.96, 0.02,
#                0.01, 0.98, 0.01,
#                0.05, 0.05, 0.95,
#                0.02, 0.02, 0.96,
#                0.01, 0.01, 0.98),
#                ncol = 3, byrow=TRUE)
# colnames(m1) <- colnames(m2) <- LETTERS[seq_len(3)]
# fuzzyPartitionMetrics(m1,m2)
fuzzyPartitionMetrics <- function(P, Q, computeWallace=TRUE, nperms=NULL,
                                  verbose=TRUE, returnElementPairAccuracy=FALSE,
                                  BPPARAM=BiocParallel::SerialParam(), 
                                  tnorm=c("product","min","lukasiewicz")){ 
  
  tnorm <- match.arg(tnorm)
  if(is.data.frame(P)) P <- as.matrix(P)
  if(is.data.frame(Q)) Q <- as.matrix(Q)
  stopifnot(is.matrix(P) && (is.numeric(P) | is.integer(P)))
  stopifnot(is.matrix(Q) && (is.numeric(Q) | is.integer(Q)))
  stopifnot(nrow(P)==nrow(Q))
  
  m <- nrow(P)
  ncomp <- m*((m-1)/2)
  if(verbose && m>=2000){
    if(returnElementPairAccuracy){
      os <- 8*(3*m^2)
    }else{
      nSim <- min(BiocParallel::bpnworkers(BPPARAM),
                  ifelse(is.null(nperms),10,nperms))
      os <- 8*((3+nSim)*m^2)
      if(computeWallace) os <- os + 8*ncomp*(1+max(ncol(Q),ncol(P)))
    }
    class(os) <- "object_size"
    message("Projected memory usage: ", format(os, units = "auto"))
  }
  
  ep <- as.matrix(1-(0.5*dist(P,method="manhattan")))
  eq <- as.matrix(1-(0.5*dist(Q,method="manhattan")))
  
  if(returnElementPairAccuracy){
    return(1-rowSums(abs(ep - eq))/(ncol(ep)-1))
  }
  
  # Hullermeier's NDC
  diff <- abs(ep[lower.tri(ep)] - eq[lower.tri(eq)])
  NDC <- 1 - ( sum( diff )/(ncomp) )
  
  membershipFn <- .membershipFn(tnorm)
  # precompute the pairs' class membership (for increased speed in permutations)
  Ppairs <- apply(P, 2, FUN=function(p) membershipFn(p)[lower.tri(ep)])
  
  getFWallace <- function(emA, emB, fuzzyClB, Bpairs=NULL, diff=NULL){
    if(is.null(diff)) diff <- abs(emA[lower.tri(emA)]-emB[lower.tri(emB)])
    a <- vapply(setNames(seq_len(ncol(fuzzyClB)),colnames(fuzzyClB)),
                FUN.VALUE=numeric(2L), FUN=function(i){
      # get the degree to which the members of each pair are of the given 
      # class in B
      if(is.null(Bpairs)){
        Bpair <- as.matrix(membershipFn(fuzzyClB[,i])[lower.tri(emA)])
      }else{
        Bpair <- Bpairs[,i]
      }
  # compute 1 - the distance between pair concordances weighted by their being 
  # of the same class in B
      c(c=sum((1-diff)*Bpair), n=sum(Bpair))
    })
    a[2,which(a[2,]==0)] <- 1  # avoid NaNs for singletons
    list( global=sum(a[1,])/sum(a[2,], na.rm=TRUE),
          perPartition=a[1,]/a[2,] )
  }
  
  if(computeWallace){
    W1 <- getFWallace(ep,eq,Q,diff=diff)
    W2 <- getFWallace(eq,ep,P,Bpairs=Ppairs,diff)
  }
  
  # get metrics for the permutations
  
  # fn for one permutation:
  onePerm <- function(col){
    p <- allp[,col]
    permutedEQ <- eq[p,p]
    diff <- abs(ep[lower.tri(ep)] - permutedEQ[lower.tri(eq)])
    NDC <- 1 -  sum( diff )/ncomp
    if(!computeWallace) return(NDC)
    W1 <- getFWallace(ep, permutedEQ, Q[p,], diff=diff)
    W2 <- getFWallace(permutedEQ, ep, P, Bpairs=Ppairs, diff=diff)
    list(NDC=NDC, W1=W1, W2=W2)
  }
  
  res1 <- NULL
  if(is.null(nperms)){
    # try few permutations to estimate if more are needed
    allp <- apply(matrix( runif(m*10), nrow=m ), 2, order)
    res1 <- bplapply(seq_len(10), BPPARAM=BPPARAM, onePerm)
    NDCs <- vapply(res1, FUN.VALUE=numeric(1L), FUN=\(x) x[[1]])
    SE <- sd(NDCs)/sqrt(length(res1))
    if(SE>0.0025) nperms <- 100
    if(SE>0.01) nperms <- 1000
    if(verbose && !is.null(nperms))
      message("Running ", nperms, " extra permutations.")
  }
  
  if(!is.null(nperms)){
    # generate more permutations
    allp <- apply(matrix( runif(m*nperms), nrow=m ), 2, order)
  
    res <- bplapply(seq_len(nperms), BPPARAM=BPPARAM, onePerm)
    if(!is.null(res1)) res <- c(res1,res)
  }else{
    res <- res1
  }
  
  NDCs <- vapply(res, FUN.VALUE=numeric(1L), FUN=\(x) x[[1]])
  SE <- sd(NDCs)/sqrt(nperms)
  
  if(verbose){
    message("Standard error of the mean NDC across permutations:", 
            format(SE, digits=3))
    if(!isFALSE(SE>0.0025))
      message("You might want to increase the number of permutations to ",
              "increase the robustness of the adjusted metrics.")
  }
  
  adj <- function(x, m) (x-m)/(1-m)
  ACI <- adj(NDC,mean(NDCs))
  if(!computeWallace) return(list(NDC=NDC, ACI=ACI))
  
  W1m <- mean(vapply(res, \(x) x[[2]][[1]], numeric(1L)))
  W2m <- mean(vapply(res, \(x) x[[3]][[1]], numeric(1L)))
  W1pm <- rowMeans(vapply(res, \(x) x[[2]][[2]], numeric(ncol(Q))))
  W2pm <- rowMeans(vapply(res, \(x) x[[3]][[2]], numeric(ncol(P))))
  names(W1pm) <- colnames(Q)
  names(W2pm) <- colnames(P)
  AW1 <- list( global=adj(W1$global,W1m),
               perPartition=mapply(x=W1$perPartition, m=W1pm, FUN=adj) )
  AW2 <- list( global=adj(W2$global,W2m),
               perPartition=mapply(x=W2$perPartition, m=W2pm, FUN=adj) )

  return(list(NDC=NDC, ACI=ACI, fuzzyWH=W1, fuzzyWC=W2,
              fuzzyAWH=AW1, fuzzyAWC=AW2))
}


.membershipFn <- function(tnorm=c("product","min","lukasiewicz")){
  tnorm <- match.arg(tnorm)
  switch(tnorm,
     product=tcrossprod,
     lukasiewicz=function(p){
       vapply(seq_along(p), FUN=function(i) pmax(0,p+p[i]-1), 
              numeric(length(p)))
     },
     min=function(p){
       vapply(seq_along(p), FUN=function(i) pmin(p,p[i]), numeric(length(p)))
     },
     stop("Undefined t-norm")
  )
}


#' Compute fuzzy-hard versions of pair-sorting partition metrics
#' 
#' Computes fuzzy-hard versions of pair-sorting partition metrics to compare a 
#' hard clustering with both a fuzzy and hard truth. This was especially 
#' designed for cases where the fuzzy truth represents an uncertainty of a hard
#' truth. Briefly put, the maximum of the pair concordance between the 
#' clustering and either the hard or the fuzzy truth is used, and the hard truth
#'  is used to compute completeness. See \code{\link{fuzzyPartitionMetrics}} for
#'  the more standard implementation of the metrics.
#'
#' @param hardPred An atomic vector coercible to a factor or integer vector 
#'  containing the predicted hard labels.
#' @param hardTrue An atomic vector coercible to a factor or integer vector 
#'  containing the true hard labels. Must have the same length as `hardPred`.
#' @param fuzzyTrue A object coercible to a numeric matrix with membership 
#'   probability of elements (rows) in clusters (columns). Must have the same 
#'   number of rows as the length of `hardTrue`. Also note that the columns of
#'   `fuzzyTrue` should be in the order of the levels (or integer values) of
#'   `hardTrue`.
#' @param nperms The number of permutations (for correction for chance). If 
#'   NULL (default), a first set of 10 permutations will be run to estimate 
#'   whether the variation across permutations is above 0.0025, in which case 
#'   more (max 1000) permutations will be run.
#' @param returnElementPairAccuracy Logical. If TRUE, returns the per-element
#'   pair accuracy instead of the various parition-level and dataset-level
#'    metrics. Default FALSE.
#' @param lowMemory Logical; whether to use the slower, low-memory algorithm.
#'   By default this is enabled if the projected memory usage is higher than 
#'   ~2GB.
#' @param verbose Logical; whether to print info and warnings, including the 
#'   standard error of the mean across permutations (giving an idea of the 
#'   precision of the adjusted metrics).
#' @param BPPARAM BiocParallel params for multithreading (default none)
#' 
#' @references Hullermeier et al. 2012; 10.1109/TFUZZ.2011.2179303;
#' @references D'Ambrosio et al. 2021; 10.1007/s00357-020-09367-0
#' 
#' @seealso [fuzzyPartitionMetrics()].
#' 
#' @author Pierre-Luc Germain
#'
#' @return A list of metrics:
#'   \item{NDC}{Hullermeier's NDC (fuzzy rand index)}
#'   \item{ACI}{Ambrosio's Adjusted Concordance Index (ACI), i.e. a 
#'     permutation-based fuzzy version of the adjusted Rand index.}
#'   \item{fuzzyWH}{Fuzzy Wallace Homogeneity index}
#'   \item{fuzzyWC}{Fuzzy Wallace Completeness index}
#'   \item{fuzzyAWH}{Adjusted fuzzy Wallace Homogeneity index}
#'   \item{fuzzyAWC}{Adjusted fuzzy Wallace Completeness index}
#' @importFrom BiocParallel SerialParam bplapply bpnworkers
#' @importFrom stats dist setNames runif
#' @importFrom Matrix sparseMatrix
#' @importFrom stats sd
#' @export
#' @examples
#' # generate a fuzzy truth:
#' fuzzyTrue <- matrix(c(
#'   0.95, 0.025, 0.025, 
#'   0.98, 0.01, 0.01, 
#'   0.96, 0.02, 0.02, 
#'   0.95, 0.04, 0.01, 
#'   0.95, 0.01, 0.04, 
#'   0.99, 0.005, 0.005, 
#'   0.025, 0.95, 0.025, 
#'   0.97, 0.02, 0.01, 
#'   0.025, 0.025, 0.95), 
#'   ncol = 3, byrow=TRUE)
#' # a hard truth:
#' hardTrue <- apply(fuzzyTrue,1,FUN=which.max)
#' # some predicted labels:
#' hardPred <- c(1,1,1,1,1,1,2,2,2)
#' fuzzyHardMetrics(hardTrue, fuzzyTrue, hardPred, nperms=3)
fuzzyHardMetrics <- function(hardTrue, fuzzyTrue, hardPred, nperms=NULL, 
                             returnElementPairAccuracy=FALSE,
                             lowMemory=NULL, verbose=TRUE,
                             BPPARAM=BiocParallel::SerialParam()){ 
  
  stopifnot(is.atomic(hardPred))
  hardPredVector <- hardPred <- as.integer(as.factor(hardPred))
  if(is.atomic(hardTrue)){
    hardTrueVector <- hardTrue <- as.integer(as.factor(hardTrue))
    hardTrue <- Matrix::sparseMatrix(seq_along(hardTrue), hardTrue, x=1L)
  }else{
    if(is.data.frame(hardTrue)) hardTrue <- as.matrix(hardTrue)
    stopifnot(is.matrix(hardTrue) && 
                (is.numeric(hardTrue) | is.integer(hardTrue)))
    hardTrueVector <- apply(hardTrue, 1, FUN=which.max)
  }
  hardPred <- Matrix::sparseMatrix(seq_along(hardPred), hardPred, x=1L)
  
  if(is.data.frame(fuzzyTrue)) fuzzyTrue <- as.matrix(fuzzyTrue)
  
  stopifnot(is.matrix(fuzzyTrue) && 
              (is.numeric(fuzzyTrue) | is.integer(fuzzyTrue)))
  stopifnot(nrow(hardTrue)==nrow(fuzzyTrue))
  stopifnot(nrow(hardPred)==nrow(hardTrue))
  stopifnot(ncol(hardTrue)==ncol(fuzzyTrue))
  
  m <- nrow(fuzzyTrue)
  ncomp <- m*((m-1)/2)
  if(isTRUE(lowMemory) || (verbose && m>=2000) ){
    nSim <- min(BiocParallel::bpnworkers(BPPARAM),
                ifelse(is.null(nperms),10,nperms))
    os <- 8*((3+nSim)*m^2)
    class(os) <- "object_size"
    if(isTRUE(lowMemory) || os > 2*10^9){
      if(verbose) message("Using the low-memory algorithm...")
      return(
        fuzzyHardMetrics2(hardTrue=hardTrueVector, fuzzyTrue=fuzzyTrue,
                          hardPred=hardPredVector, nperms=nperms, 
                          returnElementPairAccuracy=returnElementPairAccuracy,
                          verbose=verbose, BPPARAM=BPPARAM))
    }
    message("Projected memory usage: ", format(os, units = "auto"))
  }
  
  # compute pair concordance for all the three labelings
  eq <- as.matrix(1-(0.5*dist(hardPred, method="manhattan")))
  ep <- 1-(0.5*dist(fuzzyTrue, method="manhattan"))
  ep2 <- 1-(0.5*dist(hardTrue, method="manhattan"))
  
  # compute the minimum difference of the pred with either truth
  diff <- pmin( abs(as.matrix(ep) - eq),
                abs(as.matrix(ep2) - eq) )
  
  if(returnElementPairAccuracy){
    return(1-rowSums(diff)/(ncol(diff)-1))
  }
  
  NDC <- 1 - ( sum( diff )/(2*ncomp) )
  
  getFWallace <- function(membership, diff){
    a <- vapply(split(seq_along(membership),membership),
                FUN.VALUE=numeric(2L), FUN=function(i){
      c(c=sum(diff[i,i]), n=(length(i)^2-length(i)))
    })
    a[2,which(a[2,]==0)] <- 1  # avoid NaNs for singletons
    list( global=1-sum(a[1,])/sum(a[2,], na.rm=TRUE),
          perPartition=1-a[1,]/a[2,] )
  }
  
  W1 <- getFWallace(hardPredVector, diff)
  W2 <- getFWallace(hardTrueVector, diff)
  # get metrics for the permutations
  
  # fn for one permutation:
  onePerm <- function(col){
    p <- allp[,col]
    permutedEQ <- eq[p,p]
    diff <- pmin( abs(as.matrix(ep) - permutedEQ),
                  abs(as.matrix(ep2) - permutedEQ) )
    NDC <- 1 -  sum( diff )/(2*ncomp)
    W1 <- getFWallace(hardPredVector[p], diff)
    W2 <- getFWallace(hardTrueVector, diff)
    if(m>2000){
      rm(diff, permutedEQ)
      gc(verbose=FALSE)
    }
    list(NDC=NDC, W1=W1, W2=W2)
  }
  
  res1 <- NULL
  if(is.null(nperms)){
    # try few permutations to estimate if more are needed
    allp <- apply(matrix( runif(m*10), nrow=m ), 2, order)
    res1 <- bplapply(seq_len(10), BPPARAM=BPPARAM, onePerm)
    NDCs <- vapply(res1, FUN.VALUE=numeric(1L), FUN=\(x) x[[1]])
    SE <- sd(NDCs)/sqrt(length(res1))
    if(SE>0.0025) nperms <- 100
    if(SE>0.01) nperms <- 1000
    if(verbose && !is.null(nperms))
      message("Running ", nperms, " extra permutations.")
  }
  
  if(!is.null(nperms)){
    # generate more permutations
    allp <- apply(matrix( runif(m*nperms), nrow=m ), 2, order)
    
    res <- bplapply(seq_len(nperms), BPPARAM=BPPARAM, onePerm)
    if(!is.null(res1)) res <- c(res1,res)
    NDCs <- vapply(res, FUN.VALUE=numeric(1L), FUN=\(x) x[[1]])
    SE <- sd(NDCs)/sqrt(nperms)
  }else{
    res <- res1
  }
  
  if(verbose){
    message("Standard error of the mean NDC across permutations:", 
            format(SE, digits=3))
    if(!isFALSE(SE>0.0025))
      message("You might want to increase the number of permutations to ",
              "increase the robustness of the adjusted metrics.")
  }
  
  adj <- function(x, m) (x-m)/(1-m)
  
  ACI <- adj(NDC,mean(NDCs))
  
  W1m <- mean(vapply(res, \(x) x[[2]][[1]], numeric(1L)))
  W2m <- mean(vapply(res, \(x) x[[3]][[1]], numeric(1L)))
  W1pm <- rowMeans(vapply(res, \(x) x[[2]][[2]],
                          FUN.VALUE=numeric(length(unique(hardPredVector)))))
  W2pm <- rowMeans(vapply(res, \(x) x[[3]][[2]], numeric(ncol(fuzzyTrue))))

  AW1 <- list( global=adj(W1$global,W1m),
               perPartition=mapply(x=W1$perPartition, m=W1pm, FUN=adj) )
  AW2 <- list( global=adj(W2$global,W2m),
               perPartition=mapply(x=W2$perPartition, m=W2pm, FUN=adj) )
  
  return(list(NDC=NDC, ACI=ACI, fuzzyWH=W1, fuzzyWC=W2,
              fuzzyAWH=AW1, fuzzyAWC=AW2))
}


#' Compute fuzzy-hard metrics with lower memory requirement
#' 
#' This is a slightly slower, but low-memory version of 
#' \code{\link{fuzzyHardMetrics}}.
#'
#' @param hardPred An atomic vector coercible to a factor or integer vector 
#'  containing the predicted hard labels.
#' @param hardTrue An atomic vector coercible to a factor or integer vector 
#'  containing the true hard labels. Must have the same length as `hardPred`.
#' @param fuzzyTrue A object coercible to a numeric matrix with membership 
#'   probability of elements (rows) in clusters (columns). Must have the same 
#'   number of rows as the length of `hardTrue`. Also note that the columns of
#'   `fuzzyTrue` should be in the order of the levels (or integer values) of
#'   `hardTrue`.
#' @param nperms The number of permutations (for correction for chance). If 
#'   NULL (default), a first set of 10 permutations will be run to estimate 
#'   whether the variation across permutations is above 0.0025, in which case 
#'   more (max 1000) permutations will be run.
#' @param returnElementPairAccuracy Logical. If TRUE, returns the per-element
#'   pair accuracy instead of the various parition-level and dataset-level 
#'   metrics. Default FALSE.
#' @param verbose Logical; whether to print info and warnings, including the 
#'   standard error of the mean across permutations (giving an idea of the 
#'   precision of the adjusted metrics).
#' @param BPPARAM BiocParallel params for multithreading (default none)
#' 
#' @references Hullermeier et al. 2012; 10.1109/TFUZZ.2011.2179303;
#' @references D'Ambrosio et al. 2021; 10.1007/s00357-020-09367-0
#' 
#' @seealso [fuzzyHardMetrics()]
#' 
#' @author Pierre-Luc Germain
#' @keywords internal
#' @return A list of metrics:
#'   \item{NDC}{Hullermeier's NDC (fuzzy rand index)}
#'   \item{ACI}{Ambrosio's Adjusted Concordance Index (ACI), i.e. a 
#'     permutation-based fuzzy version of the adjusted Rand index.}
#'   \item{fuzzyWH}{Fuzzy Wallace Homogeneity index}
#'   \item{fuzzyWC}{Fuzzy Wallace Completeness index}
#'   \item{fuzzyAWH}{Adjusted fuzzy Wallace Homogeneity index}
#'   \item{fuzzyAWC}{Adjusted fuzzy Wallace Completeness index}
#'   
#' @importFrom BiocParallel SerialParam bplapply bpnworkers
#' @importFrom stats dist setNames runif
#' @importFrom Matrix sparseMatrix
#' @importFrom stats sd
#' @examples
#' # generate a fuzzy truth:
#' fuzzyTrue <- matrix(c(
#' 0.95, 0.025, 0.025,
#' 0.98, 0.01, 0.01,
#' 0.96, 0.02, 0.02,
#' 0.95, 0.04, 0.01,
#' 0.95, 0.01, 0.04,
#' 0.99, 0.005, 0.005,
#' 0.025, 0.95, 0.025,
#' 0.97, 0.02, 0.01,
#' 0.025, 0.025, 0.95),
#' ncol = 3, byrow=TRUE)
#' # a hard truth:
#' hardTrue <- apply(fuzzyTrue,1,FUN=which.max)
#' # some predicted labels:
#' hardPred <- c(1,1,1,1,1,1,2,2,2)
#' poem:::fuzzyHardMetrics2(hardTrue, fuzzyTrue, hardPred, nperms=3)
fuzzyHardMetrics2 <- function(hardTrue, fuzzyTrue, hardPred, nperms=10, 
                             returnElementPairAccuracy=FALSE, verbose=TRUE,
                             BPPARAM=BiocParallel::SerialParam()){ 
  
  stopifnot(is.atomic(hardPred))
  stopifnot(is.atomic(hardTrue))
  hardPred <- as.factor(hardPred)
  hardTrue <- as.factor(hardTrue)
  if(is.data.frame(fuzzyTrue)) fuzzyTrue <- as.matrix(fuzzyTrue)
  
  stopifnot(is.matrix(fuzzyTrue) && 
              (is.numeric(fuzzyTrue) | is.integer(fuzzyTrue)))
  stopifnot(length(hardTrue)==nrow(fuzzyTrue))
  stopifnot(length(hardPred)==nrow(hardTrue))
  stopifnot(length(levels(hardTrue))==ncol(fuzzyTrue))
  
  m <- nrow(fuzzyTrue)
  fuzzyTrue <- t(fuzzyTrue)
  allp <- apply(matrix( runif(m*10), nrow=m ), 2, order)
  
  if(returnElementPairAccuracy){
    return(vapply(seq_len(m), FUN.VALUE=numeric(1L), FUN=function(i){
      # compute pair agreement for all the three labelings
      eq <- as.integer(hardPred[-i]==hardPred[i])
      ep <- 1-colSums(abs(fuzzyTrue[,-i]-fuzzyTrue[,i]))/2
      ep2 <- as.integer(hardTrue[-i]==hardTrue[i])
      diff <- pmin( abs(eq-ep), abs(eq-ep2) )
      1-mean(diff)
    }))
  }
  
  doCalcs <- function(hardPred){
    a <- lapply(seq_len(m), FUN=function(i){
      # compute pair agreement for all the three labelings
      eq <- as.integer(hardPred[-i]==hardPred[i])
      ep <- 1-colSums(abs(fuzzyTrue[,-i]-fuzzyTrue[,i]))/2
      ep2 <- as.integer(hardTrue[-i]==hardTrue[i])
      diff <- pmin( abs(eq-ep), abs(eq-ep2) )
      # completeness
      ac <- vapply(split(seq_along(hardTrue[-i]), hardTrue[-i]), 
                   FUN.VALUE=numeric(2L), FUN=function(x){
        c(diff=sum(diff[x]), n=length(x))
      })
      ac[2,which(ac[2,]==0)] <- 1  # avoid NaNs for singletons
      
      # homogeneity
      a <- vapply(split(seq_along(hardPred[-i]), hardPred[-i]),
                  FUN.VALUE=numeric(2L), FUN=function(x){
        c(diff=sum(diff[x]), n=length(x))
      })

      list( NDC=1-mean(diff), aH=a, aC=ac )
    })
    
    Cdiff <- rowSums(vapply(a, FUN=\(x) x$aC[1,],
                            FUN.VALUE=numeric(nrow(fuzzyTrue))))
    Cn <- rowSums(vapply(a, FUN=function(x) x$aC[2,],
                         FUN.VALUE=numeric(nrow(fuzzyTrue))))
    Hdiff <- rowSums(vapply(a, FUN=function(x) x$aH[1,],
                            FUN.VALUE=numeric(length(levels(hardPred)))))
    Hn <- rowSums(vapply(a, FUN=function(x) x$aH[2,],
                         FUN.VALUE=numeric(length(levels(hardPred)))))
    # avoid NaNs for singletons
    Hn[which(Hn==0)] <- 1L
    Cn[which(Cn==0)] <- 1L
    
    NDC <- mean(unlist(lapply(a, FUN=function(x) x$NDC)))
    HO <- list( global=1-sum(Hdiff)/sum(Hn),
                perPartition=1-Hdiff/Hn)
    CO <- list( global=1-sum(Cdiff)/sum(Cn),
                perPartition=1-Cdiff/Cn)
    list(NDC=NDC, HO=HO, CO=CO)
  }
  
  val <- doCalcs(hardPred)

  onePerm <- function(col){
    doCalcs(hardPred[allp[,col]])
  }
    
  res1 <- NULL
  if(is.null(nperms)){
    # try few permutations to estimate if more are needed
    allp <- apply(matrix( runif(m*10), nrow=m ), 2, order)
    res1 <- bplapply(seq_len(10), BPPARAM=BPPARAM, onePerm)
    NDCs <- vapply(res1, FUN.VALUE=numeric(1L), FUN=\(x) x[[1]])
    SE <- sd(NDCs)/sqrt(length(res1))
    if(SE>0.0025) nperms <- 100
    if(SE>0.01) nperms <- 1000
    if(verbose && !is.null(nperms))
      message("Running ", nperms, " extra permutations.")
  }
  
  if(!is.null(nperms)){
    # generate more permutations
    allp <- apply(matrix( runif(m*nperms), nrow=m ), 2, order)
    res <- bplapply(seq_len(nperms), BPPARAM=BPPARAM, onePerm)
    if(!is.null(res1)) res <- c(res1,res)
    NDCs <- vapply(res, FUN.VALUE=numeric(1L), FUN=\(x) x[[1]])
    SE <- sd(NDCs)/sqrt(nperms)
  }else{
    res <- res1
  }
  
  if(verbose){
    message("Standard error of the mean NDC across permutations:", 
            format(SE, digits=3))
    if(!isFALSE(SE>0.0025))
      message("You might want to increase the number of permutations to ",
              "increase the robustness of the adjusted metrics.")
  }
  
  adj <- function(x, m) (x-m)/(1-m)
  
  ACI <- adj(val[[1]],mean(NDCs))
  
  W1m <- mean(vapply(res, \(x) x[[2]][[1]], numeric(1L)))
  W2m <- mean(vapply(res, \(x) x[[3]][[1]], numeric(1L)))
  W1pm <- rowMeans(vapply(res, \(x) x[[2]][[2]],
                          FUN.VALUE=numeric(length(unique(hardPred)))))
  W2pm <- rowMeans(vapply(res, \(x) x[[3]][[2]], numeric(nrow(fuzzyTrue))))
                          
  
  AW1 <- list( global=adj(val[[2]]$global,W1m),
               perPartition=mapply(x=val[[2]]$perPartition, m=W1pm, FUN=adj) )
  AW2 <- list( global=adj(val[[3]]$global,W2m),
               perPartition=mapply(x=val[[3]]$perPartition, m=W2pm, FUN=adj) )
  
  return(list(NDC=val[[1]], ACI=ACI, fuzzyWH=val[[2]], fuzzyWC=val[[3]],
              fuzzyAWH=AW1, fuzzyAWC=AW2))
}


#' Per-element concordance between two fuzzy partitions
#' 
#' @param P A object coercible to a numeric matrix with membership probability 
#'   of elements (rows) in clusters (columns)
#' @param Q A object coercible to a numeric matrix with membership probability 
#'   of elements (rows) in clusters (columns). Must have the same number of rows
#'   as `P`
#'
#' @return A numeric vector of concordance scores for each row of `P`.
#' @export
#' @examples
#' # generate fuzzy partitions:
#' m1 <- matrix(c(0.95, 0.025, 0.025, 
#'                0.98, 0.01, 0.01, 
#'                0.96, 0.02, 0.02, 
#'                0.95, 0.04, 0.01, 
#'                0.95, 0.01, 0.04, 
#'                0.99, 0.005, 0.005, 
#'                0.025, 0.95, 0.025, 
#'                0.97, 0.02, 0.01, 
#'                0.025, 0.025, 0.95), 
#'                ncol = 3, byrow=TRUE)
#' m2 <- matrix(c(0.95, 0.025, 0.025,  
#'                0.98, 0.01, 0.01, 
#'                0.96, 0.02, 0.02, 
#'                0.025, 0.95, 0.025, 
#'                0.02, 0.96, 0.02, 
#'                0.01, 0.98, 0.01, 
#'                0.05, 0.05, 0.95, 
#'                0.02, 0.02, 0.96, 
#'                0.01, 0.01, 0.98), 
#'                ncol = 3, byrow=TRUE)
#' colnames(m1) <- colnames(m2) <- LETTERS[seq_len(3)]
#' fuzzySpotConcordance(m1,m2)
fuzzySpotConcordance <- function(P, Q){
  if(is.data.frame(P)) P <- as.matrix(P)
  if(is.data.frame(Q)) Q <- as.matrix(Q)
  stopifnot(is.matrix(P) && (is.numeric(P) | is.integer(P)))
  stopifnot(is.matrix(Q) && (is.numeric(Q) | is.integer(Q)))
  stopifnot(nrow(P)==nrow(Q))
  
  diff <- as.matrix(0.5*abs(dist(P,method="manhattan") -
                              dist(Q,method="manhattan")))
  return(1-rowSums(diff/(ncol(diff)-1)))
}


#' Per-element maximal concordance between a hard and a fuzzy partition
#'
#' Per-element maximal concordance between a hard clustering and hard and fuzzy 
#'   ground truth labels.
#' 
#' @param hardPred A vector of predicted cluster labels
#' @param hardTrue A vector of true cluster labels
#' @param fuzzyTrue A object coercible to a numeric matrix with membership 
#'   probability of elements (rows) in clusters (columns). Must have the same 
#'   number of rows as the length of `hardTrue`.
#' @param useNegatives Logical; whether to include negative pairs in the 
#'   concordance score (tends to result in a larger overall concordance and 
#'   lower dynamic range of the score). Default TRUE.
#' @param verbose Logical; whether to print expected memory usage for large 
#'   datasets.
#'
#' @return A numeric vector of concordance scores for each element of `hardPred`
#' @export
#' @examples
#' # generate a fuzzy truth:
#' fuzzyTrue <- matrix(c(
#'   0.95, 0.025, 0.025, 
#'   0.98, 0.01, 0.01, 
#'   0.96, 0.02, 0.02, 
#'   0.95, 0.04, 0.01, 
#'   0.95, 0.01, 0.04, 
#'   0.99, 0.005, 0.005, 
#'   0.025, 0.95, 0.025, 
#'   0.97, 0.02, 0.01, 
#'   0.025, 0.025, 0.95), 
#'   ncol = 3, byrow=TRUE)
#' # a hard truth:
#' hardTrue <- apply(fuzzyTrue,1,FUN=which.max)
#' # some predicted labels:
#' hardPred <- c(1,1,1,1,1,1,2,2,2)
#' fuzzyHardSpotConcordance(hardTrue, fuzzyTrue, hardPred)
fuzzyHardSpotConcordance <- function(hardTrue, fuzzyTrue, hardPred, 
                                   useNegatives=TRUE, verbose=TRUE){
  if(useNegatives) return(fuzzyHardMetrics2(hardTrue, fuzzyTrue, hardPred,
                                            returnElementPairAccuracy=TRUE,
                                            verbose=verbose))
  stopifnot(is.atomic(hardPred))
  hardPredVector <- hardPred <- as.integer(as.factor(hardPred))
  if(is.atomic(hardTrue)){
    hardTrueVector <- hardTrue <- as.integer(as.factor(hardTrue))
    hardTrue <- Matrix::sparseMatrix(seq_along(hardTrue), hardTrue, x=1L)
  }else{
    if(is.data.frame(hardTrue)) hardTrue <- as.matrix(hardTrue)
    stopifnot(is.matrix(hardTrue) && 
                (is.numeric(hardTrue) | is.integer(hardTrue)))
    hardTrueVector <- apply(hardTrue, 1, FUN=which.max)
  }
  hardPred <- Matrix::sparseMatrix(seq_along(hardPred), hardPred, x=1L)
  
  if(is.data.frame(fuzzyTrue)) fuzzyTrue <- as.matrix(fuzzyTrue)
  
  stopifnot(is.matrix(fuzzyTrue) && 
              (is.numeric(fuzzyTrue) | is.integer(fuzzyTrue)))
  stopifnot(nrow(hardTrue)==nrow(fuzzyTrue))
  stopifnot(nrow(hardPred)==nrow(hardTrue))
  stopifnot(ncol(hardTrue)==ncol(fuzzyTrue))
  
  m <- nrow(fuzzyTrue)
  if(verbose && m>=2000){
    os <- 8 * m^2 * 4
    class(os) <- "object_size"
    message("Projected memory usage: ", format(os, units = "auto"))
  }
  
  # compute pair concordance for all the three labelings
  eq <- as.matrix(1-(0.5*dist(hardPred, method="manhattan")))
  ep <- 1-(0.5*dist(fuzzyTrue, method="manhattan"))
  ep2 <- 1-(0.5*dist(hardTrue, method="manhattan"))
  
  # compute the minimum difference of the pred with either truth
  diff <- pmin( abs(as.matrix(ep) - eq),
                abs(as.matrix(ep2) - eq) )
  rm(eq, ep, ep2)
  
  # if(useNegatives){
  #   return(1-rowSums(diff)/(ncol(diff)-1))
  # }
  
  # blend out the negative pairs
  1-vapply(seq_len(nrow(diff)), FUN.VALUE=numeric(1), FUN=function(i){
    w <- which(hardPredVector==hardPredVector[[i]] |
                 hardTrueVector==hardTrueVector[[i]])
    sum(diff[i,w])/(length(w)-1)
  })
  
}
