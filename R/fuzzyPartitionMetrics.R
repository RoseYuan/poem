#' fuzzyPartitionMetrics
#' 
#' Computes fuzzy versions of pair-sorting partition metrics. This is largely 
#' based on the permutation-based implementation by Antonio D'Ambrosio from the 
#' ConsRankClass package, modified to also compute the fuzzy versions of the 
#' adjusted Wallace indices, implement multithreading, and adjust the number of
#' permutations according to their variability.
#'
#' @param P A object coercible to a numeric matrix with membership probability 
#'   of elements (rows) in clusters (columns)
#' @param Q A object coercible to a numeric matrix with membership probability 
#'   of elements (rows) in clusters (columns). Must have the same number of rows
#'   as `P`
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
#'   pair accuracy instead of the various parition-level and global metrics.
#'   Default FALSE.
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
#' @return A list of metrics:
#'   \item NDC : Hullermeier's NDC (fuzzy rand index)
#'   \item ACI : Ambrosio's Adjusted Concordance Index (ACI), i.e. a 
#'     permutation-based fuzzy version of the adjusted Rand index.
#'   \item fuzzyWH Fuzzy Wallace Homogeneity index
#'   \item fuzzyWC Fuzzy Wallace Completeness index
#'   \item fuzzyAWH Adjusted fuzzy Wallace Homogeneity index
#'   \item fuzzyAWC Adjusted fuzzy Wallace Completeness index
#' @importFrom BiocParallel SerialParam bplapply
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
#' colnames(m1) <- colnames(m2) <- LETTERS[1:3]
#' fuzzyPartitionMetrics(m1,m2)
fuzzyPartitionMetrics <- function(P, Q, computeWallace=TRUE, nperms=NULL,
                                  verbose=TRUE, returnElementPairAccuracy=FALSE,
                                  BPPARAM=BiocParallel::SerialParam(), 
                                  tnorm=c("product","min","lukasiewicz"), ...){ 
  
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
      os <- 8*(4*m^2 + m*ifelse(is.null(nperms),10,nperms))
      if(computeWallace) os <- os + 8*ncomp*(1+max(ncol(Q),ncol(P)))
    }
    class(os) = "object_size"
    message("Projected memory usage: ", format(os, units = "auto"))
  }
  
  ep <- as.matrix(1-(0.5*dist(P,method="manhattan")))
  eq <- as.matrix(1-(0.5*dist(Q,method="manhattan")))
  
  if(returnElementPairAccuracy){
    return(1-rowMeans(abs(ep - eq)))
  }
  
  # Hullermeier's NDC
  diff <- abs(ep[lower.tri(ep)] - eq[lower.tri(eq)])
  NDC <- 1 - ( sum( diff )/(ncomp) )
  
  membershipFn <- switch(tnorm,
    product=tcrossprod,
    lukasiewicz=function(p){ sapply(seq_along(p), FUN=function(i) pmax(0,p+p[i]-1)) },
    min=function(p){ sapply(seq_along(p), FUN=function(i) pmin(p,p[i])) }
  )
  # precompute the pairs' class membership (for increased speed in permutations)
  Ppairs <- apply(P, 2, FUN=function(p) membershipFn(p)[lower.tri(ep)])
  
  getFWallace <- function(emA, emB, fuzzyClB, Bpairs=NULL, diff=NULL){
    if(is.null(diff)) diff <- abs(emA[lower.tri(emA)]-emB[lower.tri(emB)])
    a <- sapply(setNames(seq_len(ncol(fuzzyClB)),colnames(fuzzyClB)),
                FUN=function(i){
      # get the degree to which the members of each pair are of the given 
      # class in B
      if(is.null(Bpairs)){
        Bpair <- as.matrix(membershipFn(fuzzyClB[,i])[lower.tri(emA)])
      }else{
        Bpair <- Bpairs[,i]
      }
      # compute 1 - the distance between pair agreements weighted by their being 
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
    res1 <- bplapply(1:10, BPPARAM=BPPARAM, onePerm)
    NDCs <- sapply(res1, \(x) x[[1]])
    SE <- sd(NDCs)/sqrt(length(res1))
    if(SE>0.0025) nperms <- 100
    if(SE>0.01) nperms <- 1000
    if(verbose && !is.null(nperms))
      message("Running ", nperms, " extra permutations.")
  }
  
  if(!is.null(nperms)){
    # generate more permutations
    allp <- apply(matrix( runif(m*nperms), nrow=m ), 2, order)
  
    res <- bplapply(1:nperms, BPPARAM=BPPARAM, onePerm)
    if(!is.null(res1)) res <- c(res1,res)
    NDCs <- sapply(res, \(x) x[[1]])
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
  if(!computeWallace) return(list(NDC=NDC, ACI=ACI))
  
  W1m <- mean(sapply(res, \(x) x[[2]][[1]]))
  W2m <- mean(sapply(res, \(x) x[[3]][[1]]))
  W1pm <- rowMeans(sapply(res, \(x) x[[2]][[2]]))
  W2pm <- rowMeans(sapply(res, \(x) x[[3]][[2]]))
  AW1 <- list( global=adj(W1$global,W1m),
               perPartition=mapply(x=W1$perPartition, m=W1pm, FUN=adj) )
  AW2 <- list( global=adj(W2$global,W2m),
               perPartition=mapply(x=W2$perPartition, m=W2pm, FUN=adj) )

  return(list(NDC=NDC, ACI=ACI, fuzzyWH=W1, fuzzyWC=W2,
              fuzzyAWH=AW1, fuzzyAWC=AW2))
}



#' fuzzyHardMetrics
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
#' @param hardTruth An atomic vector coercible to a factor or integer vector 
#'  containing the true hard labels. Must have the same length as `hardPred`.
#' @param fuzzyTruth A object coercible to a numeric matrix with membership 
#'   probability of elements (rows) in clusters (columns). Must have the same 
#'   number of rows as the length of `hardTruth`.
#' @param nperms The number of permutations (for correction for chance). If 
#'   NULL (default), a first set of 10 permutations will be run to estimate 
#'   whether the variation across permutations is above 0.0025, in which case 
#'   more (max 1000) permutations will be run.
#' @param verbose Logical; whether to print info and warnings, including the 
#'   standard error of the mean across permutations (giving an idea of the 
#'   precision of the adjusted metrics).
#' @param BPPARAM BiocParallel params for multithreading (default none)
#' 
#' @references Hullermeier et al. 2012; 10.1109/TFUZZ.2011.2179303;
#' @references D'Ambrosio et al. 2021; 10.1007/s00357-020-09367-0
#' 
#' @seealso [fuzzyPartitionMetrics()]
#' 
#' @author Pierre-Luc Germain
#'
#' @return A list of metrics:
#'   \item NDC : Hullermeier's NDC (fuzzy rand index)
#'   \item ACI : Ambrosio's Adjusted Concordance Index (ACI), i.e. a 
#'     permutation-based fuzzy version of the adjusted Rand index.
#'   \item fuzzyWH Fuzzy Wallace Homogeneity index
#'   \item fuzzyWC Fuzzy Wallace Completeness index
#'   \item fuzzyAWH Adjusted fuzzy Wallace Homogeneity index
#'   \item fuzzyAWC Adjusted fuzzy Wallace Completeness index
#' @importFrom BiocParallel SerialParam bplapply bpnworkers
#' @importFrom stats dist
#' @importFrom Matrix sparseMatrix
#' @export
#' @examples
#' # generate a fuzzy truth:
#' fuzzyTruth <- matrix(c(
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
#' hardTruth <- apply(fuzzyTruth,1,FUN=which.max)
#' # some predicted labels:
#' hardPred <- c(1,1,1,1,1,1,2,2,2)
#' fuzzyHardMetrics(hardPred, hardTruth, fuzzyTruth, nperms=3)
fuzzyHardMetrics <- function(hardPred, hardTruth, fuzzyTruth, nperms=NULL, 
                             verbose=TRUE, BPPARAM=BiocParallel::SerialParam(),
                             ...){ 
  
  stopifnot(is.atomic(hardPred))
  hardPredVector <- hardPred <- as.integer(as.factor(hardPred))
  if(is.atomic(hardTruth)){
    hardTruthVector <- hardTruth <- as.integer(as.factor(hardTruth))
    hardTruth <- Matrix::sparseMatrix(seq_along(hardTruth), hardTruth, x=1L)
  }else{
    if(is.data.frame(hardTruth)) hardTruth <- as.matrix(hardTruth)
    stopifnot(is.matrix(hardTruth) && 
                (is.numeric(hardTruth) | is.integer(hardTruth)))
    hardTruthVector <- apply(hardTruth, 1, FUN=which.max)
  }
  hardPred <- Matrix::sparseMatrix(seq_along(hardPred), hardPred, x=1L)
  
  if(is.data.frame(fuzzyTruth)) fuzzyTruth <- as.matrix(fuzzyTruth)
  
  stopifnot(is.matrix(fuzzyTruth) && 
              (is.numeric(fuzzyTruth) | is.integer(fuzzyTruth)))
  stopifnot(nrow(hardTruth)==nrow(fuzzyTruth))
  stopifnot(nrow(hardPred)==nrow(hardTruth))
  stopifnot(ncol(hardTruth)==ncol(fuzzyTruth))
  
  m <- nrow(fuzzyTruth)
  ncomp <- m*((m-1)/2)
  if(verbose && m>=2000){
    nSim <- min(BiocParallel::bpnworkers(BPPARAM),
                ifelse(is.null(nperms),10,nperms))
    os <- 8*((3+nSim)*m^2)
    class(os) = "object_size"
    message("Projected memory usage: ", format(os, units = "auto"))
  }
  
  # compute pair agreement for all the three labelings
  eq <- as.matrix(1-(0.5*dist(hardPred, method="manhattan")))
  ep <- 1-(0.5*dist(fuzzyTruth, method="manhattan"))
  ep2 <- 1-(0.5*dist(hardTruth, method="manhattan"))
  
  # compute the minimum difference of the pred with either truth
  diff <- pmin( abs(as.matrix(ep) - eq),
                abs(as.matrix(ep2) - eq) )
  NDC <- 1 - ( sum( diff )/(2*ncomp) )
  
  getFWallace <- function(hardTruthVector, diff){
    a <- sapply(split(seq_along(hardTruthVector),hardTruthVector),
                FUN=function(i){
      c(c=sum(diff[i,]), n=length(i)*ncol(diff))
    })
    a[2,which(a[2,]==0)] <- 1  # avoid NaNs for singletons
    list( global=1-sum(a[1,])/sum(a[2,], na.rm=TRUE),
          perPartition=1-a[1,]/a[2,] )
  }
  
  W1 <- getFWallace(hardPredVector, diff)
  W2 <- getFWallace(hardTruthVector, diff)
  # get metrics for the permutations
  
  # fn for one permutation:
  onePerm <- function(col){
    p <- allp[,col]
    permutedEQ <- eq[p,p]
    diff <- pmin( abs(as.matrix(ep) - permutedEQ),
                  abs(as.matrix(ep2) - permutedEQ) )
    NDC <- 1 -  sum( diff )/(2*ncomp)
    W1 <- getFWallace(hardPredVector[p], diff)
    W2 <- getFWallace(hardTruthVector, diff)
    # rm(diff, permutedEQ)
    # gc(verbose=FALSE)
    list(NDC=NDC, W1=W1, W2=W2)
  }
  
  res1 <- NULL
  if(is.null(nperms)){
    # try few permutations to estimate if more are needed
    allp <- apply(matrix( runif(m*10), nrow=m ), 2, order)
    res1 <- bplapply(1:10, BPPARAM=BPPARAM, onePerm)
    NDCs <- sapply(res1, \(x) x[[1]])
    SE <- sd(NDCs)/sqrt(length(res1))
    if(SE>0.0025) nperms <- 100
    if(SE>0.01) nperms <- 1000
    if(verbose && !is.null(nperms))
      message("Running ", nperms, " extra permutations.")
  }
  
  if(!is.null(nperms)){
    # generate more permutations
    allp <- apply(matrix( runif(m*nperms), nrow=m ), 2, order)
    
    res <- bplapply(1:nperms, BPPARAM=BPPARAM, onePerm)
    if(!is.null(res1)) res <- c(res1,res)
    NDCs <- sapply(res, \(x) x[[1]])
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
  
  W1m <- mean(sapply(res, \(x) x[[2]][[1]]))
  W2m <- mean(sapply(res, \(x) x[[3]][[1]]))
  W1pm <- rowMeans(sapply(res, \(x) x[[2]][[2]]))
  W2pm <- rowMeans(sapply(res, \(x) x[[3]][[2]]))
  AW1 <- list( global=adj(W1$global,W1m),
               perPartition=mapply(x=W1$perPartition, m=W1pm, FUN=adj) )
  AW2 <- list( global=adj(W2$global,W2m),
               perPartition=mapply(x=W2$perPartition, m=W2pm, FUN=adj) )
  
  return(list(NDC=NDC, ACI=ACI, fuzzyWH=W1, fuzzyWC=W2,
              fuzzyAWH=AW1, fuzzyAWC=AW2))
}
