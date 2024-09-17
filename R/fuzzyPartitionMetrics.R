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
#'   \item fuzzyWallace1 Fuzzy Wallace index for each partition of `Q`
#'   \item fuzzyWallace2 Fuzzy Wallace index for each partition of `P`
#'   \item fuzzyAdjW1 Adjusted fuzzy Wallace index for each partition of `Q`
#'   \item fuzzyAdjW2 Adjusted fuzzy Wallace index for each partition of `P`
#'   
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
                                  BPPARAM=BiocParallel::SerialParam()){ 
  
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
  NDC <- 1 - ( sum(abs(ep[lower.tri(ep)] - eq[lower.tri(eq)]) )/(ncomp) )
  
  # precompute the pairs' class membership (for increased speed in permutations)
  Ppairs <- apply(P, 2, FUN=function(p) tcrossprod(p)[lower.tri(ep)])
  
  getFWallace <- function(emA, emB, fuzzyClB, Bpairs=NULL){
    a <- sapply(setNames(seq_len(ncol(fuzzyClB)),colnames(fuzzyClB)),
                FUN=function(i){
      # get the degree to which the members of each pair are of the given 
      # class in B
      if(is.null(Bpairs)){
        Bpair <- as.matrix(tcrossprod(fuzzyClB[,i])[lower.tri(emA)])
      }else{
        Bpair <- Bpairs[,i]
      }
      # compute 1 - the distance in A between elements weighted by their being 
      # of the same class in B
      Btot <- sum(Bpair)
      c(c=sum(emA[lower.tri(emA)]*Bpair), n=Btot)
    })
    list( global=sum(a[1,])/sum(a[2,]),
          perPartition=a[1,]/a[2,] )
  }
  
  if(computeWallace){
    W1 <- getFWallace(ep,eq,Q)
    W2 <- getFWallace(eq,ep,P,Bpairs=Ppairs)
  }
  
  
  # get metrics for the permutations
  
  # fn for one permutation:
  onePerm <- function(col){
    p <- allp[,col]
    permutedEQ <- eq[p,p]
    NDC <- 1 - ( sum(abs(ep-permutedEQ) )/(m*((m-1))) )
    if(!computeWallace) return(NDC)
    W1 <- getFWallace(ep,permutedEQ,Q[p,])
    W2 <- getFWallace(permutedEQ,ep,P,Bpairs=Ppairs)
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
    if(SE>0.0025)
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

  return(list(NDC=NDC, ACI=ACI, fuzzyWallace1=W1, fuzzyWallace2=W2,
              fuzzyAdjW1=AW1, fuzzyAdjW2=AW2))
}
