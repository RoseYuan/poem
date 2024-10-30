#' @importFrom igraph subgraph vertex_attr is_directed
.igraphFunPerClass <- function(g, classAttr="class", directed=is_directed(g), 
                               FUN, ...){
  stopifnot(is(g,"igraph"))
  if(is.null(directed)) directed <- FALSE
  if(is.character(FUN)) FUN <- get(FUN)
  if(!is.null(classAttr)){
    va <- vertex_attr(g,classAttr)
    stopifnot(is.factor(va) || is.character(va) || is.integer(va))
    names(cls) <- cls <- unique(va)
    return(sapply(cls, FUN=function(cl){
      FUN(subgraph(g, vids=which(va==cl)), ...)
    }))
  }
  FUN(g, ...)
}

#' @importFrom igraph mean_distance decompose.graph is_directed V
.adjMeanShortestPath <- function(g, directed=FALSE, normalize=TRUE){
  stopifnot(is(g,"igraph"))
  if(is.null(directed)) directed <- FALSE
  gc <- decompose.graph(g)
  msp <- sapply(gc, directed=directed, FUN=mean_distance)
  amsp <- length(gc)+sum(pmax(1L,msp,na.rm=TRUE))
  if(normalize) amsp <- amsp/length(V(g))
  return(amsp)
}

#' @importFrom utils relist
#' @importFrom igraph vertex_attr
.simpsonIndex <- function(knn, labels=NULL, directed=TRUE){
  if(is.null(directed)) directed <- TRUE
  # for undirected, first transform directed nn to graph
  if(!is(knn,"igraph") && !directed) knn <- .nn2graph(knn, labels)
  if(is(knn,"igraph")){
    if(is.null(labels)) labels <- vertex_attr(knn, "class")
    knn <- .igraph2nn(knn, labels, directed=directed)
  }
  stopifnot(!is.null(labels))
  if(.isNNlist(knn)){ # list of varying number of neighbor indices
    knn <- relist(labels[unlist(knn)],knn)
    return(sapply(knn, FUN=function(nn){
        sum((sapply(unique(labels), FUN=function(x) sum(nn==x))/length(nn))^2)
      }))
  }
  .checkInputs(knn, labels)
  k <- ncol(knn$index)
  p <- sapply(unique(labels), FUN=function(x){
    rowSums(knn$nncl==x)/k
  })
  rowSums(p^2)
}

# neighborhood purity (i.e. proportion with same labels as target node)
#' @importFrom igraph as_adj_list vertex_attr
.nPurity <- function(knn, labels=NULL, directed=TRUE){
  if(is.null(directed)) directed <- TRUE
  # for undirected, first transform directed nn to graph
  if(!is(knn,"igraph") && !directed) knn <- .nn2graph(knn, labels)
  if(is(knn,"igraph")){
    if(is.null(labels)) labels <- vertex_attr(knn, "class")
    knn <- .igraph2nn(knn, labels, directed=directed)
  }
  stopifnot(!is.null(labels))
  if(.isNNlist(knn)){ # list of varying number of neighbor indices
    return(mapply(nn=relist(labels[unlist(knn)],knn),
                  label=labels, FUN=function(nn,label){
      sum(nn==label)/length(nn)
    }))
  }
  .checkInputs(knn, labels)
  rowSums(knn$nncl==labels)/ncol(knn$index)
}

#' @importFrom igraph vertex_attr
# neighborhood class over-representation
.nlog2Enrichment <- function(knn, labels, directed=TRUE, pseudoCount=1){
  if(is.null(directed)) directed <- TRUE
  # for undirected, first transform directed nn to graph
  if(!is(knn,"igraph") && !directed) knn <- .nn2graph(knn, labels)
  if(is(knn,"igraph")){
    if(is.null(labels)) labels <- vertex_attr(knn, "class")
    knn <- .igraph2nn(knn, labels, directed=directed)
  }
  stopifnot(!is.null(labels))
  if(.isNNlist(knn)){ # list of varying number of neighbor indices
    expect <- table(labels)/length(labels)
    return(mapply(nn=relist(labels[unlist(knn)],knn),
                  label=labels, FUN=function(nn,label){
                    k <- ncol(knn$index)
                    expected <- length(nn)*expect
                    log2((pseudoCount+sum(nn==label))/
                           (pseudoCount+length(nn)*expect))
                  }))
  }
  .checkInputs(knn, labels)
  k <- ncol(knn$index)
  expected <- k*as.numeric(table(labels))/length(labels)
  log2((pseudoCount+rowSums(knn$nncl==labels))/
         (pseudoCount+expected[as.integer(labels)]))
}

#' Calculate F measure
#' 
#' Compute the F measure between two clustering results. This is directly copied
#' from the package `FlowSOM`.
#'
#' @param true Array containing real cluster labels for each sample
#' @param pred Array containing predicted cluster labels for each
#'                          sample
#' @param silent    Logical, if FALSE, print some information about 
#'                  precision and recall
#' @keywords internal
#' @return  F measure score
FMeasure <- function(true, pred, silent=TRUE){
  if (sum(pred)==0)
    return(0);
  a <- table(true, pred);
  p <- t(apply(a,1,function(x)x/colSums(a)))
  r <- apply(a,2,function(x)x/rowSums(a))
  if(!silent) message("Precision: ",
                      sum(apply(p,1,max) * (rowSums(a)/sum(a))),
                      "\nRecall: ",sum(apply(r,1,max) * (rowSums(a)/sum(a))),"\n")
  f <- 2*r*p / (r+p)
  f[is.na(f)] <- 0
  sum(apply(f,1,max) * (rowSums(a)/sum(a)))
}

.NCR <- function(true, pred){
  co <- table(pred, true)
  a <- mean(colSums(choose(co,2))/choose(table(true),2),na.rm=TRUE)
  dmax <- tcrossprod(table(true))
  d <- (dmax-crossprod(co))/dmax
  d <- mean(d[lower.tri(d)])
  (a+d)/2
}


#' Per-element agreement score
#' 
#' Per-element agreement between a clustering and a ground truth
#'
#' @param pred A vector of predicted clusters
#' @param true A vector of true class labels
#' @param usePairs Logical; whether to compute over pairs instead of elements 
#'   Recommended and TRUE by default.
#' @param useNegatives Logical; whether to include the consistency of negative
#'   pairs in the score (default FALSE).
#' @param adjust Logical; whether to adjust for chance. Only implemented for
#'   `useNegatives=FALSE` (doesn't make sense on a element-level otherwise).
#'
#' @return A vector of agreement scores
getAgreement <- function(true, pred, usePairs=TRUE, useNegatives=FALSE, adjust=FALSE){
  if(useNegatives & adjust)
    stop("Adjustment for chance is only implemented for useNegatives=FALSE")
  # contingency table
  co <- table(true, pred)
  pairs <- choose(co,2)
  # best possible contingency table
  bco <- diag(table(true))
  bpairs <- choose(bco,2)
  
  # number of spots in the union between any class and any cluster:
  tot <- matrix(rep(rowSums(co),ncol(co)),nrow=nrow(co))+
    matrix(rep(colSums(co),each=nrow(co)),nrow=nrow(co))-co
  
  if(usePairs){
    if(useNegatives){
      # if(adjust){
      #   allp <- apply(matrix(runif(length(pred) * 20), nrow=length(pred)), 2, order)
      #   return(sapply(seq_along(pred), FUN=function(i){
      #     true.co <- (true==true[[i]])
      #     pa <- 1-sum(abs( (pred==pred[[i]]) - true.co ))/(length(pred)-1)
      #     rpa <- mean(apply(allp,2,FUN=function(x){
      #       pred <- pred[x]
      #       1-sum(abs( (pred==pred[[i]]) - true.co ))/(length(pred)-1)
      #     }))
      #     (pa - rpa)/(1 - rpa)
      #   }))
      # }

      return(sapply(seq_along(pred), FUN=function(i){
        1-sum(abs( (pred==pred[[i]]) - (true==true[[i]]) ))/(length(pred)-1)
      }))
      
      # # faster version on CM, but wrong at the moment:
      # p <- sapply(seq_len(ncol(co)),FUN=function(j){
      #   sapply(seq_len(nrow(co)), FUN=function(i){
      #     if(co[i,j]==0) return(0)
      #     pos <- choose(co[i,j],2)
      #     neg <- co[i,j]*sum(co[-i,-j])/2
      #     max <- pos+(co[i,j]*(nrow(co)-co[i,j]))/2
      #     (pos+neg)/max
      #   })
      # })
      # dimnames(p) <- dimnames(co)
    }else{
      truePairsPerCluster <- matrix(rep(colSums(pairs),each=nrow(co)),nrow=nrow(co))
      truePairs <- truePairsPerCluster + #per cluster
        rowSums(pairs) - pairs # per class minus double-counting
      p <- unclass(truePairs/choose(tot,2))
      if(adjust){
        # expected random contingency table
        rco <- as.numeric(table(true)/length(true)) %*% t(table(pred)/length(pred))
        rco <- rco*length(true)
        # number of spots in the expected union between any class and any cluster:
        rtot <- matrix(rep(rowSums(rco),ncol(rco)),nrow=nrow(rco))+
          matrix(rep(colSums(rco),each=nrow(rco)),nrow=nrow(rco))-rco
        pairs <- choose(rco,2)
        truePairsPerCluster <- matrix(rep(colSums(pairs),each=nrow(rco)),nrow=nrow(rco))
        truePairs <- truePairsPerCluster + #per cluster
          rowSums(pairs) - pairs # per class minus double-counting
        rp <- unclass(truePairs/choose(rtot,2))
        p <- (p - rp)/(1 - rp)
      }
    }
  }else{
    if(adjust || useNegatives) stop("Only implemented for usePairs=TRUE")
    # intersection over union, i.e. proportion of spots in the class-or-cluster
    # that agree:
    p <- unclass(co/tot)
  }
  # assign each spot its score:
  p <- setNames(as.numeric(p), paste(rep(row.names(p),ncol(p)),rep(colnames(p),each=nrow(p))))
  as.numeric(p[paste(true, pred)])
}

#' The non-spatially-weighted counterpart of nnWeightedAccuracy
#' 
#' @param true True class labels (vector coercible to factor)
#' @param pred Predicted labels (vector coercible to factor)
#' @keywords internal
#' @return A scalar representing the weighted accuracy.
setMatchingAccuracy <- function(true, pred){
  pred <- as.factor(pred)
  matching <- matchSets(pred, true, returnIndices=TRUE)
  pred <- matching[as.integer(pred)]
  true <- as.integer(as.factor(true))
  sum(pred==true,na.rm=TRUE)/length(pred)
}
