

#' Find the k nearest spatial neighbors
#' @description For a given dataset, find the `k` nearest neighbors for each 
#' object based on their spatial locations, with the option of handling ties.
#' @param location A numeric data matrix containing location information, where 
#' rows are points and columns are location dimensions.
#' @param k The number of nearest neighbors to look at.
#' @param keep_ties A Boolean indicating if ties are counted once or not. If 
#'  TRUE, neighbors of the same distances will be included even if it means 
#'  returning more than `k` neighbors.
#' @param useMedianDist Use the median distance of the k nearest neighbor as
#'  maximum distance to be included. Ignored if `keep_ties=FALSE`.
#' @param BNPARAM BNPARAM object passed to BiocNeighbors::findKNN specifying the
#'  kNN approximation method to use. Defaults to exact for small datasets, and 
#'  Annoy for larger ones.
#' @param ... Ignored
#'
#' @return A list of indices.
#' @importFrom BiocNeighbors findKNN
#' @export
findSpatialKNN <- function(location, k, keep_ties=TRUE, useMedianDist=FALSE,
                           BNPARAM=NULL, ...){
  BNPARAM <- .decideBNPARAM(nrow(location), BNPARAM)
  if(keep_ties){
    nn <- BiocNeighbors::findKNN(location, k=k*3, warn.ties=FALSE, BNPARAM=BNPARAM)
    mkd <- median(nn$distance[,k])
    nn <- lapply(seq_len(nrow(nn[[1]])), FUN=function(i){
      d <- nn$distance[i,]
      if(!useMedianDist) mkd <- d[order(d)[k]]
      nn$index[i,sort(which(d<=mkd))]
    })
  }else{
    nn <- BiocNeighbors::findKNN(location, k=k, warn.ties=FALSE, BNPARAM=BNPARAM)$index
    nn <- split(nn, seq_len(nrow(nn)))
  }
  return(nn)
}

#' Compute neighborhood composition
#' @description For a given dataset with locations and labels, compute the label
#' composition of the neighborhood for each sample.
#' @param alpha The parameter to control to what extend the spot itself contribute 
#' to the class composition calculation. `"equal"` means it is weighted the 
#' same as other neighbors. A numeric value between `0` and `1` means the weight of the 
#' frequency contribution for the spot itself, and the frequency contribution 
#' for its knn is then `1-alpha`. By default `0.5`.
#' @inheritParams findSpatialKNN 
#' @param label A vector containing the label for the dataset.
#' @export
#' @return A numerical matrix indicating the composition, where rows correspond 
#' to samples and columns correspond to the classes in `label`. 
knnComposition <- function(location, k=6, label, alpha=0.5, ...){
  label <- factor(label)
  ind <- findSpatialKNN(location, k, ...)
  knnLabels <- lapply(ind, function(x){label[x[1:length(x)]]})
  if(alpha=="equal"){ 
    alpha <- 1/(k+1) 
  }else{
    if(!(is.numeric(alpha) & alpha<=1 & alpha>=0)){
      stop("alpha must be either 'equal', or a numeric between 0 and 1.")
    }
  }
  knn_weights <- lapply(knnLabels, function(x){x<-factor(x, levels=levels(label)); as.vector(table(x)/length(x)) * (1-alpha)})
  knn_weights <- do.call(rbind, knn_weights)
  i_weights <-  as.data.frame.matrix(table(seq_along(label), label)) * (alpha)
  return(knn_weights + i_weights)
}

#' Get fuzzy representation of labels
#' Get fuzzy representation of labels according to the spatial neighborhood 
#' label composition.
#' @param label An anomic vector of cluster labels
#' @param location A matrix or data.frame of coordinates
#' @param k The wished number of nearest neighbors
#' @inheritParams knnComposition
#' @param ... Passed to [findSpatialKNN()].
#'
#' @return A matrix of fuzzy memberships.
#' @export

getFuzzyLabel <- function(label, location, k=6, alpha=0.5, ...){
  label <- factor(label)
  stopifnot(!any(is.na(label)))
  res <- knnComposition(location=location, k=k, label=label, alpha=alpha, ...)
  return(res)
}

#' matchSets
#' 
#' Match sets from a partitions to a reference partition using the Hungarian
#' algorithm to optimize F1 scores.
#'
#' @param pred An integer or factor of cluster labels
#' @param true An integer or factor of reference labels
#' @param forceMatch Logical; whether to enforce a match for every set of `pred`
#' @param returnIndices Logical; whether to return indices rather than levels
#'
#' @return A vector of matching sets (i.e. level) from `true` for every set 
#'   (i.e. level) of `pred`.
#' @importFrom clue solve_LSAP
#' @export
matchSets <- function(pred, true, forceMatch=TRUE,
                      returnIndices=is.integer(true)){
  true <- as.factor(true)
  pred <- as.factor(pred)
  co <- unclass(table(true, pred))
  recall <- co/rowSums(co)
  prec <- t(t(co)/colSums(co))
  F1 <- 2*(prec*recall)/(prec+recall)
  F1[is.na(F1)] <- 0
  if(nrow(F1)>ncol(F1)){
    # using alternative dependencies:
    # match <- RcppHungarian::HungarianSolver(-t(F1))$pairs[,2]
    match <- as.integer(clue::solve_LSAP(t(F1), maximum=TRUE))
  }else{
    # using alternative dependencies:
    # match1 <- RcppHungarian::HungarianSolver(-F1)$pairs[,2]
    # match1[which(match1==0)] <- NA_integer_ 
    match1 <- as.integer(clue::solve_LSAP(F1, maximum=TRUE))
    match <- rep(NA_integer_, length(levels(pred)))
    nonMatched <- setdiff(seq_along(match), match1)
    match[match1] <- seq_along(match1)
    if(forceMatch){
      for(i in nonMatched) match[i] <- which.max(F1[,i])
    }
  }
  if(!returnIndices) match <- levels(true)[match]
  setNames(match, levels(pred))
}