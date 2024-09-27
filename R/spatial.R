# Functions to turn neighborhood class distribution into fuzzy clusterings


#' findSpatialKNN
#' @param location A numeric data matrix containing location information, where 
#' rows are points and columns are dimensions.
#' @param k The number of nearest neighbors to look at.
#' @param keep_ties A Boolean indicating if ties of neighbors are counted once
#'  or not. If TRUE, neighbors of the same distances will be counted only once,
#'  and the resulting KNN will be more than k if there are any ties exist.
#' @param n 
#' @param BNPARAM 
#' @param ... 
#'
#' @description For a dataset, return the indices of knn for each object
findSpatialKNN <- function(location, k, keep_ties=TRUE, n=5, BNPARAM=NULL, ...){
  BNPARAM <- .decideBNPARAM(nrow(location), BNPARAM)
  if(keep_ties){
    nn <- BiocNeighbors::findKNN(location, k=k*n, warn.ties=FALSE, BNPARAM=BNPARAM)
    nn <- lapply(seq_len(nrow(nn[[1]])), FUN=function(i){
      d <- nn$distance[i,]
      nn$index[i,sort(which(d<=d[order(d)[k]]))]
    })
  }else{
    nn <- BiocNeighbors::findKNN(location, k=k, warn.ties=FALSE, BNPARAM=BNPARAM)$index
    nn <- split(nn, seq_len(nrow(nn)))
  }
  return(nn)
}

#' @alpha the parameter to control to what extend the spot itself contribute 
#' to the class composition calculation. "equal" means it is weighted the 
#' same as other NNs. A numeric value between 0 and 1 means the weight of the 
#' frequency contribution for the spot itself, and the frequency contribution 
#' for its knn is then 1-alpha.
knnComposition <- function(location, k=6, label, alpha="equal", ...){
  label <- factor(label)
  ind <- findSpatialKNN(location, k, ...)
  knnLabels <- lapply(ind, function(x){label[x[2:length(x)]]})
  if(alpha=="equal"){ 
    alpha <- 1/(k+1) 
  }else{
    if(!(is.numeric(alpha) & alpha<=1 & alpha>=0)){
      stop("alpha must be either 'equal', or a numeric between 0 and 1.")
    }
  }
  knn_weights <- lapply(knnLabels, function(x){x<-factor(x, levels=levels(label)); as.vector(table(x)/length(x)) * (1-alpha)})
  knn_weights <- do.call(rbind, knn_weights)
  i_weights <-  as.matrix(table(seq_along(label), label)) * (alpha)
  return(knn_weights + i_weights)
}

getFuzzyLabel <- function(label, location, k=6, alpha="equal", ...){
  label <- factor(label)
  NAs <- which(is.na(label))
  if(length(NAs>0)){
    label <- label[-NAs]
    location <- location[-NAs,]
  }
  res <- knnComposition(location=location, k=k, label=label, alpha=alpha, ...)
  return(res)
}

library(clue)
library(mclust)
library(RcppHungarian)


getPredLabels <- function(ref_labels, pred_clusters) {
  cost_matrix <- .computeCostMatrix(ref_labels, pred_clusters)
  cluster_map <- .getClusterMapping(cost_matrix)
  
  pred_labels <- unlist(cluster_map[pred_clusters])
  names(pred_labels) <- NULL
  
  return(pred_labels)
}

.computeCostMatrix <- function(ref_labels, pred_clusters) {
  # Create a matrix to store the cost
  unique_ref_labels <- unique(ref_labels)
  unique_pred_clusters <- unique(pred_clusters)
  count_matrix <- matrix(0, nrow = length(unique_ref_labels), ncol = length(unique_pred_clusters), dimnames = list(unique_ref_labels, unique_pred_clusters))
  
  # Iterate over the indices and update the matrix
  for (i in seq_along(ref_labels)) {
    count_matrix[ref_labels[i], pred_clusters[i]] <- count_matrix[ref_labels[i], pred_clusters[i]] + 1
  }
  
  if (ncol(count_matrix) > 1) {
    cost_matrix <- apply(count_matrix, 1, function(row) max(row) - row)
  } else {
    cost_matrix <- t(max(count_matrix) - count_matrix)
  }
  
  return(cost_matrix)
}

.getClusterMapping <- function(cost_matrix) {
  solved <- HungarianSolver(cost_matrix)
  
  cluster_map <- list()
  for (i in 1:nrow(solved$pairs)) {
    from <- rownames(cost_matrix)[solved$pairs[i, 1]]
    to <- if(solved$pairs[i, 2] == 0) from else colnames(cost_matrix)[solved$pairs[i, 2]] 
    cluster_map[[from]] <- to
  }
  
  return(cluster_map)
}

#' matchSets
#' 
#' Match sets from a partitions to a reference partition using the Hungarian
#' algorithm.
#'
#' @param pred An integer or factor of cluster labels
#' @param true An integer or factor of reference labels
#' @param forceMatch Logical; whether to enforce a match for every set of `pred`
#' @param returnIndices Logical; whether to return indices rather than levels
#'
#' @return A vector of matching sets (i.e. level) from `true` for every set 
#'   (i.e. level) of `pred`.
#' @importFrom clue solve_LSAP
matchSets <- function(pred, true, forceMatch=TRUE, returnIndices=is.integer(true)){
  true <- as.factor(true)
  pred <- as.factor(pred)
  co <- unclass(table(true, pred))
  recall <- co/rowSums(co)
  prec <- t(t(co)/colSums(co))
  F1 <- 2*(prec*recall)/(prec+recall)
  F1[is.na(F1)] <- 0
  if(nrow(F1)>ncol(F1)){
    match <- as.integer(clue::solve_LSAP(t(F1), maximum=TRUE))
  }else{
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