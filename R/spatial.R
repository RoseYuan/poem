# Functions to turn neighborhood class distribution into fuzzy clusterings

#' @description For the i-th object, return the indices of its knn
findSpatialKNN <- function(i, location_in, k){
  require(pdist)
  line_i <- rep(0,dim(location_in)[1])
  line_i <- pdist(location_in[i,],location_in)@dist
  ind <- order(line_i)[1:(k+1)]
  return(ind)
}

#' @alpha the parameter to control to what extend the spot itself contribute 
#' to the class composition calculation. "equal" means it is weighted the 
#' same as other NNs. A numeric value between 0 and 1 means the weight of the 
#' frequency contribution for the spot itself, and the frequency contribution 
#' for its knn is then 1-alpha.
knnComposition <- function(i, location, k=6, label, alpha="equal"){
  label <- factor(label)
  ind <- findSpatialKNN(i, location, k)
  knnLabels <- label[ind[2:length(ind)]]
  if(alpha=="equal"){ 
    alpha <- 1/(k+1) 
  }else{
    if(!(is.numeric(alpha) & alpha<=1 & alpha>=0)){
      stop("alpha must be either 'equal', or a numeric between 0 and 1.")
    }
  }
  knn_weights <- matrix(table(knnLabels)/k) * (1-alpha)
  i_weights <-  matrix(table(label[i])) * (alpha)
  return(knn_weights + i_weights)
}

getFuzzyLabel <- function(label, location, k=6, alpha="equal"){
  require(parallel)
  label <- factor(label)
  NAs <- which(is.na(label))
  if(length(NAs>0)){
    label <- label[-NAs]
    location <- location[-NAs,]
  }
  res <- mclapply(1:dim(location)[1], knnComposition, location=location, k=k, label=label, alpha=alpha, mc.cores = 5)
  return(do.call(cbind, res))
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
