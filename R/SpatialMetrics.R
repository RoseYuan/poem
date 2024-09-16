#' nnWeightedAccuracy
#' 
#' Computes an accuracy score which weighs elements/spots that are misclassified
#' by the proportion of their (spatial) neighborhood that is not of the 
#' element/spot's predicted class. This reduces the weight of misclassifications
#' happening at the boundary between domains.
#'
#' @param location The spatial coordinates to compute the nearest neighbors.
#' @param truth True class labels (vector coercible to factor)
#' @param pred Predicted labels (vector coercible to factor)
#' @param k Number of nearest neighbors
#'
#' @return A scalar representing the weighted accuracy.
nnWeightedAccuracy <- function(location, truth, pred, k=5, ...){
  pred <- as.factor(pred)
  levels(pred) <- matchSets(pred, truth, returnIndices=TRUE)
  pred <- as.integer(pred)
  truth <- as.integer(as.factor(truth))
  nn <- findSpatialKNN(location, k, ...)
  nn_dis <- relist(truth[unlist(nn)]!=rep(pred[w],lengths(nn)),nn)
  1-sum(sapply(nn_dis, mean))/nrow(location)
}

#' @title Calculate PAS score to measure clustering performance.
#' @description PAS score measures the randomness of the spots that located outside of the spatial region where it was clustered to.
#' Lower PAS score indicates better spatial domian clustering performance.
#' @param label Cluster labels.
#' @param location A n by k matrix of spatial locations.
#' @param k size of the neighborhood.
#' @return A numeric value for PAS score, and a boolean vector about the abnormal spots.
#' @export
PAS <- function(location, label, k=10, ...){
  matched_location=location
  NAs = which(is.na(label))
  if(length(NAs>0)){
    label=label[-NAs]
    matched_location = matched_location[-NAs,]
  }
  comp <- knnComposition(matched_location, k=k, label, alpha=0, ...)
  prop <- unlist(lapply(seq_along(1:dim(comp)[1]), function(i){comp[i,label[i]]}))
  results <- prop < 0.5
  return(list(PAS=sum(results)/length(results), abnormalty=results))
}

#' @title Calculate CHAOS score to measure clustering performance.
#' @description CHAOS score measures the mean length of the graph edges in the 
#' 1-nearest neighbor (1NN) graph for each cluster, averaged across clusters.
#' Lower CHAOS score indicates better spatial domain clustering performance.
#' @param location A n by k matrix of spatial locations.
#' @param label Cluster labels.
#' @param BNPARAM 
#' @return A numeric value for CHAOS score.
#' @export
CHAOS <- function(location, label, BNPARAM=NULL) {
  BNPARAM <- .decideBNPARAM(nrow(location), BNPARAM)
  label <- as.vector(label)
  location <- as.matrix(location)
  matched_location <- scale(location)
  label_unique <- unique(label)
  # Initialize a vector to hold the distance values
  dist_val <- numeric(length(label_unique))
  
  count <- 1
  for (k in label_unique) {
    # Get the locations belonging to the current cluster
    location_cluster <- matched_location[label == k, ]
    # Skip clusters with fewer than or equal to 2 points
    if (nrow(location_cluster) <= 2) {
      next
    }
    # Find the nearest neighbors within the cluster using BiocNeighbors::findKNN
    knn_result <- BiocNeighbors::findKNN(location_cluster, k = 1, warn.ties=FALSE, get.index = FALSE, BNPARAM=BNPARAM)
    # The distances to the nearest neighbors are stored in the knn_result$distance
    # We sum the distances to the nearest neighbor for each point in the cluster
    dist_val[count] <- sum(knn_result$distance[, 1])  # 2nd column for the nearest neighbor (not itself)
    count <- count + 1
    return(sum(dist_val) / length(label))
  }
}
  
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom spdep knn2nb knearneigh nbdists
ELSA <- function(location, label, k=10){
  require(elsa)
  spdf <- SpatialPointsDataFrame(location, data=data.frame(label=label))
  k1 <- knn2nb(knearneigh(location, k=k))
  all.linked <- max(unlist(nbdists(k1, location)))
  elc <- elsa(x=spdf, d=dneigh(spdf, d1=0, d2=all.linked, longlat=FALSE), zcol=label)
  df <- as.data.frame(elc@data)
  colnames(df) <- c("ELSA.Ea", "ELSA.Ec","ELSA")
  return(df)
}

.decideBNPARAM <- function(ncells, BNPARAM){
  if(is.null(BNPARAM)){
    if(nrow(location)>500){
      BNPARAM <- BiocNeighbors::AnnoyParam()
    }else{
      BNPARAM <- BiocNeighbors::ExhaustiveParam()
    }
  }
  return(BNPARAM)
}
