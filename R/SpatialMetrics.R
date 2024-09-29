#' nnWeightedAccuracy
#' 
#' Computes an accuracy score which weighs elements/spots that are misclassified
#' by the proportion of their (spatial) neighborhood that is not of the 
#' element/spot's predicted class. This reduces the weight of misclassifications
#' happening at the boundary between domains.
#'
#' @param location The spatial coordinates to compute the nearest neighbors.
#' @param true True class labels (vector coercible to factor)
#' @param pred Predicted labels (vector coercible to factor)
#' @param k Number of nearest neighbors
#' @export
#' @return A scalar representing the weighted accuracy.
nnWeightedAccuracy <- function(true, pred, location, k=5, ...){
  pred <- as.factor(pred)
  matching <- matchSets(pred, true, returnIndices=TRUE)
  pred <- matching[as.integer(pred)]
  true <- as.integer(as.factor(true))
  nn <- findSpatialKNN(location, k, ...)
  w <- which(pred!=true)
  nn <- nn[w]
  nn_dis <- relist(true[unlist(nn)]!=rep(pred[w],lengths(nn)),nn)
  1-sum(sapply(nn_dis, mean))/nrow(location)
}

#' Just the non-spatially-weighted counterpart of nnWeightedAccuracy
#' 
#' @param true True class labels (vector coercible to factor)
#' @param pred Predicted labels (vector coercible to factor)
#'
#' @return A scalar representing the weighted accuracy.
.setMatchingAccuracy <- function(true, pred){
  pred <- as.factor(pred)
  matching <- matchSets(pred, true, returnIndices=TRUE)
  pred <- matching[as.integer(pred)]
  true <- as.integer(as.factor(true))
  sum(pred==true,na.rm=TRUE)/length(pred)
}

#' @title Calculate PAS score
#' @description PAS score measures the clustering performance by calculating 
#' the randomness of the spots that located outside of the spatial region where it was clustered to.
#' Lower PAS score indicates better spatial domian clustering performance.
#' @param label Cluster labels.
#' @param k Number of nearest neighbors.
#' @inheritParams getSpatialGlobalInternalMetrics
#' @return A numeric value for PAS score, and a boolean vector about the abnormal spots.
#' @export
#' @examples 
PAS <- function(label, location, k=10, ...){
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

#' @title Calculate CHAOS score
#' @description CHAOS score measures the clustering performance by calculating 
#' the mean length of the graph edges in the 1-nearest neighbor (1NN) graph 
#' for each cluster, averaged across clusters.
#' Lower CHAOS score indicates better spatial domain clustering performance.
#' @param label Cluster labels.
#' @inheritParams findSpatialKNN
#' @return A numeric value for CHAOS score.
#' @examples 
#' @export
CHAOS <- function(label, location, BNPARAM=NULL) {
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
  
#' Calculate ELSA scores
#' @description
#' Calculating the Entropy-based Local indicator of Spatial Association (ELSA) 
#' scores, which consist of Ea, Ec and the overall ELSA.
#' 
#' @inheritParams PAS
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom spdep knn2nb knearneigh nbdists
#' @importFrom elsa elsa dneigh
#' @references Naimi, Babak, et al., 2019; 10.1016/j.spasta.2018.10.001
#' @return A dataframe containing the Ea, Ec and ELSA for all samples in the dataset.
#' @examples
#' # example code
#' 
#' @export
ELSA <- function(label, location, k=10){
  spdf <- sp::SpatialPointsDataFrame(location, data=data.frame(label=label))
  k1 <- spdep::knn2nb(spdep::knearneigh(location, k=k))
  all.linked <- max(unlist(spdep::nbdists(k1, location)))
  elc <- elsa::elsa(x=spdf, d=elsa::dneigh(spdf, d1=0, d2=all.linked, longlat=FALSE), 
              zcol="label")
  df <- as.data.frame(elc@data)
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

#' getAgreement
#' 
#' Per-spot agreement between a clustering and a ground truth
#'
#' @param pred A vector of predicted clusters
#' @param true A vector of true class labels
#' @param usePairs Logical; whether to compute over pairs instead of elements
#'
#' @return A vector of agreement scores
#'
#' @examples
#' # Give the SPE:
#' spe$agreement <- getAgreement(spe$BayesSpace_default, true=spe$ground_truth)
getAgreement <- function(pred, true, usePairs=TRUE){
  co <- table(true, pred)
  # number of spots in the union between any class and any cluster:
  tot <- matrix(rep(rowSums(co),ncol(co)),nrow=nrow(co))+
    matrix(rep(colSums(co),each=nrow(co)),nrow=nrow(co))-co
  if(usePairs){
    pairs <- choose(co,2)
    truePairsPerCluster <- matrix(rep(colSums(pairs),each=nrow(co)),nrow=nrow(co))
    # wrongPairsPerCluster <- matrix(rep(choose(colSums(co),2),each=nrow(co)),nrow=nrow(co))-
    #   truePairsPerCluster
    truePairs <- truePairsPerCluster + #per cluster
      rowSums(pairs) - pairs # per class minus double-counting
    p <- unclass(truePairs/choose(tot,2))
  }else{
    # intersection over union, i.e. proportion of spots in the class-or-cluster
    # that agree:
    p <- unclass(co/tot)
  }
  # assign each spot its score:
  p <- setNames(as.numeric(p), paste(rep(row.names(p),ncol(p)),rep(colnames(p),each=nrow(p))))
  p[paste(true, pred)]
}