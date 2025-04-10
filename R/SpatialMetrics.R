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
#' @param ... Optional params passed to [findSpatialKNN()].
#' @keywords internal
#' @return A scalar representing the weighted accuracy.
nnWeightedAccuracy <- function(true, pred, location, k=5, ...){
  pred <- as.factor(pred)
  matching <- matchSets(pred, true, returnIndices=TRUE)
  pred <- matching[as.integer(pred)]
  true <- as.integer(as.factor(true))
  nn <- findSpatialKNN(location, k, ...)
  w <- which(pred!=true)
  if(length(which(pred!=true))==0){
    return(0)
  }
  nn <- nn[w]
  nn_dis <- relist(true[unlist(nn)]!=rep(pred[w],lengths(nn)),nn)
  1-sum(vapply(nn_dis, mean, FUN.VALUE=numeric(1L)))/nrow(location)
}

#' @title Calculate PAS score
#' @description PAS score measures the clustering performance by calculating 
#' the randomness of the spots that located outside of the spatial region where it was clustered to.
#' Lower PAS score indicates better spatial domian clustering performance.
#' @param labels Cluster labels.
#' @param k Number of nearest neighbors.
#' @param ... Optional params for [findSpatialKNN()].
#' @inheritParams getSpatialGlobalInternalMetrics
#' @return A numeric value for PAS score, and a boolean vector about the abnormal spots.
#' @export
#' @examples 
#' data(sp_toys)
#' data <- sp_toys
#' PAS(data$label, data[,c("x", "y")], k=6)
#' PAS(data$p1, data[,c("x", "y")], k=6)
#' PAS(data$p2, data[,c("x", "y")], k=6)
PAS <- function(labels, location, k=10, ...){
  stopifnot(!any(is.na(labels)))
  comp <- knnComposition(location, k=k, labels, alpha=0, ...)
  prop <- unlist(lapply(seq_len(dim(comp)[1]), function(i){comp[i,labels[i]]}))
  results <- prop < 0.5
  return(list(PAS=sum(results)/length(results), abnormalty=results))
}

#' @title Calculate CHAOS score
#' @description CHAOS score measures the clustering performance by calculating 
#' the mean length of the graph edges in the 1-nearest neighbor (1NN) graph 
#' for each cluster, averaged across clusters.
#' Lower CHAOS score indicates better spatial domain clustering performance.
#' @param labels Cluster labels.
#' @inheritParams findSpatialKNN
#' @return A numeric value for CHAOS score.
#' @examples 
#' data(sp_toys)
#' data <- sp_toys
#' CHAOS(data$label, data[,c("x", "y")])
#' CHAOS(data$p1, data[,c("x", "y")])
#' CHAOS(data$p2, data[,c("x", "y")])
#' @export
CHAOS <- function(labels, location, BNPARAM=NULL) {
  BNPARAM <- .decideBNPARAM(nrow(location), BNPARAM)
  labels <- as.vector(labels)
  location <- as.matrix(location)
  matched_location <- scale(location)
  label_unique <- unique(labels)
  # Initialize a vector to hold the distance values
  dist_val <- numeric(length(label_unique))
  
  for (count in seq_along(label_unique)) {
    k <- label_unique[count]
    # Get the locations belonging to the current cluster
    location_cluster <- matched_location[labels == k, ]
    # Skip clusters with fewer than or equal to 2 points
    if (nrow(location_cluster) <= 2) {
      next
    }
    # Find the nearest neighbors within the cluster using BiocNeighbors::findKNN
    knn_result <- BiocNeighbors::findKNN(as.matrix(location_cluster), k = 1,
                                         warn.ties=FALSE, get.index = FALSE, BNPARAM=BNPARAM)
    # The distances to the nearest neighbors are stored in the knn_result$distance
    # We sum the distances to the nearest neighbor for each point in the cluster
    dist_val[count] <- sum(knn_result$distance[, 1])  # 2nd column for the nearest neighbor (not itself)
  }
  res_class <- dist_val/unlist(lapply(label_unique, function(x){sum(labels==x)}))
  names(res_class) <- label_unique
  return(list(CHAOS=sum(dist_val) / length(labels), CHAOS_class=res_class))
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
#' data(sp_toys)
#' data <- sp_toys
#' ELSA(data$label, data[,c("x", "y")], k=6)
#' ELSA(data$p1, data[,c("x", "y")], k=6)
#' ELSA(data$p2, data[,c("x", "y")], k=6)
#' @export
ELSA <- function(labels, location, k=10){
  labels <- as.numeric(labels)
  spdf <- sp::SpatialPointsDataFrame(location, data=data.frame(label=labels))
  k1 <- spdep::knn2nb(spdep::knearneigh(location, k=k))
  all.linked <- max(unlist(spdep::nbdists(k1, location)))
  elc <- elsa::elsa(x=spdf, d=elsa::dneigh(spdf, d1=0, d2=all.linked, longlat=FALSE), 
              zcol="label")
  df <- as.data.frame(elc@data)
  return(df)
}

#' Per-element local concordance between a clustering and a ground truth
#' 
#' Per-element local concordance between a clustering and a ground truth
#'
#' @param location A matrix or data.frame with spatial dimensions as columns.
#'   Alternatively, a nearest neighbor object as produced by 
#'   \code{\link[BiocNeighbors]{findKNN}}.
#' @param pred A vector of predicted clusters
#' @param true A vector of true class labels
#' @param k Approximate number of nearest neighbors to consider
#' @param useNegatives Logical; whether to include the concordance of negative
#'   pairs in the score (default FALSE).
#' @param distWeights Logical; whether to weight concordance by distance 
#'   (default TRUE).
#' @param BNPARAM A BiocNeighbors parameter object to compute kNNs. Ignored 
#'   unless the input is a matrix or data.frame. If omitted, the Annoy 
#'   approximation will be used if there are more than 500 elements.
#' @return A vector of concordance scores
#' @examples
#' data(sp_toys)
#' data <- sp_toys
#' getNeighboringPairConcordance(data$label, data$p1, data[,c("x", "y")], k=6)
#' @export
getNeighboringPairConcordance <- function(true, pred, location, k=20L,
                                       useNegatives=FALSE, distWeights=TRUE,
                                       BNPARAM=NULL){
  if(.isKnn(location, checkNNcl=FALSE, triggerError=FALSE)){
    stopifnot(ncol(location$index)>=k)
    nn <- location
  }else{
    BNPARAM <- .decideBNPARAM(nrow(location), BNPARAM)
    nn <- BiocNeighbors::findKNN(as.matrix(location), k=k*2, warn.ties=FALSE,
                                 BNPARAM=BNPARAM)
  }
  mkd <- median(nn$distance[,k])
  vapply(seq_len(nrow(nn[[1]])), FUN.VALUE=numeric(1), FUN=function(i){
    d <- nn$distance[i,]
    w <- which(d<=mkd & d>0)
    d <- d[w]
    true_conc <- true[nn$index[i,w]]==true[i]
    pred_conc <- pred[nn$index[i,w]]==pred[i]
    if(!useNegatives){
      w <- which(true_conc | pred_conc)
      true_conc <- true_conc[w]
      pred_conc <- pred_conc[w]
      d <- d[w]
    }
    if(distWeights){
      w <- (1/d)/sum(1/d)
    }else{
      w <- 1/length(d)
    }
    sum((1-abs(true_conc-pred_conc))*w)
  })
}


#' spatialARI
#' 
#' Computes the spatial Rand Index and spatial ARI (Yan, Feng and Luo, 2025).
#' See details for more information.
#'
#' @param pred A vector of predicted clusters
#' @param true A vector of true class labels
#' @param coords A matrix of spatial coordinates, with dimensions as columns
#' @param normCoords Logical; whether to normalize the coordinates to 0-1.
#' @param alpha The alpha used in the `f` and `h` functions (default 0.8).
#' @param fbeta,hbeta Additional factors used in the exponential decay functions 
#'   (see details). A higher value means a faster decay. These are ignored if 
#'   `original=TRUE`.
#' @param spotWise Logical; whether to return the spot-wise spatial concordance
#'   (not adjusted for chance).
#' @param nChunks The number of processing chunks. If NULL, this will be 
#'   determined automatically based on the size of the dataset.
#' @param original Logical; whether to use the original h/f functions from Yan, 
#'   Feng and Luo (default FALSE).
#' @param f The f function, which determines the positive contribution of pairs
#'   that are in different partitions in the reference, but grouped together in
#'   the clustering, based on the distance between mates.
#' @param h The h function, which determines the positive contribution of pairs
#'   that are in the same partition in the reference, but different ones in
#'   the clustering, based on the distance between mates.
#' @param BPPARAM Optional BiocParallel BPPARAM object for multithreading. This
#'   only makes sense if the dataset is very large and the data is split into
#'   many chunks.
#'
#' @author Pierre-Luc Germain
#' @references
#' Yan, Feng and Luo, biorxiv 2025, https://doi.org/10.1101/2025.03.25.645156
#' 
#' @details
#' This is a reimplementation of the method from the `spARI` package, made more
#' scalable through chunk-based processing, and with some additional options. 
#' Note that by default, this will not produce the same results as the original 
#' method: to do so, set `original=TRUE`. In our exploration of the method and 
#' its behavior, we found the decay to be too slow, and we therefore 1) do not 
#' square the distances, and 2) introduced a beta parameter in each function 
#' which allows to scale it (a higher beta parameter means a faster decay).
#' 
#' @return A vector containing the spatial Rand Index (spRI) and spatial 
#'   adjusted Rand Index (spARI). Alternatively, if `spotWise=TRUE`, a vector 
#'   of spatial pair concordances for each spot.
#' @importFrom matrixStats rowMins rowMaxs
#' @importFrom pdist pdist
#' @export
spatialARI <- function(true, pred, coords, normCoords=TRUE, alpha=0.8, fbeta=4,
                       hbeta=1, spotWise=FALSE,  nChunks=NULL, original=FALSE,
                       f=function(x){ alpha*exp(-x*fbeta) },
                       h=function(x){ alpha*(1-exp(-x*hbeta)) },
                       BPPARAM=SerialParam()){
  if(isTRUE(original)){
    f <- function(x){ alpha*exp(-x^2) }
    h <- function(x){ alpha*(1-exp(-x^2)) }
  }
  stopifnot(is.function(h) && is.function(f))
  N <- length(true)
  stopifnot(length(pred)==N && nrow(coords)==N)
  if(normCoords){
    coords <- t(coords)
    coords <- t((coords-matrixStats::rowMins(coords))/
                  (matrixStats::rowMaxs(coords)-matrixStats::rowMins(coords)))
  }

  n_choose = choose(N, 2)
  p <- sum(choose(table(true), 2))/n_choose
  q <- sum(choose(table(pred), 2))/n_choose
  
  if(is.null(nChunks)) nChunks <- ceiling(length(true)^2/2e+7)
  
  sp <- split(seq_len(N), head(rep(seq_len(nChunks), ceiling(N/nChunks)), N))

  o <- bplapply(sp, BPPARAM=BPPARAM, FUN=function(i){
    if(nChunks==1){
      di <- as.matrix(dist(coords))
      same_in_pred <- outer(pred, pred, `==`)
      same_in_true <- outer(true, true, `==`)
    }else{
      di <- as.matrix(pdist::pdist(coords, indices.A=i, indices.B=seq_len(N)))
      same_in_pred <- outer(pred[i], pred, `==`)
      same_in_true <- outer(true[i], true, `==`)
    }
    hv <- h(di)
    fv <- f(di)
    w <- ifelse(same_in_pred==same_in_true, 1, ifelse(same_in_pred, fv, hv))
    spc <- (rowSums(w)-1)/(ncol(w)-1L)
    if(spotWise) return(spc)
    list(spc, sum(hv), sum(fv))
  })
  
  if(spotWise) return(unlist(o))
  spRI <- mean(unlist(lapply(o, \(x) x[[1]])))
  o <- matrix(unlist(lapply(o, \(x) x[2:3])), ncol=2, byrow = TRUE)
  sH <- sum(o[,1])/2 - h(0)*length(true)
  sF <- sum(o[,2])/2

  expSpRI <- 1 + 2*p*q - (p+q) + (sF*q + sH*p - (sH+sF)*p*q)/n_choose
  spARI <- (spRI-expSpRI)/(1-expSpRI)
  c(spRI=spRI, spARI=spARI)
}
