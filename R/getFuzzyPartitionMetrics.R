
#' getFuzzyPartitionMetrics
#' 
#' Computes a selection of external fuzzy clustering evaluation metrics for spatial data.
#'
#' @param true A vector containing the labels of the true classes. Must be a 
#'  vector of characters, integers, numerics, or a factor, but not a list.
#' @param pred A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param location
#' @param fuzzy_true A boolean indicating if the `true` label will be transformed 
#' into fuzzy cluster membership representation.
#' @param fuzzy_pred A boolean indicating if the `pred` label will be transformed 
#' into fuzzy cluster membership representation.
#' @param k the number of neighbors to consider when transforming the location
#' information into fuzzy class memberships.
#' @param alpha The parameter to control to what extend the spot itself 
#' contribute to the class composition calculation. "equal" means it is 
#' weighted the same as other NNs. A numeric value between 0 and 1 means the 
#' weight of the frequency contribution for the spot itself, and the frequency 
#' contribution for its knn is then 1-alpha.

#' @return A list of metric results.
#' @details
#' Additional details...
#' 
#' @importFrom aricode sortPairs AMI
#' @importFrom clevr mutual_info variation_info homogeneity completeness v_measure
#' @importFrom mclustcomp mclustcomp
#' @importFrom FlowSOM FMeasure
#' @export
getFuzzyPartitionMetrics <-function(true, pred, location, fuzzy_true=TRUE, fuzzy_pred=TRUE, k=6, alpha="equal", ...){
  if(fuzzy_true){
    P <- getFuzzyLabel(true, location, k=k, alpha=alpha)
  }else{
    P <- as.matrix(Matrix::sparseMatrix(i=seq_along(true), j=as.integer(true), x = 1L))
  }
  if(fuzzy_pred){
    Q <- getFuzzyLabel(pred, location, k=k, alpha=alpha)
  }else{
    Q <- as.matrix(Matrix::sparseMatrix(i=seq_along(pred), j=as.integer(pred), x = 1L))
  }
  res <- fuzzyPartitionMetrics(P, Q, ...)
  return(res)
}

