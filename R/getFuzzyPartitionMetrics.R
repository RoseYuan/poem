
#' getFuzzyPartitionMetrics
#' 
#' Computes a selection of external fuzzy clustering evaluation metrics for spatial data.
#'
#' @param true A vector containing the labels of the true classes. Must be a 
#'  vector of characters, integers, numerics, or a factor, but not a list.
#' @param pred A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param location A numeric data matrix containing location information, where 
#' rows are points and columns are dimensions. 
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
getFuzzyPartitionMetrics <-function(true, pred, location, fuzzy_true=TRUE, 
                                    fuzzy_pred=FALSE, k=6, alpha="equal", ...){
  if(fuzzy_true & fuzzy_pred){
    message("Comparing between a fuzzy truth and a fuzzy prediction...")
    P <- getFuzzyLabel(true, location, k=k, alpha=alpha)
    Q <- getFuzzyLabel(pred, location, k=k, alpha=alpha)
    res <- fuzzyPartitionMetrics(P, Q, ...)
  }
  if(fuzzy_true & (!fuzzy_pred)){
    message("Comparing between a fuzzy truth and a hard prediction...")
    P <- getFuzzyLabel(true, location, k=k, alpha=alpha)
    res <- fuzzyHardMetrics(pred, true, P, ...)
  }
  if((!fuzzy_true) & fuzzy_pred){
    message("Comparing between a hard truth and a fuzzy prediction...")
    Q <- getFuzzyLabel(pred, location, k=k, alpha=alpha)
    res <- fuzzyHardMetrics(true, pred, Q, ...)
    res <- .switchListItem(mylist, "fuzzyWH", "fuzzyWC")
    res <- .switchListItem(mylist, "fuzzyAWH", "fuzzyAWC")
  }
  if((!fuzzy_true) & (!fuzzy_pred)){
    warning("You are comparing between two hard clusterings! Use function
            `getPartitionMetrics()` to access more metrics for this.")
    res <- getPartitionMetrics(true, pred, 
                               metrics=c("RI","WC","WH","ARI","AWC","AWH"), ...)
  }
  return(res)
}

