#' Compute partition-based metrics
#' 
#' Computes a selection of external evaluation metrics for partition.
#'
#' @param true A vector containing the labels of the true classes. Must be a 
#'  vector of characters, integers, numerics, or a factor, but not a list.
#' @param pred A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param metrics The metrics to compute. If omitted, main metrics will be 
#'   computed. See details.
#' @param level The level to calculate the metrics. Options include "element",
#' `"class"` and `"dataset"`.
#' @inheritParams getPartitionGlobalMetrics
#' @return A data.frame of metrics.
#' @details
#' The allowed values for `metrics` depend on the value of `level`:
#'   - If `level = "element"`, the allowed `metrics` are: 
#'      - `"SPC"`: Spot-wise Pair Concordance.
#'      - `"ASPC"`: Adjusted Spot-wise Pair Concordance.
#'   - If `level = "class"`, the allowed `metrics` are: `"WC"`,`"WH"`,`"AWC"`,
#'   `"AWH"`,`"FM"` (see below for details).
#'   - If `level = "dataset"`, the allowed `metrics` are:
#'      - `"RI"`: Rand Index 
#'      - `"WC"`: Wallace Completeness
#'      - `"WH"`: Wallace Homogeneity
#'      - `"ARI"`: Adjusted Rand Index
#'      - `"AWC"`: Adjusted Wallace Completeness
#'      - `"AWH"`: Adjusted Wallace Homogeneity
#'      - `"NCR"`: Normalized class size Rand index
#'      - `"MI"`: Mutual Information
#'      - `"AMI"`: Adjusted Mutual Information
#'      - `"VI"`: Variation of Information
#'      - `"EH"`: (Entropy-based) Homogeneity
#'      - `"EC"`: (Entropy-based) Completeness
#'      - `"VM"`: V-measure
#'      - `"FM"`: F-measure/weighted average F1 score
#'      - `"VDM"`: Van Dongen Measure
#'      - `"MHM"`: Meila-Heckerman Measure
#'      - `"MMM"`: Maximum-Match Measure
#'      - `"Mirkin"`: Mirkin Metric
#'      - `"Accuracy"`: Set Matching Accuracy
#' @export
#' @examples
#' true <- rep(LETTERS[seq_len(3)], each=10)
#' pred <- c(rep("A", 8), rep("B", 9), rep("C", 3), rep("D", 10))
#' getPartitionMetrics(true, pred, level="class")
#' getPartitionMetrics(true, pred, level="dataset")
getPartitionMetrics <-function(true, pred, metrics=NULL, level="class", ...){
  level <- match.arg(level, c("dataset","class","element"))
  # Map level to the corresponding function
  level_functions <- list(
    "element" = getPartitionElementMetrics,
    "class" = getPartitionClassMetrics,
    "dataset" = getPartitionGlobalMetrics
  )
  if(is.null(metrics))
    metrics <- switch(level,
                      "dataset"=c("RI","WC","WH","ARI","NCR","AWC","AWH","MI",
                                  "AMI","VI","EH","EC","VM","FM"),
                      "class"=c("WC","WH","AWC","AWH","FM"),
                      "element"=c("SPC"),
                      stop("Unknown `level` specified.")
    )
  .checkMetricsLevel(metrics, level, level_functions, use_default=FALSE, 
                     use_attribute=TRUE, attr_name="allowed_metrics")
  # Collect all arguments into a list
  args <- list(true = true, pred = pred, metrics = metrics, ...)
  do.call(level_functions[[level]], args)
}

#' getPartitionElementMetrics
#'
#' Computes a selection of external evaluation metrics for partition. The 
#' metrics are reported per element.
#' 
#' @inheritParams getPairConcordance
#' @param metrics The metrics to compute.
#' @keywords internal
#' @return A dataframe of metrics.
getPartitionElementMetrics <- function(true, pred, metrics=c("SPC"), 
                                       usePairs=TRUE, useNegatives=TRUE){
  if (anyNA(true) | anyNA(pred))
    stop("NA are not supported.")
  if (is.character(true)) true <- as.factor(true)
  if (is.character(pred)) pred <- as.factor(pred)
  if (!is.atomic(true) || (!is.factor(true) && !is.integer(true)) ||
      !is.atomic(pred) || (!is.factor(pred) && !is.integer(pred)) )
    stop("true and pred must be vectors or factors but not lists.")
  if(length(true) != length(pred)){
    stop("The two input vectors should have the same length.")
  }

    res <- as.data.frame(lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           SPC = getPairConcordance(true, pred, usePairs=usePairs,
                                    useNegatives=useNegatives, adjust=FALSE),
           ASPC = getPairConcordance(true, pred, usePairs=usePairs,
                                     useNegatives=useNegatives, adjust=TRUE),
           stop("Unknown metric.")
           )})
    )
  return(res)
}
attr(getPartitionElementMetrics, "allowed_metrics") <- c("SPC","ASPC")