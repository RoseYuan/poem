#' getPartitionMetrics
#' 
#' Computes a selection of external evaluation metrics for partition.
#'
#' @param true A vector containing the labels of the true classes. Must be a 
#'  vector of characters, integers, numerics, or a factor, but not a list.
#' @param pred A vector containing the labels of the predicted clusters. Must 
#'  be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param metrics The metrics to compute. If omitted, main metrics will be 
#'   computed. See details.
#' @param level The level to calculate the metrics. Options include 
#' `"class"` and `"dataset"`.
#' @param ... 
#' @return A data.frame of metrics.
#' @details
#' The allowed values for `metrics` depend on the value of `level`:
#'   - If `level = "class"`, the allowed `metrics` are: `"WC"`,`"WH"`,`"AWC"`,`"AWH"`,`"FM"` (see below for details).
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
#' @export
#' @examples
#' true <- rep(LETTERS[1:3], each=10)
#' pred <- c(rep("A", 8), rep("B", 9), rep("C", 3), rep("D", 10))
#' getPartitionMetrics(true, pred, level="class")
#' getPartitionMetrics(true, pred, level="dataset")
getPartitionMetrics <-function(true, pred, metrics=c("WC","WH","AWC","AWH","FM"), 
                           level="class", ...){
  # Map level to the corresponding function
  level_functions <- list(
    "class" = getPartitionClassMetrics,
    "dataset" = getPartitionGlobalMetrics
  )
  .checkMetricsLevel(metrics, level, level_functions, use_default=FALSE, 
                     use_attribute=TRUE, attr_name="allowed_metrics")
  # Collect all arguments into a list
  args <- list(true = true, pred = pred, metrics = metrics, ...)
  do.call(level_functions[[level]], args)
}