
#' Compute internal metrics for spatial data
#' 
#' A generic function to compute a selection of internal clustering evaluation
#'  metrics for spatial data. It can be applied to raw components 
#' (`labels`, `location`) or directly to a `SpatialExperiment` object.
#' @param object The main input. Can be a `SpatialExperiment` object or missing
#'   (when using `labels`, and `location` directly).
#' @param labels When `object` is missing: a vector containing the labels of the 
#' predicted clusters. Must be a vector of characters, integers, numerics, or a 
#' factor, but not a list. When `object` is a `SpatialExperiment` object: the 
#' column name in `colData(object)` containing the labels. 
#' @inheritParams getSpatialElementInternalMetrics
#' @param metrics The metrics to compute. See details.
#' @param level The level to calculate the metrics. Options include `"element"`,
#' `"class"` and `"dataset"`.
#' @return A data.frame of metrics.
#' @importFrom SpatialExperiment spatialCoords SpatialExperiment
#' @export
#' @details
#' The allowed values for `metrics` depend on the value of `level`:
#'   - If `level = "element"`, the allowed `metrics` are: `"PAS"`, `"ELSA"`.
#'   - If `level = "class"`, the allowed `metrics` are: `"CHAOS"`, `"PAS"`, 
#'   `"ELSA"`.
#'   - If `level = "dataset"`, the allowed `metrics` are:
#'      - `"PAS"`: Proportion of abnormal spots (PAS score)
#'      - `"ELSA"`: Entropy-based Local indicator of Spatial Association 
#'      (ELSA score)
#'      - `"CHAOS"`: Spatial Chaos Score.
#'      - `"MPC"`: Modified partition coefficient
#'      - `"PC"`: Partition coefficient
#'      - `"PE"`: Partition entropy
#' @examples
#' # Example with individual components
#' data(sp_toys)
#' data <- sp_toys
#' getSpatialInternalMetrics(labels=data$label, location=data[,c("x", "y")], 
#'                           k=6, level="class")
#' 
#' # Example with SpatialExperiment object
#' se_object <- SpatialExperiment::SpatialExperiment(assays=matrix(NA, 
#'                                              ncol = nrow(data[,c("x", "y")]), 
#'                                              nrow = ncol(data[,c("x", "y")])), 
#'                                spatialCoords=as.matrix(data[,c("x", "y")]))
#' colData(se_object) <- cbind(colData(se_object), data.frame(label=data$label))
#' getSpatialInternalMetrics(object=se_object, labels="label", k=6, 
#'                           level="class")
setGeneric("getSpatialInternalMetrics", signature="object",
           def=function(object=NULL, labels, location=NULL, 
           k=6, alpha=0.5, level="class",
           metrics=c("CHAOS", "PAS", "ELSA"), ...) {
  standardGeneric("getSpatialInternalMetrics")
})


setMethod("getSpatialInternalMetrics", signature(object="missing"), 
          function(labels, location, k, level, metrics, ...) {
                    # input validation
                    if (anyNA(labels)) stop("NA are not supported.")
                    if (is.character(labels)) labels <- as.factor(labels)
                    if (!is.atomic(labels) || 
                        (!is.factor(labels) && !is.integer(labels)))
                      stop("labels must be vector or factor but not list.")
  # Map level to the corresponding function
  level_functions <- list(
    "element" = getSpatialElementInternalMetrics,
    "class" = getSpatialClassInternalMetrics,
    "dataset" = getSpatialGlobalInternalMetrics
  )
  .checkMetricsLevel(metrics, level, level_functions, use_default=FALSE, 
                     use_attribute=TRUE, attr_name="allowed_metrics")
  # Collect all arguments into a list
  args <- list(labels = labels, location=location, k=k, metrics = metrics, ...)
  do.call(level_functions[[level]], args)
          })

setMethod("getSpatialInternalMetrics", signature(object="SpatialExperiment"), 
          function(object, labels, k, level, metrics, ...) {
            if (!labels %in% colnames(colData(object))) {
              stop(paste("The column", labels, "is not present."))
            }
            # Extract true, pred, and location from the SpatialExperiment object
            labels <- colData(object)[, labels]
            location <- data.frame(spatialCoords(object))
            # Call the main function
            getSpatialInternalMetrics(labels=labels, 
                                      location=location, 
                                       k=k, alpha=alpha, level=level, 
                                       metrics=metrics, ...)
          })

#' Compute dataset-level internal evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of internal clustering evaluation metrics for spatial 
#' data at the dataset level. MPC, PC and PE are internal metrics for fuzzy 
#' clustering, and their implementations in package `fclust` are used.
#' 
#' @param labels A vector containing the labels to be evaluated.
#' @param location A numerical matrix containing the location information, with
#' rows as samples and columns as location dimensions.
#' @param k The size of the spatial neighborhood to look at for each spot. 
#' This is used for calculating PAS and ELSA scores.
#' @param metrics The metrics to compute. See below for more details.
#' @param ... Optional arguments for [PAS()].
#' @keywords internal
#' @importFrom fclust MPC PC PE
#' 
#' @references Yuan, Zhiyuan, et al., 2024; 10.1038/s41592-024-02215-8
#' @references Naimi, Babak, et al., 2019; 10.1016/j.spasta.2018.10.001
#' @references Wang, et al., 2022; 10.1016/j.ins.2022.11.010
#' 
#' @return A named vector containing metric values. Possible metrics are:
#'   \item{PAS}{Proportion of abnormal spots (PAS score).}
#'   \item{ELSA}{Entropy-based Local indicator of Spatial Association 
#'   (ELSA score).}
#'   \item{CHAOS}{Spatial Chaos Score.}
#'   \item{MPC}{Modified partition coefficient} 
#'   \item{PC}{Partition coefficient} 
#'   \item{PE}{Partition entropy} 
#' 
getSpatialGlobalInternalMetrics <- function(labels, location, k=6,
                                            metrics=c("PAS", "ELSA", "CHAOS"),
                                            ...){
  res <- lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           PAS = PAS(labels, location, k=k, ...)$PAS,
           ELSA = colMeans(ELSA(labels, location, k=k), na.rm = TRUE),
           CHAOS = CHAOS(labels, location, BNPARAM=NULL)$CHAOS,
           MPC = MPC(getFuzzyLabel(labels, location)),
           PC = PC(getFuzzyLabel(labels, location)),
           PE = PE(getFuzzyLabel(labels, location)),
           stop("Unknown metric.")
    )}
    )
  res <- unlist(res)
  res <- data.frame(t(res), row.names = NULL)
  names(res)[names(res) == "ELSA.ELSA"] <- "ELSA"
  return(res)
}
attr(getSpatialGlobalInternalMetrics, "allowed_metrics") <- c("PAS","ELSA",
                                                              "CHAOS","MPC",
                                                              "PC","PE")

#' Compute spot-level internal evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of internal clustering evaluation metrics for spatial 
#' data at each spot level.
#'
#' @inheritParams getSpatialGlobalInternalMetrics
#' @param metrics Possible metrics: "PAS" and "ELSA".
#' @param ... Optional params for [PAS()].
#' @keywords internal
#' @return A dataframe containing the metric values for all samples in the 
#' dataset.
#' If PAS is calculated, the value is a Boolean about the abnormality of a spot.
#' If ELSA is calculated, Ea, Ec and ELSA for all spots will be returned.
getSpatialElementInternalMetrics <- function(labels, location, k=6, 
                                      metrics=c("PAS", "ELSA"), ...){
  res <- as.data.frame(lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           PAS = PAS(labels, location, k=k, ...)$abnormalty,
           ELSA = ELSA(labels, location, k=k),
           stop("Unknown metric.")
           )})
    )
  colnames(res)[colnames(res) == "ELSA.ELSA"] <- "ELSA"
  return(res)
}
attr(getSpatialElementInternalMetrics, "allowed_metrics") <- c("PAS","ELSA")

#' Compute class-level internal evaluation metrics for spatially-resolved data
#' 
#' Computes a selection of internal clustering evaluation metrics for spatial 
#' data for each class.
#'
#' @inheritParams getSpatialElementInternalMetrics
#' @param metrics Possible metrics: "CHAOS", "PAS" and "ELSA".
#' @param ... Optional params for [PAS()].
#' @keywords internal
#' @return A dataframe of metric values.
getSpatialClassInternalMetrics <- function(labels, location, k=6, 
                                      metrics=c("CHAOS", "PAS", "ELSA"), ...){
  res <- data.frame(class=sort(unique(labels)))
  PAS <- .element2class(data.frame(PAS=PAS(labels, location, 
                                           k=k, ...)$abnormalty, class=labels))
  ELSA <- .element2class(data.frame(ELSA=ELSA(labels, location, 
                                              k=k), class=labels))
  CHAOS <- data.frame(CHAOS=CHAOS(labels, location, BNPARAM=NULL)$CHAOS_class, 
                              class=names(CHAOS(labels, location, 
                                                BNPARAM=NULL)$CHAOS_class))
  if("PAS" %in% metrics){res <- merge(res, PAS, by="class")}
  if("ELSA" %in% metrics){res <- merge(res, ELSA, by="class")}
  if("CHAOS" %in% metrics){res <- merge(res, CHAOS, by="class")}
  colnames(res)[colnames(res) == "ELSA.ELSA"] <- "ELSA"
  return(res)
}
attr(getSpatialClassInternalMetrics, "allowed_metrics") <- c("CHAOS","PAS",
                                                             "ELSA")
