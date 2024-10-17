#' getFuzzyPartitionMetrics
#' 
#' Computes a selection of external fuzzy clustering evaluation metrics.
#' @param hardTrue An atomic vector coercible to a factor or integer vector 
#' containing the true hard labels.
#' @param fuzzyTrue A object coercible to a numeric matrix with membership 
#' probability of elements (rows) in clusters (columns). 
#' @param fuzzyPred A object coercible to a numeric matrix with membership 
#'   probability of elements (rows) in clusters (columns).
#' @param metrics The metrics to compute. See details.
#' @param level The level to calculate the metrics. Options include `"element"`, 
#' `"class"` and `"dataset"`.
#' @inheritParams fuzzyPartitionMetrics
#' @inheritParams fuzzyHardMetrics
#' @inheritParams fuzzyHardSpotAgreement
#' @inheritParams getAgreement
#' @param ... Optional arguments for [fuzzyPartitionMetrics()]: `tnorm`. Only 
#' useful when `fuzzy_true=TRUE` and `fuzzy_pred=TRUE`.
#' @details
#' The allowed values for `metrics` depend on the value of `level`:
#'   - If `level = "element"`, the allowed `metrics` are: `"spotAgreement"`.
#'   - If `level = "class"`, the allowed `metrics` are: `"fuzzyWH"`, `"fuzzyAWH"`, `"fuzzyWC"`, `"fuzzyAWC"`.
#'   - If `level = "dataset"`, the allowed `metrics` are: `"fuzzyRI"`, `"fuzzyARI"`, `"fuzzyWH"`, `"fuzzyAWH"`, `"fuzzyWC"`, `"fuzzyAWC"`.

#' @return A dataframe of metric results.
#' @export
#' @examples
#'# generate fuzzy partitions:
#' m1 <- matrix(c(0.95, 0.025, 0.025, 
#'                0.98, 0.01, 0.01, 
#'                0.96, 0.02, 0.02, 
#'                0.95, 0.04, 0.01, 
#'                0.95, 0.01, 0.04, 
#'                0.99, 0.005, 0.005, 
#'                0.025, 0.95, 0.025, 
#'                0.97, 0.02, 0.01, 
#'                0.025, 0.025, 0.95), 
#'                ncol = 3, byrow=TRUE)
#' m2 <- matrix(c(0.95, 0.025, 0.025,  
#'                0.98, 0.01, 0.01, 
#'                0.96, 0.02, 0.02, 
#'                0.025, 0.95, 0.025, 
#'                0.02, 0.96, 0.02, 
#'                0.01, 0.98, 0.01, 
#'                0.05, 0.05, 0.95, 
#'                0.02, 0.02, 0.96, 
#'                0.01, 0.01, 0.98), 
#'                ncol = 3, byrow=TRUE)
#' colnames(m1) <- colnames(m2) <- LETTERS[1:3]
#' getFuzzyPartitionMetrics(fuzzyTrue=m1,fuzzyPred=m2, level="class")
#' 
#' # generate a fuzzy truth:
#' fuzzyTrue <- matrix(c(
#'   0.95, 0.025, 0.025, 
#'   0.98, 0.01, 0.01, 
#'   0.96, 0.02, 0.02, 
#'   0.95, 0.04, 0.01, 
#'   0.95, 0.01, 0.04, 
#'   0.99, 0.005, 0.005, 
#'   0.025, 0.95, 0.025, 
#'   0.97, 0.02, 0.01, 
#'   0.025, 0.025, 0.95), 
#'   ncol = 3, byrow=TRUE)
#' # a hard truth:
#' hardTrue <- apply(fuzzyTrue,1,FUN=which.max)
#' # some predicted labels:
#' hardPred <- c(1,1,1,1,1,1,2,2,2)
#' getFuzzyPartitionMetrics(hardPred=hardPred, hardTrue=hardTrue, fuzzyTrue=fuzzyTrue, nperms=3, level="class")
#' 
getFuzzyPartitionMetrics <- function(hardTrue=NULL, fuzzyTrue=NULL, 
                                     hardPred=NULL, fuzzyPred=NULL, 
                                     metrics=c("fuzzyWH", "fuzzyAWH",
                                               "fuzzyWC", "fuzzyAWC"),
                                     level="class",
                                     nperms=NULL, verbose=TRUE, 
                                     returnElementPairAccuracy=FALSE,
                                     BPPARAM=BiocParallel::SerialParam(), 
                                     useNegatives=TRUE, usePairs=NULL, ...){
  
  if(verbose){
    mc <- match.call()
    mc <- mc[intersect(names(mc), c("hardTrue","fuzzyTrue","hardPred","fuzzyPred"))]
    paste0(names(mc),"=",unlist(mc), collapse=", ")
  }
  if(!is.null(fuzzyTrue)){
    fuzzy_true=TRUE
  }else{
    stopifnot(!is.null(hardTrue))
    fuzzy_true=FALSE
  }
  if(!is.null(fuzzyPred)){
    fuzzy_pred=TRUE
  }else{
    stopifnot(!is.null(hardPred))
    fuzzy_pred=FALSE
  }
  level_functions <- list(
    "element" = getFuzzyPartitionElementMetrics,
    "class" = getFuzzyPartitionClassMetrics,
    "dataset" = getFuzzyPartitionGlobalMetrics
  )
  .checkMetricsLevel(metrics, level, level_functions, use_default=TRUE, 
                     use_attribute=FALSE)
  
  # Collect all arguments into a list
  args <- list(
    "element" = list(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue, 
                     hardPred=hardPred, fuzzyPred=fuzzyPred, 
                     fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                     metrics=metrics,
                     useNegatives=useNegatives, verbose=verbose,
                     usePairs=usePairs),
    "class" = list(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue, 
                   hardPred=hardPred, fuzzyPred=fuzzyPred, 
                   fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                   metrics=metrics, nperms=nperms, verbose=verbose, 
                   returnElementPairAccuracy=returnElementPairAccuracy,
                   BPPARAM=BPPARAM, ...),
    "dataset" = list(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue, 
                    hardPred=hardPred, fuzzyPred=fuzzyPred, 
                    fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                    metrics=metrics, nperms=nperms, verbose=verbose, 
                    returnElementPairAccuracy=returnElementPairAccuracy,
                    BPPARAM=BPPARAM, ...)
  )
  do.call(level_functions[[level]], args[[level]])
}



.cal_fuzzyPartitionMetrics <- function(hardTrue=NULL, fuzzyTrue=NULL, 
                                       hardPred=NULL, fuzzyPred=NULL, 
                                       fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                       nperms=NULL, verbose=TRUE, 
                                       returnElementPairAccuracy=FALSE,
                                       BPPARAM=BiocParallel::SerialParam(), ...){
  if(fuzzy_true & fuzzy_pred){
    stopifnot(!(is.null(fuzzyTrue)|is.null(fuzzyPred)))
    message("Comparing between a fuzzy truth and a fuzzy prediction...")
    res <- fuzzyPartitionMetrics(fuzzyTrue, fuzzyPred, nperms=nperms, verbose=verbose, 
                                 returnElementPairAccuracy=returnElementPairAccuracy,
                                 BPPARAM=BPPARAM, ...)
  }else if(fuzzy_true & (!fuzzy_pred)){
    stopifnot(!(is.null(hardTrue)|is.null(fuzzyTrue)|is.null(hardPred)))
    message("Comparing between a fuzzy truth and a hard prediction...")
    res <- fuzzyHardMetrics(hardTrue, fuzzyTrue, hardPred, nperms=nperms, verbose=verbose, 
                            returnElementPairAccuracy=returnElementPairAccuracy,
                            BPPARAM=BPPARAM)
  }else if((!fuzzy_true) & fuzzy_pred){
    stopifnot(!(is.null(hardTrue)|is.null(fuzzyPred)|is.null(hardPred)))
    message("Comparing between a hard truth and a fuzzy prediction...")
    res <- fuzzyHardMetrics(hardPred, fuzzyPred, hardTrue, nperms=nperms, verbose=verbose, 
                            returnElementPairAccuracy=returnElementPairAccuracy,
                            BPPARAM=BPPARAM)
    res <- .switchListItem(res, "fuzzyWH", "fuzzyWC")
    res <- .switchListItem(res, "fuzzyAWH", "fuzzyAWC")
  }else if((!fuzzy_true) & (!fuzzy_pred)){
    stop("You are comparing between two hard clusterings! Use function
              `getPartitionMetrics()` to do this.")
    
  }
  return(res)
}
  
getFuzzyPartitionGlobalMetrics <- function(hardTrue=NULL, fuzzyTrue=NULL, 
                                           hardPred=NULL, fuzzyPred=NULL, 
                                           fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                           metrics=c("fuzzyRI", "fuzzyARI",
                                                     "fuzzyWH", "fuzzyAWH",
                                                     "fuzzyWC", "fuzzyAWC"),
                                           nperms=NULL, verbose=TRUE, 
                                           returnElementPairAccuracy=FALSE,
                                           BPPARAM=BiocParallel::SerialParam(), ...){

  res <- .cal_fuzzyPartitionMetrics(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue, 
                                    hardPred=hardPred, fuzzyPred=fuzzyPred, 
                                    fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                                    nperms=nperms, verbose=verbose, 
                                    returnElementPairAccuracy=returnElementPairAccuracy,
                                    BPPARAM=BPPARAM, ...)
  fuzzyMetrics <- res
  fuzzyRI <- fuzzyMetrics$NDC
  fuzzyARI <- fuzzyMetrics$ACI
  fuzzyWH <- fuzzyMetrics$fuzzyWH$global
  fuzzyAWH <- fuzzyMetrics$fuzzyAWH$global
  fuzzyWC <- fuzzyMetrics$fuzzyWC$global
  fuzzyAWC <- fuzzyMetrics$fuzzyAWC$global
  res <- lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           fuzzyRI = fuzzyRI,
           fuzzyARI = fuzzyARI,
            fuzzyWH =  fuzzyWH,
           fuzzyAWH = fuzzyAWH,
           fuzzyWC = fuzzyWC,
           fuzzyAWC = fuzzyAWC
    )
  })
  res <- unlist(res)
  return(data.frame(t(res)))
}

getFuzzyPartitionClassMetrics <- function(hardTrue=NULL, fuzzyTrue=NULL, 
                                          hardPred=NULL, fuzzyPred=NULL, 
                                          fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                          metrics=c("fuzzyWH", "fuzzyAWH",
                                                    "fuzzyWC", "fuzzyAWC"),
                                          nperms=NULL, verbose=TRUE, 
                                          returnElementPairAccuracy=FALSE,
                                          BPPARAM=BiocParallel::SerialParam(), ...){
  
  res <- .cal_fuzzyPartitionMetrics(hardTrue=hardTrue, fuzzyTrue=fuzzyTrue, 
                                    hardPred=hardPred, fuzzyPred=fuzzyPred, 
                                    fuzzy_true=fuzzy_true, fuzzy_pred=fuzzy_pred,
                                    nperms=nperms, verbose=verbose, 
                                    returnElementPairAccuracy=returnElementPairAccuracy,
                                    BPPARAM=BPPARAM, ...)
  fuzzyMetrics <- res
  fuzzyWH <- fuzzyMetrics$fuzzyWH$perPartition
  fuzzyAWH <- fuzzyMetrics$fuzzyAWH$perPartition
  fuzzyWC <- fuzzyMetrics$fuzzyWC$perPartition
  fuzzyAWC <- fuzzyMetrics$fuzzyAWC$perPartition
  
  res_class <- lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           fuzzyWC = fuzzyWC,
           fuzzyAWC = fuzzyAWC
    )
  })[metrics %in% c("fuzzyWC", "fuzzyAWC")]
  
  res_class <- as.data.frame(res_class)
  res_class$class <-seq_along(fuzzyWC)
  
  res_cluster <- lapply(setNames(metrics, metrics), FUN=function(m){
    switch(m,
           fuzzyWH = fuzzyWH,
           fuzzyAWH = fuzzyAWH
    )
  })[metrics %in% c("fuzzyWH", "fuzzyAWH")]
  res_cluster <- as.data.frame(res_cluster)
  res_cluster$cluster <- seq_along(fuzzyWH)
  
  res <- .rbind_na(res_class, res_cluster)
  rownames(res) <- NULL
  return(res)
}

#' getFuzzyPartitionElementMetrics
#'
#' Computes a selection of external fuzzy clustering evaluation metrics at the element level.
#' @param metrics The metrics to compute. Currently only `"spotAgreement"` is included at the element level.
#' @inheritParams fuzzyHardSpotAgreement
#' @inheritParams getAgreement
#' @inheritParams getFuzzyPartitionMetrics
#' @param fuzzy_true Logical; whether the truth is fuzzy.
#' @param fuzzy_pred Logical; whether the prediction is fuzzy.
#' @param usePairs Logical; whether to compute over pairs instead of elements. 
#' Only useful when `fuzzy_true=TRUE` and `fuzzy_pred=FALSE`.
#'
#' @return A dataframe of metric values.
#' @examples
#'# generate fuzzy partitions:
#' m1 <- matrix(c(0.95, 0.025, 0.025, 
#'                0.98, 0.01, 0.01, 
#'                0.96, 0.02, 0.02, 
#'                0.95, 0.04, 0.01, 
#'                0.95, 0.01, 0.04, 
#'                0.99, 0.005, 0.005, 
#'                0.025, 0.95, 0.025, 
#'                0.97, 0.02, 0.01, 
#'                0.025, 0.025, 0.95), 
#'                ncol = 3, byrow=TRUE)
#' m2 <- matrix(c(0.95, 0.025, 0.025,  
#'                0.98, 0.01, 0.01, 
#'                0.96, 0.02, 0.02, 
#'                0.025, 0.95, 0.025, 
#'                0.02, 0.96, 0.02, 
#'                0.01, 0.98, 0.01, 
#'                0.05, 0.05, 0.95, 
#'                0.02, 0.02, 0.96, 
#'                0.01, 0.01, 0.98), 
#'                ncol = 3, byrow=TRUE)
#' colnames(m1) <- colnames(m2) <- LETTERS[1:3]
#' getFuzzyPartitionElementMetrics(fuzzyTrue=m1,fuzzyPred=m2, fuzzy_true=TRUE, fuzzy_pred=TRUE)
getFuzzyPartitionElementMetrics <- function(hardTrue=NULL, fuzzyTrue=NULL, 
                                            hardPred=NULL, fuzzyPred=NULL, 
                                            fuzzy_true=TRUE, fuzzy_pred=FALSE,
                                            metrics=c("spotAgreement"),
                                            useNegatives=TRUE, verbose=TRUE,
                                            usePairs=TRUE){
  if(fuzzy_true & fuzzy_pred){
    stopifnot(!(is.null(fuzzyTrue)|is.null(fuzzyPred)))
    message("Comparing between a fuzzy truth and a fuzzy prediction...")
    res <- fuzzySpotAgreement(fuzzyTrue, fuzzyPred)
  }else if(fuzzy_true & (!fuzzy_pred)){
    stopifnot(!(is.null(hardTrue)|is.null(fuzzyTrue)|is.null(hardPred)))
    message("Comparing between a fuzzy truth and a hard prediction...")
    res <- fuzzyHardSpotAgreement(hardTrue, fuzzyTrue, hardPred, 
                                  useNegatives=useNegatives, verbose=verbose)
  }else if((!fuzzy_true) & fuzzy_pred){
    stopifnot(!(is.null(hardTrue)|is.null(fuzzyPred)|is.null(hardPred)))
    message("Comparing between a hard truth and a fuzzy prediction...")
    res <- fuzzyHardSpotAgreement(hardPred, fuzzyPred, hardTrue,
                                  useNegatives=useNegatives, verbose=verbose)
  }else if((!fuzzy_true) & (!fuzzy_pred)){
    message("Comparing between a hard truth and a hard prediction...")
    res <- getAgreement(hardTrue, hardPred, usePairs=usePairs, useNegatives=useNegatives)
  }
  res <- data.frame(res)
  colnames(res) <- "spotAgreement"
  return(res)
}
