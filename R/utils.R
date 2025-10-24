# check whether a given input object is a kNN as produced by BiocNeighbors
.isKnn <- function(x, checkNNcl=TRUE, triggerError=TRUE){
  res <- is.list(x) && all(c("index","distance") %in% names(x)) && 
    all(vapply(x[c("index","distance")],FUN.VALUE=logical(1),FUN=is.matrix)) &&
    all(dim(x$index)==dim(x$distance)) && max(x$index) <= nrow(x$index)
  if(res && checkNNcl){
    if(!is.matrix(x$nncl) || !all(dim(x$nncl)==dim(x$index))) res <- FALSE 
  }
  if(res) return(TRUE)
  if(triggerError) stop("The object is not a set of kNN or is not in a ",
                        "BiocNeighbors-like format.")
  FALSE
}

.isNNlist <- function(x){
  !is.matrix(x[[1]]) && is.integer(x[[1]])
}

.checkInputs <- function(knn, labels, ...){
  .isKnn(knn, ...)
  stopifnot(is.character(labels) || is.factor(labels) || is.integer(labels))
  stopifnot(length(labels)==nrow(knn$index))
}

#' Computes k nearest neighbors from embedding
#' 
#' Computes k nearest neighbors from embedding.
#' 
#' @param x A numeric matrix (with features as columns and items as 
#'   rows) from which nearest neighbors will be computed.
#' @param k The number of nearest neighbors.
#' @param BNPARAM A BiocNeighbors parameter object to compute kNNs. Ignored 
#'   unless the input is a matrix or data.frame. If omitted, the Annoy 
#'   approximation will be used if there are more than 500 elements.
#' @return A knn list.
#' @importFrom BiocNeighbors AnnoyParam ExhaustiveParam findKNN
#' @export
#' @examples
#' d1 <- mockData()
#' emb2knn(as.matrix(d1[,seq_len(2)]),k=5)
emb2knn <- function(x, k, BNPARAM=NULL){
  stopifnot(is.matrix(x) && is.numeric(x))
  stopifnot(is.numeric(k) && length(k)==1 && k>0 && (k %/% 1)==k)
  if(is.null(BNPARAM)){
    if(nrow(x)>500){
      BNPARAM <- BiocNeighbors::AnnoyParam()
    }else{
      BNPARAM <- BiocNeighbors::ExhaustiveParam()
    }
  }
  findKNN(as.matrix(x), k=k, BNPARAM=BNPARAM)
}

#' Computes shared nearest neighbors from embedding
#' 
#' computes shared nearest neighbors from embedding.
#' 
#' @param x A numeric matrix (with features as columns and items as 
#'   rows) from which nearest neighbors will be computed.
#' @param k The number of nearest neighbors.
#' @param type A string specifying the type of weighting scheme to use for 
#' shared neighbors.
#' Possible choices include "rank", "number", and "jaccard". See `type` in 
#' [bluster::neighborsToSNNGraph()] for details.
#' @param BNPARAM A BiocNeighbors parameter object to compute kNNs. Ignored 
#'   unless the input is a matrix or data.frame. If omitted, the Annoy 
#'   approximation will be used if there are more than 500 elements.
#' @return An igraph object.
#' @importFrom BiocNeighbors AnnoyParam ExhaustiveParam findKNN
#' @export
#' @examples
#' d1 <- mockData()
#' emb2snn(as.matrix(d1[,seq_len(2)]),k=5)
#' @importFrom bluster neighborsToSNNGraph
emb2snn <- function(x, k, type="rank", BNPARAM=NULL){
  knn <- emb2knn(x, k, BNPARAM=BNPARAM)
  bluster::neighborsToSNNGraph(knn$index, type = type)
}

# computes nearest neighbors from pairwise distance matrix
#' @importFrom methods is
.dist2knn <- function(x, k){
  stopifnot(is(x,"dist"))
  x <- as.matrix(x)
  n <- dim(x)[1]
  knn_distance <- matrix(NA, n, k)
  knn_index <- matrix(NA, n, k)
  
  for(i in seq_len(n)){
    distances <- x[i,]    
    # Exclude the distance to itself
    distances[i] <- Inf
    # Get the indices of the k-nearest neighbors
    neighbors <- order(distances)[seq_len(k)]
    # Store the indices and distances of the k-nearest neighbors
    knn_index[i, ] <- neighbors
    knn_distance[i, ] <- distances[neighbors]
  }
  return(list(index = knn_index, distance = knn_distance))
}

# computes shared nearest neighbors from pairwise distance
#' @importFrom bluster neighborsToSNNGraph
.dist2snn <- function(x, k, type="rank"){
  knn <- .dist2knn(x, k)
  bluster::neighborsToSNNGraph(knn$index, type = type)
}

#' @importFrom igraph as_adj_list
.igraph2nn <- function(x, labels=NULL, directed=TRUE){
  if(is.null(labels)) labels <- vertex_attr(x,"class")
  nn <- lapply(as_adj_list(x, loops="ignore", 
                           mode=ifelse(directed,"out","total")),
               as.integer)
  if(!directed || !all(lengths(nn)==length(nn[[1]]))) return(nn)
  knn <- matrix(unlist(nn),nrow=length(nn),byrow=TRUE)
  knn <- list(index=knn,
              distance=matrix(NA_integer_,nrow=nrow(knn),ncol=ncol(knn)))
  if(!is.null(labels)){
    knn$nncl <- matrix(labels[as.integer(knn$index)], nrow=length(nn))
  }
  knn
}

#' @importFrom bluster neighborsToKNNGraph
#' @importFrom igraph set_vertex_attr
.nn2graph <- function(x, labels=NULL){
  g <- bluster::neighborsToKNNGraph(x$index, directed=TRUE)
  if(!is.null(labels)) g <- set_vertex_attr(g, "class", value=labels)
  g
}

# Adapted from https://github.com/cran/mclustcomp/blob/master/R/auxiliary.R
.aux.conversion <- function(x){
  if (is.character(x)){
    x <- as.numeric(as.factor(unlist(strsplit(x,split=""))))
  } else if (is.factor(x)){
    x <- as.numeric(x)
  } else {
    x <- as.numeric(as.factor(x))
  }
  return(round(x))
}

# switch the values between two named items in a list
.switchListItem <- function(my_list, name1, name2){
  # Switching the values
  temp <- my_list[[name1]]
  my_list[[name1]] <- my_list[[name2]]
  my_list[[name2]] <- temp
  return(my_list)
}

# Check for unrecognized arguments and filter arguments for each function
# example usage: function1(!!!.checkEllipsisArgs(fnList = list(function1, 
# function2), a = 1, b = 2, c = 3)[[1]])
# fnList: A list of functions (defaults to empty).
# ...: Arbitrary named arguments to validate/split.
.checkEllipsisArgs <- function(fnList=list(), ...) {
  args <- list(...)
  formal_args <- lapply(fnList, FUN=\(x) names(formals(x)))
  
  if(length(unknown <- setdiff(names(args), unlist(formal_args)))>0)
    stop("Some unrecognized arguments were given: ", paste(unknown, 
                                                           collapse=", "))
  
  lapply(formal_args, FUN=\(x){
    args[names(args) %in% x]
  })
}

# Check the argument where multiple values can be inputed once
# difference with match.arg: if the input arguments contain anything that is not
# recognised, this will throw an error or warning (depend on the argument 
# "warning").
.checkInvalidArgs <- function(args, allowed_args, arg_name, warning=TRUE){
  valid_args <- match.arg(args, allowed_args, several.ok = TRUE)
  invalid_args <- setdiff(args, valid_args)
  if (length(invalid_args) > 0) {
    if (warning){
      warning("Invalid ", arg_name, ": ", paste(invalid_args, collapse = ", "),
           ". Allowed ", arg_name, " are: ", paste(allowed_args, 
                                                   collapse = ", "))
    }else{
      stop("Invalid ", arg_name, ": ", paste(invalid_args, collapse = ", "),
           ". Allowed ", arg_name, " are: ", paste(allowed_args, 
                                                   collapse = ", "))
    }
  }
}

# for a function and an argument, get all the possibilities for this argument,
# either by looking at the default, or looking at the attribute of the function.
.get_allowed_args <- function(function_with_args, arg_name, use_default=TRUE, 
                                 use_attribute=FALSE, attr_name=NULL){
  if(use_default & (!use_attribute)){
    # Get the allowed args from the default args settings
    allowed_args <- eval(formals(function_with_args)[[arg_name]])
  }else if((!use_default) & use_attribute){
    # Get the allowed args from the function's attributes
    allowed_args <- attr(function_with_args, attr_name)
  }else{stop("Either use default, or use the function attribute.")}
  if (is.null(allowed_args)) {
    stop("The candidate args are not defined for the function.")
  }
  return(allowed_args)
}

# check if the "metrics" argument is valide for the specified "level" of 
# calculation
.checkMetricsLevel <- function(metrics, level, level_functions, ...) {
  # Check if level is valid by looking it up in level_functions
  if (!level %in% names(level_functions)) {
    stop(paste("Invalid level:", level, ". Allowed levels are:", 
               paste(names(level_functions), collapse = ", ")))
  }
  
  # Retrieve the formal argument list of the function for the given level
  function_for_level <- level_functions[[level]]
  # Extract allowed metrics from function
  allowed_metrics <- .get_allowed_args(function_for_level, "metrics",...)  
  # Check that all provided metrics are valid for this level
  .checkInvalidArgs(metrics, allowed_metrics, "metrics", warning=FALSE)
}


.class2global <- function(class_res, summarize_fun=base::mean){
  stopifnot(is.data.frame(class_res))
  class_res$class <- 1
  # use NA-insensitive aggregation by default
  res <- aggregate(. ~ class, data = class_res, 
                   FUN = function(x) summarize_fun(x, na.rm = TRUE))
  subset(res, select = -class)
}

.element2class <- function(element_res, summarize_fun=base::mean){
  stopifnot(is.data.frame(element_res))
  aggregate(. ~ class, data = element_res, FUN = summarize_fun)
}
# row-bind two dataframes with no common columns
.rbind_na <- function(df1, df2){
  if(nrow(df1)==0){return(df2)}
  if(nrow(df2)==0){return(df1)}
  # Find the columns that are missing in each data frame
  missing_cols_df1 <- setdiff(names(df2), names(df1))
  missing_cols_df2 <- setdiff(names(df1), names(df2))
  
  # Add the missing columns with NA values to each data frame
  df1[missing_cols_df1] <- NA
  df2[missing_cols_df2] <- NA
  
  # Now bind the two data frames using rbind
  return(rbind(df1, df2))
}

.decideBNPARAM <- function(ncells, BNPARAM=NULL){
  if(is.null(BNPARAM)){
    if(ncells>500){
      BNPARAM <- BiocNeighbors::AnnoyParam()
    }else{
      BNPARAM <- BiocNeighbors::ExhaustiveParam()
    }
  }
  return(BNPARAM)
}
