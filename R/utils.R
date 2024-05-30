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

.checkInputs <- function(knn, labels, ...){
  .isKnn(knn, ...)
  stopifnot(is.character(labels) || is.factor(labels))
  stopifnot(length(labels)==nrow(knn$index))
}


# computes nearest neighbors 
#' @importFrom BiocNeighbors AnnoyParam ExhaustiveParam findKNN
.emb2knn <- function(x, k, BNPARAM=NULL){
  stopifnot(is.matrix(x) && is.numeric(x))
  stopifnot(is.numeric(k) && length(k)==1 && k>0 && (k %/% 1)==k)
  if(is.null(BNPARAM)){
    if(nrow(x)>500){
      BNPARAM <- BiocNeighbors::AnnoyParam()
    }else{
      BNPARAM <- BiocNeighbors::ExhaustiveParam()
    }
  }
  findKNN(x, k=k, BNPARAM=BNPARAM)
}