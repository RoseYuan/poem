emb2knn <- function(x, k, BNPARAM=NULL){
  stopifnot(is.matrix(x) && is.numeric(x))
  stopifnot(is.numeric(k) && k>0 && (x %/% 1)==x)
  if(is.null(BNPARAM)){
    if(nrow(x)>500){
      BNPARAM <- BiocNeighbors::AnnoyParam()
    }else{
      BNPARAM <- BiocNeighbors::ExhaustiveParam()
    }
  }
  findKNN(x, k=k, BNPARAM=BNPARAM)
}