#' mockData
#' 
#' Generates mock multidimensional data of a given number of classes of points,
#' for testing.
#'
#' @param Ns A vector of more than one positive integers specifying the number 
#'   of elements of each class.
#' @param classDiff The distances between the classes. If there are more than 2
#'   classes, this can be a `dist` object or a symmetric matrix of `length(Ns)-1`
#'   columns/rows where the lower triangle indicates the desired distances
#'   between classes.
#' @param Sds The standard deviation. Can either be a fixed value, a value per 
#'   class, or a matrix of values for each class (rows) and dimension (column).
#' @param ndims The number of dimensions to generate (default 2).
#' @param spread The spread of the points. Can either be a fixed value, a value 
#'   per class, or a matrix of values for each class (rows) and dimension (col).
#' @param rndFn The random function, by default `rnorm`, but should also work
#'   for `rlnorm` and similar.
#'
#' @return A data.frame with coordinates and a `class` column.
#' @export
#'
#' @examples
#' d <- mockData()
mockData <- function(Ns=c(25,15), classDiff=2, Sds=1, ndims=2, spread=c(1,2), rndFn=rnorm){
  stopifnot(ndims>1 & ndims<=pmax(2,length(Ns)-1))
  stopifnot(length(Ns)>1 && length(Ns)<26 && all(Ns>1))
  ncomp <- lower.tri(matrix(nrow=length(Ns),ncol=length(Ns)), diag=FALSE)
  stopifnot(length(classDiff) %in% c(sum(ncomp),1L) && 
              length(spread) %in% c(length(Ns),1L,ndims*length(Ns)))
  if(length(spread)==1) spread <- rep(spread,length(Ns))
  if(length(spread)< (ndims*length(Ns)) ){
    spread <- matrix(rep(spread,ndims),ncol=ndims)
  }
  if(length(Sds)==1) Sds <- rep(Sds,length(Ns))
  if(length(Sds)< (ndims*length(Ns)) ){
    Sds <- matrix(rep(Sds,ndims),ncol=ndims)
  }
  if(!is(classDiff,"dist") && !is.matrix(classDiff)){
    ncomp[which(ncomp)] <- classDiff
    classDiff <- ncomp
  }
  classDiff <- as.dist(classDiff)
  if(length(classDiff)<ndims){
    # only two classes
    means <- rbind(c(0,as.numeric(classDiff)),
                   c(0,0))
  }else{
    means <- suppressWarnings(cmdscale(classDiff, k=ndims))
    while(ncol(means)<ndims){
      means <- cbind(means,0)
    } 
  }
  d <- as.data.frame(lapply(seq_len(ndims), FUN=function(x){
    unlist(lapply(seq_along(Ns), FUN=function(i){
      m <- means[i,x]+scale(seq_len(spread[i,x]), scale=FALSE)
      rndFn(Ns[i], m, Sds[i,x])
    }))
  }))
  colnames(d) <- paste0("V",seq_len(ncol(d)))
  d$class <- factor(rep(LETTERS[seq_along(Ns)],Ns))
  d
}
