#' @importFrom igraph subgraph vertex_attr is_directed
.igraphFunPerClass <- function(g, classAttr="class", directed=is_directed(g), 
                               FUN, ...){
  stopifnot(is(g,"igraph"))
  if(is.null(directed)) directed <- FALSE
  if(is.character(FUN)) FUN <- get(FUN)
  if(!is.null(classAttr)){
    va <- vertex_attr(g,classAttr)
    stopifnot(is.factor(va) || is.character(va) || is.integer(va))
    names(cls) <- cls <- unique(va)
    return(sapply(cls, FUN=function(cl){
      FUN(subgraph(g, v=which(va==cl)), ...)
    }))
  }
  FUN(g, ...)
}

#' @importFrom igraph mean_distance decompose.graph is_directed V
.adjMeanShortestPath <- function(g, directed=FALSE, normalize=TRUE){
  stopifnot(is(g,"igraph"))
  if(is.null(directed)) directed <- FALSE
  gc <- decompose.graph(g)
  msp <- sapply(gc, directed=directed, FUN=mean_distance)
  amsp <- length(gc)+sum(pmax(1L,msp,na.rm=TRUE))
  if(normalize) amsp <- amsp/length(V(g))
  return(amsp)
}

.simpsonIndex <- function(knn, labels=NULL, directed=TRUE){
  if(is.null(directed)) directed <- TRUE
  # for undirected, first transform directed nn to graph
  if(!is(knn,"igraph") && !directed) knn <- .nn2graph(knn, labels)
  if(is(knn,"igraph")){
    if(is.null(labels)) labels <- get_vertex_attr(knn, "class")
    knn <- .igraph2nn(knn, labels, directed=directed)
  }
  stopifnot(!is.null(labels))
  if(.isNNlist(knn)){ # list of varying number of neighbor indices
    knn <- relist(labels[unlist(knn)],knn)
    return(sapply(knn, FUN=function(nn){
        sum((sapply(unique(labels), FUN=function(x) sum(nn=x))/length(nn))^2)
      }))
  }
  .checkInputs(knn, labels)
  k <- ncol(knn$index)
  p <- sapply(unique(labels), FUN=function(x){
    rowSums(knn$nncl==x)/k
  })
  rowSums(p^2)
}

# neighborhood purity (i.e. proportion with same labels as target node)
#' @importFrom igraph as_adj_list
.nPurity <- function(knn, labels=NULL, directed=TRUE){
  if(is.null(directed)) directed <- TRUE
  # for undirected, first transform directed nn to graph
  if(!is(knn,"igraph") && !directed) knn <- .nn2graph(knn, labels)
  if(is(knn,"igraph")){
    if(is.null(labels)) labels <- get_vertex_attr(knn, "class")
    knn <- .igraph2nn(knn, labels, directed=directed)
  }
  stopifnot(!is.null(labels))
  if(.isNNlist(knn)){ # list of varying number of neighbor indices
    return(mapply(nn=relist(labels[unlist(knn)],knn),
                  label=labels, FUN=function(nn,label){
      sum(nn==label)/length(nn)
    }))
  }
  .checkInputs(knn, labels)
  rowSums(knn$nncl==labels)/ncol(knn$index)
}

# neighborhood class over-representation
.nlog2Enrichment <- function(knn, labels, directed=TRUE, pseudoCount=1){
  if(is.null(directed)) directed <- TRUE
  # for undirected, first transform directed nn to graph
  if(!is(knn,"igraph") && !directed) knn <- .nn2graph(knn, labels)
  if(is(knn,"igraph")){
    if(is.null(labels)) labels <- get_vertex_attr(knn, "class")
    knn <- .igraph2nn(knn, labels, directed=directed)
  }
  stopifnot(!is.null(labels))
  if(.isNNlist(knn)){ # list of varying number of neighbor indices
    expect <- table(labels)/length(labels)
    return(mapply(nn=relist(labels[unlist(knn)],knn),
                  label=labels, FUN=function(nn,label){
                    k <- ncol(knn$index)
                    expected <- length(nn)*expect
                    log2((pseudoCount+sum(nn==label))/
                           (pseudoCount+length(nn)*expect))
                  }))
  }
  .checkInputs(knn, labels)
  k <- ncol(knn$index)
  expected <- k*as.numeric(table(labels))/length(labels)
  log2((pseudoCount+rowSums(knn$nncl==labels))/
         (pseudoCount+expected[as.integer(labels)]))
}

#' F measure
#' 
#' Compute the F measure between two clustering results. This is directly copied
#' from the package `FlowSOM`.
#'
#' @param true Array containing real cluster labels for each sample
#' @param pred Array containing predicted cluster labels for each
#'                          sample
#' @param silent    Logical, if FALSE, print some information about 
#'                  precision and recall
#' 
#' @return  F measure score
#' @examples
#' # Generate some random data as an example
#' true <- sample(1:5,100,replace = TRUE)
#' pred <- sample(1:6, 100, replace = TRUE)
#' 
#' # Calculate the FMeasure
#' FMeasure(true,pred)
#' @export
.FMeasure <- function(true, pred, silent=TRUE){
  if (sum(pred)==0)
    return(0);
  a <- table(true, pred);
  p <- t(apply(a,1,function(x)x/colSums(a)))
  r <- apply(a,2,function(x)x/rowSums(a))
  if(!silent) message("Precision: ",
                      sum(apply(p,1,max) * (rowSums(a)/sum(a))),
                      "\nRecall: ",sum(apply(r,1,max) * (rowSums(a)/sum(a))),"\n")
  f <- 2*r*p / (r+p)
  f[is.na(f)] <- 0
  sum(apply(f,1,max) * (rowSums(a)/sum(a)))
}