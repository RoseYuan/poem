#' @importFrom igraph subgraph vertex_attr
.igraphFunPerClass <- function(g, classAttr="class", FUN, ...){
  stopifnot(is(g,"igraph"))
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

#' @importFrom igraph mean_distance decompose.graph
.adjMeanShortestPath <- function(g, classAttr="class", directed=FALSE){
  stopifnot(is(g,"igraph"))
  gc <- decompose.graph(g)
  msp <- sapply(gc, directed=directed, FUN=mean_distance)
  return(length(gc)+sum(pmax(1L,msp,na.rm=TRUE)))
}

.simpsonIndex <- function(knn, labels){
  .checkInputs(knn, labels)
  k <- ncol(knn$index)
  p <- sapply(unique(labels), FUN=function(x){
    rowSums(knn$nncl==x)/k
  })
  rowSums(p^2)
}

# neighborhood purity (i.e. proportion with same labels as target node)
.nPurity <- function(knn, labels){
  .checkInputs(knn, labels)
  rowSums(knn$nncl==labels)/ncol(knn$index)
}

# neighborhood class over-representation
.nlog2Enrichment <- function(knn, labels, pseudoCount=1){
  .checkInputs(knn, labels)
  k <- ncol(knn$index)
  expected <- k*as.numeric(table(labels))/length(labels)
  log2((pseudoCount+rowSums(knn$nncl==labels))/
         (pseudoCount+expected[as.integer(labels)]))
}
