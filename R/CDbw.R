#' Calculate CDbw index
#'
#' Computes the CDbw-index (Halkidi and Vazirgiannis 2008; Halkidi, 
#' Vazirgiannis and Hennig, 2015). This function is directly copied from the 
#' `fpc` CRAN package and was written by Christian Hennig. It is included here
#' to reduce the package dependencies (since `fpc` has a few not-so-light 
#' dependencies that aren't required here).
#'
#' @param x Something that can be coerced into a numerical matrix, with elements
#'   as rows.
#' @param labels A vector of integers with length `=nrow(x)` indicating the 
#'   cluster for each observation.
#' @param r Number of cluster border representatives.
#' @param s Vector of shrinking factors.
#' @param clusterstdev Logical. If `TRUE`, the neighborhood radius for 
#'   intra-cluster density is the within-cluster estimated squared distance from
#'   the mean of the cluster; otherwise it is the average of these over all
#'   clusters.
#' @param trace Logical; whether to print processing info.
#'
#' @return A vector with the following values (see refs for details):
#'   \item{cdbw}{value of CDbw index (the higher the better).}
#'   \item{cohesion}{cohesion.}
#'   \item{compactness}{compactness.}
#'   \item{sep}{separation.} 
#' @export
#' @importFrom stats mahalanobis dist
#' @references 
#' Halkidi, M. and Vazirgiannis, M. (2008) A density-based cluster validity 
#' approach using multi-representatives. Pattern Recognition Letters 29, 
#' 773-786.
#' @references Halkidi, M., Vazirgiannis, M. and Hennig, C. (2015) 
#' Method-independent indices for cluster validation. In C. Hennig, M. Meila, 
#' F. Murtagh, R. Rocci (eds.) Handbook of Cluster Analysis, CRC Press/Taylor 
#' \code{&} Francis, Boca Raton.
#'
#' @author Christian Hennig
#' @examples
#' d1 <- mockData()
#' CDbw(d1[,seq_len(2)], d1[,3])
CDbw <- function(x, labels, r=10, s=seq(0.1,0.8,by=0.1),
                 clusterstdev=TRUE, trace=FALSE){
  labels <- as.integer(as.factor(labels))
  p <- ncol(x)
  n <- nrow(x)
  x <- as.matrix(x)
  cn <- max(labels)
  repc <- list()
  repr <- rep(0,cn)
  mrepr <- vrepc <- numeric(0)
  minrepr <- 1
  xcc <- matrix(0,nrow=cn,ncol=p)
  nc <- numeric(0)
  wvar <- numeric(0)
  repx <- list()
  if (trace)
    message("Find representatives")
  for (i in seq_len(cn)){
    nc[i] <- sum(labels==i)
    xcc[i,] <- colMeans(x[labels==i,,drop=FALSE])
    rrx <- findrep(x,xcc[i,],labels,i,r,p,n,nc[i])
    repc[[i]] <- rrx$repc
    repx[[i]] <- rrx$repx
    repr[i] <- rrx$maxr
    wvar[i] <- rrx$wvar
    mrepr[i] <- sum(repr)
    if(i>1)
      minrepr[i] <- mrepr[i-1]+1
    vrepc <- c(vrepc,repc[[i]]) 
  }
  if (trace)
    print(repc)
  stdev <- mean(sqrt(wvar))
# Find the closest representatives for all representatives and the sets rcr
# for all cluster pairs
  dv <- as.matrix(dist(x[vrepc,]))
  dijmin <- list()
  rcr <- list()
  for (i in seq_len(cn)){
    dijmin[[i]] <- list()
    rcr[[i]] <- list()
  }
#  browser()
  for (i in seq_len(cn-1)){
    for (j in seq(from=i+1L,to=cn)){
      dij <- dv[minrepr[i]:mrepr[i],minrepr[j]:mrepr[j],drop=FALSE]
      ii <- ij <- numeric(0)
      dijmin[[i]][[j]] <- dijmin[[j]][[i]] <- numeric(0)
      for (k in seq_len(repr[i])){
        ii[k] <- which.min(dij[k,])
        dijmin[[i]][[j]][k] <- repc[[j]][ii[k]]
      }
      for (k in seq_len(repr[j])){
        ij[k] <- which.min(dij[,k])
        dijmin[[j]][[i]][k] <- repc[[i]][ij[k]]
      }
      rcr[[i]][[j]] <- numeric(0)
      for (k in seq_len(repr[i])){
        if (k==ij[ii[k]])
          rcr[[i]][[j]] <- rbind(rcr[[i]][[j]],
                                 c(dijmin[[i]][[j]][k],dijmin[[j]][[i]][ii[k]]))
      }
    }
  }
  if(trace){
    print("rcr")
    print(rcr)
    print("wvar")
    print(wvar)
    print("stdev")
    print(stdev)
  }
# Find dens(C_i,C_j), dist
  dens <- dk <- matrix(0,ncol=cn,nrow=cn)
  for (i in seq_len(cn-1)){
    for (j in seq(from=i+1L,to=cn)){      
      nrcr <- nrow(rcr[[i]][[j]])
      wsdij <- mean(sqrt(wvar[i]),sqrt(wvar[j])) 
      for (k in seq_len(nrcr)){
        u <- (x[rcr[[i]][[j]][k,1],]+x[rcr[[i]][[j]][k,2],])/2
#        browser()
        ud <- sqrt(mahalanobis(x[labels==i | labels==j,,drop=FALSE],
                               u,diag(p)))
        dkd <-  sqrt(sum((x[rcr[[i]][[j]][k,1],]-x[rcr[[i]][[j]][k,2],])^2))
        dk[i,j] <- dk[i,j]+dkd
        dens[i,j] <- dens[i,j]+dkd*sum(ud<wsdij)/(2*wsdij*(nc[i]+nc[j]))
      }
      dens[j,i] <- dens[i,j] <- dens[i,j]/nrcr
      dk[j,i] <- dk[i,j] <- dk[i,j]/nrcr
    }
  }
  if(trace){
    print("dens")
    print(dens)
  }
# Interdens and Sep
  maxd <- mind <- numeric(0)
  for (i in seq_len(cn)){
    maxd[i] <- max(dens[i,])
    mind[i] <- min(dk[i,-i])
  }
  interdens <- mean(maxd)
  sep <- mean(mind/(1+interdens))
  if (trace)
    message("sep= ",sep," interdens=",interdens," mind=",mind,"\n")
# Intradens and compactness
  ns <- length(s)
  intradens <- numeric(0)
  denscl <- matrix(0,nrow=cn,ncol=ns) 
  for (i in seq_len(ns)){
    for (j in seq_len(cn)){
#     browser()
      xcj <- x[labels==j,,drop=FALSE]
      if (clusterstdev)
        stdevj <- sqrt(wvar[j])
      else
        stdevj <- stdev
#      dxj <- as.matrix(dist(x))[labels==j,labels==j]
      for (k in seq_len(repr[j])){
        srep <- (1-s[i])*xcj[repx[[j]][k],]+s[i]*xcc[j,]
        dsjk <- mahalanobis(xcj,srep,diag(p))
        denscl[j,i] <- denscl[j,i]+sum(dsjk<stdevj)/nc[j]
      }
      denscl[j,i] <- denscl[j,i]/repr[j]
      if (trace)
        message("denscl cluster ",j," s ",i,": ",denscl[j,i],"\n")
    }
    intradens[i] <- sum(denscl[,i])/(cn*stdev)
  }
  compactness <- mean(intradens)
  if (trace){
    #print(intradens)
    message("compactness= ",compactness,"\n")
  }
# Intrachange and Cohesion
  ic <- intradens[2:ns]-intradens[seq_len(ns-1)]
  intrachange <- sum(ic)/(ns-1)
  cohesion <- compactness/(1+intrachange)
  sc <- sep*compactness
  cdbw <- cohesion*sc
  if (trace){
    #print(intrachange)
    message("cohesion= ",cohesion," sc=",sc," cdbw=",cdbw,"\n")
  }
  c(cdbw=cdbw,cohesion=cohesion,compactness=compactness,sep=sep)
}


# Find representative objects and within-cluster variance

findrep <- function(x,xcen,labels,cluster,r,p=ncol(x),n=nrow(x),
                    nc=sum(labels==cluster)){
  repxx <- matrix(0,nrow=r,ncol=p)
  xc <- x[labels==cluster,,drop=FALSE]
  drx <- matrix(0,nrow=r,ncol=nc)
  drx[1,] <- mahalanobis(xc,xcen,diag(p))
  wvar <- sum(drx[1,])/(nc-1)
  if (is.na(wvar)) wvar <- 0
  #   browser()
  repxi <- which.max(drx[1,])
  repxx[1,] <- xc[repxi,]
  maxr <- r
  #   browser()
  if (r>1){
    for (ri in 2:r){
      drx[ri,] <- mahalanobis(xc,repxx[ri-1,],diag(p))
      di <- numeric(0)
      for (i in seq_len(nc))
        di[i] <- min(drx[2:ri,i])
      if (max(di)>0){
        repxi[ri] <- which.max(di)
        repxx[ri,] <- xc[repxi[ri],]
      }
      else{
        maxr <- ri-1
        break
      }
    }
  }
  list(repc=(seq_len(n))[labels==cluster][repxi],repx=repxi,
       maxr=maxr,wvar=wvar)
}
