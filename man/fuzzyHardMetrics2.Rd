% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuzzyPartitionMetrics.R
\name{fuzzyHardMetrics2}
\alias{fuzzyHardMetrics2}
\title{Compute fuzzy-hard metrics with lower memory requirement}
\usage{
fuzzyHardMetrics2(
  hardTrue,
  fuzzyTrue,
  hardPred,
  nperms = 10,
  returnElementPairAccuracy = FALSE,
  verbose = TRUE,
  BPPARAM = BiocParallel::SerialParam()
)
}
\arguments{
\item{hardTrue}{An atomic vector coercible to a factor or integer vector
containing the true hard labels. Must have the same length as \code{hardPred}.}

\item{fuzzyTrue}{A object coercible to a numeric matrix with membership
probability of elements (rows) in clusters (columns). Must have the same
number of rows as the length of \code{hardTrue}. Also note that the columns of
\code{fuzzyTrue} should be in the order of the levels (or integer values) of
\code{hardTrue}.}

\item{hardPred}{An atomic vector coercible to a factor or integer vector
containing the predicted hard labels.}

\item{nperms}{The number of permutations (for correction for chance). If
NULL (default), a first set of 10 permutations will be run to estimate
whether the variation across permutations is above 0.0025, in which case
more (max 1000) permutations will be run.}

\item{returnElementPairAccuracy}{Logical. If TRUE, returns the per-element
pair accuracy instead of the various parition-level and dataset-level
metrics. Default FALSE.}

\item{verbose}{Logical; whether to print info and warnings, including the
standard error of the mean across permutations (giving an idea of the
precision of the adjusted metrics).}

\item{BPPARAM}{BiocParallel params for multithreading (default none)}
}
\value{
A list of metrics:
\item{NDC}{Hullermeier's NDC (fuzzy rand index)}
\item{ACI}{Ambrosio's Adjusted Concordance Index (ACI), i.e. a
permutation-based fuzzy version of the adjusted Rand index.}
\item{fuzzyWH}{Fuzzy Wallace Homogeneity index}
\item{fuzzyWC}{Fuzzy Wallace Completeness index}
\item{fuzzyAWH}{Adjusted fuzzy Wallace Homogeneity index}
\item{fuzzyAWC}{Adjusted fuzzy Wallace Completeness index}
}
\description{
This is a slightly slower, but low-memory version of
\code{\link{fuzzyHardMetrics}}.
}
\examples{
# generate a fuzzy truth:
fuzzyTrue <- matrix(c(
0.95, 0.025, 0.025,
0.98, 0.01, 0.01,
0.96, 0.02, 0.02,
0.95, 0.04, 0.01,
0.95, 0.01, 0.04,
0.99, 0.005, 0.005,
0.025, 0.95, 0.025,
0.97, 0.02, 0.01,
0.025, 0.025, 0.95),
ncol = 3, byrow=TRUE)
# a hard truth:
hardTrue <- apply(fuzzyTrue,1,FUN=which.max)
# some predicted labels:
hardPred <- c(1,1,1,1,1,1,2,2,2)
poem:::fuzzyHardMetrics2(hardTrue, fuzzyTrue, hardPred, nperms=3)
}
\references{
Hullermeier et al. 2012; 10.1109/TFUZZ.2011.2179303;

D'Ambrosio et al. 2021; 10.1007/s00357-020-09367-0
}
\seealso{
\code{\link[=fuzzyHardMetrics]{fuzzyHardMetrics()}}
}
\author{
Pierre-Luc Germain
}
\keyword{internal}
