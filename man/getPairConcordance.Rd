% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{getPairConcordance}
\alias{getPairConcordance}
\title{Per-element pair concordance score}
\usage{
getPairConcordance(
  true,
  pred,
  usePairs = TRUE,
  useNegatives = FALSE,
  adjust = FALSE
)
}
\arguments{
\item{true}{A vector of true class labels}

\item{pred}{A vector of predicted clusters}

\item{usePairs}{Logical; whether to compute over pairs instead of elements
Recommended and TRUE by default.}

\item{useNegatives}{Logical; whether to include the consistency of negative
pairs in the score (default FALSE).}

\item{adjust}{Logical; whether to adjust for chance. Only implemented for
\code{useNegatives=FALSE} (doesn't make sense on a element-level otherwise).}
}
\value{
A vector of concordance scores
}
\description{
Per-element pair concordance between a clustering and a ground truth. Note
that by default, negative pairs (i.e. that are split in both the predicted
and true groupings) are not counted. To count it (as in the standard Rand
Index), use \code{useNegatives=TRUE}.
}
