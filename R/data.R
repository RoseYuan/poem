#' A simple toy dataset consists of two interleaving half circles.
#'
#' @format ## `noisy_moon`
#' A data frame with 100 rows and 5 columns:
#' \describe{
#'   \item{x, y}{Coordinates of each observations.}
#'   \item{label}{Ground truth labels. Either 1 or 2.}
#'   \item{kmeans_label}{Predicted clustering labels using kmeans with 2 centers.}
#'   \item{hdbscan_label}{Predicted clustering labels using hdbscan with `minPts = 5`.}
#' }
"noisy_moon"