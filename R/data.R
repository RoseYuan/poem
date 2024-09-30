#' @title The noisy moon dataset
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


#' @title Toy examples of spatial data
#' @description
#' Toy examples of spatial data. 
#' 
#' @format ## `sp_toys`
#' A data frame with 240 rows and 11 columns, representing a 16 x 15 array of spots:
#' \describe{
#'   \item{x, y}{Coordinates of each spots.}
#'   \item{row, col}{The row and column index of each spots.}
#'   \item{label}{Ground truth labels. Either 1 or 2.}
#'   \item{p1-p6}{Hypothetical predicted spatial clustering labels.}
#' }
"sp_toys"