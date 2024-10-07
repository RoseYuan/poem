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

#' @title Toy embedding examples
#' @description
#' Toy example 2D embeddings of elements of different classes, with varying mixing and 
#' spread. Graphs 1-3 all have 20 elements of each of 4 classes, but that are mixed in
#' different fashion in the embedding space. Graphs 4-7 all have 100 elements of class1
#' and 60 of class2, and the class1 elements vary in their spread.
#' 
#' @format ## `toyExamples`
#' A data frame.
#' \describe{
#'   \item{graph}{The name of the embedding to which the element belongs.}
#'   \item{x, y}{Coordinates in the 2D embedding.}
#'   \item{class}{The class to which the element belongs.}
#' }
"toyExamples"


#' @title Metrics Information
#' @description A dataframe storing the information of all metrics

#' @format ## `metric_info`
#' A data frame.
"metric_info"