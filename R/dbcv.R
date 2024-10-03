#' @title Compute Pair to Pair Distances
#' @description Compute the pairwise distances between points in matrix `X`.
#' @param X Numeric matrix.
#' @param metric String specifying the metric to compute the distances.
#' @return Numeric matrix of pairwise distances with self-distances set to `Inf`.
.compute_pair_to_pair_dists <- function(X, metric="euclidean") {
  if(metric == "sqeuclidean"){
    dists <- as.matrix(dist(X, method = "euclidean"))
    dists <- dists^2
  }else{
    dists <- as.matrix(dist(X, method = metric))
  }
  dists[dists < 1e-12] <- 1e-12
  diag(dists) <- Inf
  return(dists)
}


#' @title Get Sub matrix
#' @description Extract a sub matrix from a matrix based on optional row and column indices.
#' @param arr Numeric matrix.
#' @param inds_a Optional integer vector for row indices.
#' @param inds_b Optional integer vector for column indices.
#' @return Numeric matrix representing the sub matrix.
.get_submatrix <- function(arr, inds_a = NULL, inds_b = NULL) {
  if (is.null(inds_a)) {
    return(arr)
  }
  if (is.null(inds_b)) {
    inds_b <- inds_a
  }
  return(arr[inds_a, inds_b, drop = FALSE])
}

#' @title Get Internal Objects
#' @description Computes the internal nodes and edges using Minimum Spanning Tree.
#' @param mutual_reach_dists Numeric matrix representing mutual reachability distances.
#' @param use_scipy_mst Logical flag to use MST implementation 
#' in scipy. If TRUE, python is required.
#' @return A list containing the indices of internal nodes and their edge weights.
.get_internal_objects <- function(mutual_reach_dists, 
                                 use_scipy_mst=TRUE) {
  rownames(mutual_reach_dists) <- NULL
  colnames(mutual_reach_dists) <- NULL
  if (use_scipy_mst) {
    # Convert the R matrix to a Python object
    mutual_reach_dists_py <- reticulate::r_to_py(mutual_reach_dists)
    # Run the Python code
    scipy <- reticulate::import("scipy.sparse.csgraph")
    np <- reticulate::import("numpy")
    mst <- scipy$minimum_spanning_tree(mutual_reach_dists_py)
    mst <- as.matrix(mst %*% Diagonal(n = ncol(mst)))
    mst <- mst + t(mst)
  }else{
    mst_g <- igraph::mst(igraph::graph_from_adjacency_matrix(
      mutual_reach_dists, mode = "undirected", weighted = TRUE), 
      algorithm = "prim")
    mst <- as.matrix(igraph::as_adjacency_matrix(mst_g, attr = "weight", 
                                                 sparse = FALSE))
  }
  is_mst_edges <- (mst > 0)
  internal_node_inds <- which(colSums(is_mst_edges) > 1)
  internal_edge_weights <- .get_submatrix(mst, inds_a = internal_node_inds)
  return(list(internal_node_inds = internal_node_inds, 
              internal_edge_weights = internal_edge_weights))
}

#' @title Compute Cluster Core Distance
#' @description Computes the core distance for each point in a cluster.
#' @param dists Numeric matrix of distances.
#' @param d Integer, the dimensionality.
#' @return Numeric vector of core distances for each point.
.compute_cluster_core_distance <- function(dists, d) {
  n <- nrow(dists)
  core_dists <- rowSums(dists^(-d)) / (n - 1)
  core_dists <- core_dists^(-1/d)
  return(as.numeric(core_dists))
}

#' @title Compute Mutual Reachability Distances
#' @description Computes the mutual reachability distances between points.
#' @param dists Numeric matrix of distances.
#' @param d Float, the dimensionality.
#' @return A list containing core distances and mutual reachability distances.
.compute_mutual_reach_dists <- function(dists, d) {
  core_dists <- .compute_cluster_core_distance(dists = dists, d = d)
  bc_core_dists <- matrix(rep(core_dists, dim(dists)[2]), nrow = dim(dists)[1], ncol = dim(dists)[2], byrow=TRUE)
  mutual_reach_dists <- pmax(dists, bc_core_dists)
  mutual_reach_dists <- pmax(mutual_reach_dists, t(bc_core_dists))
  return(list(core_dists = core_dists, mutual_reach_dists = mutual_reach_dists))
}

#' @title Density Sparseness of a Cluster
#' @description Computes the density sparseness for a given cluster.
#' @param cls_inds Integer vector of cluster indices.
#' @param dists Numeric matrix of distances.
#' @param d Integer, the dimensionality.
#' @param use_scipy_mst Logical flag to use MST implementation 
#' in scipy. If TRUE, python is required.
#' @return A list containing the density sparseness, internal core distances, and internal node indices.
.fn_density_sparseness <- function(cls_inds, dists, d, 
                                  use_scipy_mst) {
  mutual_reach <- .compute_mutual_reach_dists(dists = dists, d = d)
  internal_objects <- .get_internal_objects(mutual_reach$mutual_reach_dists, 
                                           use_scipy_mst)
  dsc <- max(internal_objects$internal_edge_weights)
  internal_core_dists <- mutual_reach$core_dists[internal_objects$internal_node_inds]
  internal_node_inds <- cls_inds[internal_objects$internal_node_inds]
  return(list(dsc = dsc, internal_core_dists = internal_core_dists, internal_node_inds = internal_node_inds))
}

#' @title Density Separation of a Pair of Clusters
#' @description Computes the density separation between two clusters.
#' @param cls_i Integer, first cluster index.
#' @param cls_j Integer, second cluster index.
#' @param dists Numeric matrix of distances.
#' @param internal_core_dists_i Numeric vector of core distances for cluster i.
#' @param internal_core_dists_j Numeric vector of core distances for cluster j.
#' @return A list containing the cluster indices and their density separation.
.fn_density_separation <- function(cls_i, cls_j, dists, internal_core_dists_i, internal_core_dists_j) {
  bc_internal_core_dists_i <- matrix(rep(internal_core_dists_i, dim(dists)[2]), nrow = dim(dists)[1], ncol = dim(dists)[2], byrow=FALSE)
  sep <- pmax(dists, bc_internal_core_dists_i)
  bc_internal_core_dists_j <- matrix(rep(internal_core_dists_j, dim(dists)[1]), nrow = dim(dists)[1], ncol = dim(dists)[2], byrow=TRUE)
  sep <- pmax(sep, bc_internal_core_dists_j)
  dspc_ij <- ifelse(length(sep) > 0, min(sep), Inf)
  return(list(cls_i = cls_i, cls_j = cls_j, dspc_ij = dspc_ij))
}

#' @title Check Duplicated Samples
#' @description Checks for duplicated samples in matrix `X`.
#' @param X Numeric matrix of samples.
#' @param threshold Numeric, the distance threshold to consider samples as duplicates.
#' @return None
.check_duplicated_samples <- function(X, threshold = 1e-9) {
  if (nrow(X) <= 1) {
    return()
  }
  dists <- as.matrix(dist(X))
  if (any(dists < threshold)) {
    stop("Duplicated samples have been found in X.")
  }
}

#' @title Convert Singleton Clusters to Noise
#' @description Converts clusters containing a single instance to noise.
#' @param y Integer vector of cluster IDs.
#' @param noise_id Integer, the ID for noise.
#' @return Integer vector with singleton clusters converted to noise.
.convert_singleton_clusters_to_noise <- function(y, noise_id) {
  cluster_sizes <- table(y)
  singleton_clusters <- names(cluster_sizes[cluster_sizes == 1])
  if (length(singleton_clusters) == 0) {
    return(y)
  }
  return(ifelse(y %in% singleton_clusters, noise_id, y))
}

#' @title DBCV Metric Calculation
#' @description Compute the DBCV (Density-Based Clustering Validation) metric.
#' 
#' @param X Numeric matrix of samples.
#' @param y Integer vector of cluster IDs.
#' @param metric String specifying the distance metric. `"sqeuclidean"`, or 
#' possible `method` in `stats::dist()`. By default `"euclidean"`.
#' @param noise_id Integer, the cluster ID in `y` for noise (default `-1`). 
#' @param check_duplicates Logical flag to check for duplicate samples.
#' @param use_scipy_mst Logical flag to use `scipy`'s Kruskal's MST 
#' implementation. If `TRUE`, python is required, and this will reproduce the same results as
#' this python implementation of DBCV at \url{https://github.com/FelSiq/DBCV}.
#' If `FALSE`, use MST implementation in `igraph`.
#' @param BPPARAM BiocParallel params for multithreading (default none)
#' 
#' @importFrom BiocParallel SerialParam bplapply
#' 
#' @return A list:
#' \item{vcs}{Numeric vector of validity index for each cluster.}
#' \item{dbcv}{Numeric value representing the overall DBCV metric.}
#'
#' @importFrom utils combn
#' 
#' @references Davoud Moulavi, et al. 2014; 10.1137/1.9781611973440.96.
#' 
#' @examples
#' data <- noisy_moon
#' dbcv(data[, c("x", "y")], data$kmeans_label)
#' dbcv(data[, c("x", "y")], data$hdbscan_label)
dbcv <- function(X, y, metric = "euclidean", noise_id = -1, check_duplicates = FALSE,
                 use_scipy_mst = FALSE, BPPARAM=BiocParallel::SerialParam(), ...) {
  X <- as.matrix(X)
  y <- as.integer(y)
  n <- dim(X)[1]
  if (nrow(X) != length(y)) {
    stop("Mismatch in number of samples and cluster labels.")
  }
  
  y <- .convert_singleton_clusters_to_noise(y, noise_id = noise_id)
  
  non_noise_inds <- y != noise_id
  X <- X[non_noise_inds, ]
  y <- y[non_noise_inds]
  
  if (length(y) == 0) {
    return(0.0)
  }
  
  cluster_ids <- sort(unique(y))
  
  if (check_duplicates) {
    .check_duplicated_samples(X)
  }
  
  dists <- .compute_pair_to_pair_dists(X, metric = metric)
  
  dscs <- numeric(length(cluster_ids))
  min_dspcs <- rep(Inf, length(cluster_ids))
  
  internal_objects_per_cls <- list()
  internal_core_dists_per_cls <- list()
  
  cls_inds <- lapply(cluster_ids, function(cls_id) which(y == cls_id))
  
  if (n_processes == "auto") {
    n_processes <- ifelse(length(y) > 500, 4, 1)
  }
  
  density_sparseness_results <- bplapply(seq_along(cls_inds), function(cls_id) {
    .fn_density_sparseness(cls_inds[[cls_id]], .get_submatrix(dists, inds_a = cls_inds[[cls_id]]), 
                          d = ncol(X), use_scipy_mst = use_scipy_mst)
  }, BPPARAM)
  
  for (cls_id in seq_along(density_sparseness_results)) {
    internal_objects_per_cls[[cls_id]] <- density_sparseness_results[[cls_id]]$internal_node_inds
    internal_core_dists_per_cls[[cls_id]] <- density_sparseness_results[[cls_id]]$internal_core_dists
    dscs[cls_id] <- density_sparseness_results[[cls_id]]$dsc
  }
  
  if (length(cluster_ids) > 1) {
    density_separation_results <- bplapply(combn(length(cluster_ids), 2, simplify = FALSE), function(pair) {
      .fn_density_separation(pair[1], pair[2], .get_submatrix(dists, inds_a = internal_objects_per_cls[[pair[1]]], inds_b = internal_objects_per_cls[[pair[2]]]),
                            internal_core_dists_per_cls[[pair[1]]], internal_core_dists_per_cls[[pair[2]]])
    }, BPPARAM)
    
    for (result in density_separation_results) {
      min_dspcs[result$cls_i] <- min(min_dspcs[result$cls_i], result$dspc_ij)
      min_dspcs[result$cls_j] <- min(min_dspcs[result$cls_j], result$dspc_ij)
    }
  }
  
  vcs <- (min_dspcs - dscs) / pmax(min_dspcs, dscs)
  vcs[is.nan(vcs)] <- 0.0
  dbcv <- sum(vcs * table(y)) / n
  return(list("vcs"=vcs ,"dbcv"=dbcv))
}
