#' @title Compute Pair to Pair Distances
#' @description Compute the pairwise distances between points in matrix `X`.
#' @param X Numeric matrix.
#' @param distance String specifying the metric to compute the distances.
#' @return Numeric matrix of pairwise distances with self-distances set to 
#' `Inf`.
#' @keywords internal
.compute_pair_to_pair_dists <- function(X, distance="euclidean") {
  if(distance == "sqeuclidean"){
    dists <- as.matrix(dist(X, method = "euclidean"))
    dists <- dists^2
  }else{
    dists <- as.matrix(dist(X, method = distance))
  }
  dists[dists < 1e-12] <- 1e-12
  diag(dists) <- Inf
  return(dists)
}


#' @title Get Sub matrix
#' @description Extract a sub matrix from a matrix based on optional row and 
#' column indices.
#' @param arr Numeric matrix.
#' @param inds_a Optional integer vector for row indices.
#' @param inds_b Optional integer vector for column indices.
#' @return Numeric matrix representing the sub matrix.
#' @keywords internal
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
#' @description Computes the internal nodes and edges using Minimum Spanning 
#' Tree.
#' @param mutual_reach_dists Numeric matrix representing mutual reachability 
#' distances.
#' @param use_igraph_mst Logical flag to use MST implementation 
#' in igraph. Currently only mst from igraph is implemented.
#' @return A list containing the indices of internal nodes and their edge 
#' weights.
#' @importFrom Matrix Diagonal
#' @keywords internal
.get_internal_objects <- function(mutual_reach_dists, 
                                 use_igraph_mst=TRUE) {
  rownames(mutual_reach_dists) <- NULL
  colnames(mutual_reach_dists) <- NULL
  if (use_igraph_mst) {
    mst_g <- igraph::mst(igraph::graph_from_adjacency_matrix(
      mutual_reach_dists, mode = "undirected", weighted = TRUE), 
      algorithm = "prim")
    mst <- as.matrix(igraph::as_adjacency_matrix(mst_g, attr = "weight", 
                                                 sparse = FALSE))
  }else{
    stop("Currently only igraph's MST is implemented.")
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
#' @keywords internal
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
#' @keywords internal
.compute_mutual_reach_dists <- function(dists, d) {
  core_dists <- .compute_cluster_core_distance(dists = dists, d = d)
  bc_core_dists <- matrix(rep(core_dists, dim(dists)[2]), nrow = dim(dists)[1], 
                          ncol = dim(dists)[2], byrow=TRUE)
  mutual_reach_dists <- pmax(dists, bc_core_dists)
  mutual_reach_dists <- pmax(mutual_reach_dists, t(bc_core_dists))
  return(list(core_dists = core_dists, mutual_reach_dists = mutual_reach_dists))
}

#' @title Density Sparseness of a Cluster
#' @description Computes the density sparseness for a given cluster.
#' @param cls_inds Integer vector of cluster indices.
#' @param dists Numeric matrix of distances.
#' @param d Integer, the dimensionality.
#' @param use_igraph_mst Logical flag to use MST implementation 
#' in igraph. Currently only mst from igraph is implemented.
#' @keywords internal
#' @return A list containing the density sparseness, internal core distances, 
#' and internal node indices.
.fn_density_sparseness <- function(cls_inds, dists, d, 
                                  use_igraph_mst) {
  mutual_reach <- .compute_mutual_reach_dists(dists = dists, d = d)
  internal_objects <- .get_internal_objects(mutual_reach$mutual_reach_dists, 
                                           use_igraph_mst)
  dsc <- max(internal_objects$internal_edge_weights)
  internal_core_dists <- 
    mutual_reach$core_dists[internal_objects$internal_node_inds]
  internal_node_inds <- cls_inds[internal_objects$internal_node_inds]
  return(list(dsc = dsc, internal_core_dists = internal_core_dists, 
              internal_node_inds = internal_node_inds))
}

#' @title Density Separation of a Pair of Clusters
#' @description Computes the density separation between two clusters.
#' @param cls_i Integer, first cluster index.
#' @param cls_j Integer, second cluster index.
#' @param dists Numeric matrix of distances.
#' @param internal_core_dists_i Numeric vector of core distances for cluster i.
#' @param internal_core_dists_j Numeric vector of core distances for cluster j.
#' @keywords internal
#' @return A list containing the cluster indices and their density separation.
.fn_density_separation <- function(cls_i, cls_j, dists, internal_core_dists_i, 
                                   internal_core_dists_j) {
  bc_internal_core_dists_i <- matrix(rep(internal_core_dists_i, dim(dists)[2]),
                                     nrow = dim(dists)[1],ncol = dim(dists)[2], 
                                     byrow=FALSE)
  sep <- pmax(dists, bc_internal_core_dists_i)
  bc_internal_core_dists_j <- matrix(rep(internal_core_dists_j, dim(dists)[1]),
                                     nrow = dim(dists)[1],ncol = dim(dists)[2], 
                                     byrow=TRUE)
  sep <- pmax(sep, bc_internal_core_dists_j)
  dspc_ij <- ifelse(length(sep) > 0, min(sep), Inf)
  return(list(cls_i = cls_i, cls_j = cls_j, dspc_ij = dspc_ij))
}

#' @title Check Duplicated Samples
#' @description Checks for duplicated samples in matrix `X`.
#' @param X Numeric matrix of samples.
#' @param threshold Numeric, the distance threshold to consider samples as 
#' duplicates.
#' @keywords internal
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
#' @param labels Integer vector of cluster IDs.
#' @param noise_id Integer, the ID for noise.
#' @keywords internal
#' @return Integer vector with singleton clusters converted to noise.
.convert_singleton_clusters_to_noise <- function(labels, noise_id) {
  cluster_sizes <- table(labels)
  singleton_clusters <- names(cluster_sizes[cluster_sizes == 1])
  if (length(singleton_clusters) == 0) {
    return(labels)
  }
  return(ifelse(labels %in% singleton_clusters, noise_id, labels))
}

#' @title Compute Cross Distances
#' @description Computes the cross distances between two sets of indices in matrix `X`.
#' @param X Numeric matrix of samples.
#' @param inds_a Integer vector of indices for the first set.
#' @param inds_b Integer vector of indices for the second set.
#' @param distance String specifying the distance metric.
#' @return Numeric matrix of cross distances.
#' @keywords internal
.compute_cross_dists <- function(X, inds_a, inds_b, distance) {
  if (length(inds_a) == 0 || length(inds_b) == 0) {
    return(matrix(numeric(0), nrow = length(inds_a), ncol = length(inds_b)))
  }
  
  if (distance == "euclidean") {
    a <- X[inds_a, , drop = FALSE]
    b <- X[inds_b, , drop = FALSE]
    aa <- rowSums(a^2)
    bb <- rowSums(b^2)
    ab <- tcrossprod(a, b)
    dists <- sqrt(pmax(outer(aa, bb, "+") - 2 * ab, 0))
  } else if (distance == "sqeuclidean") {
    a <- X[inds_a, , drop = FALSE]
    b <- X[inds_b, , drop = FALSE]
    aa <- rowSums(a^2)
    bb <- rowSums(b^2)
    ab <- tcrossprod(a, b)
    dists <- pmax(outer(aa, bb, "+") - 2 * ab, 0)
  } else {
    # General method for other distances
    dists <- as.matrix(dist(rbind(X[inds_a, , drop = FALSE], 
                            X[inds_b, , drop = FALSE]), 
                         method = distance))
  }
  dists[dists < 1e-12] <- 1e-12
  return(dists)
}

#' @title Calculate DBCV Metric 
#' @description Compute the DBCV (Density-Based Clustering Validation) metric.
#' 
#' @param X Numeric matrix of samples.
#' @param labels Integer vector of cluster IDs.
#' @param distance String specifying the distance metric. `"sqeuclidean"`, or 
#' possible `method` in `stats::dist()`. By default `"euclidean"`.
#' @param noise_id Integer, the cluster ID in `y` for noise (default `-1`). 
#' @param check_duplicates Logical flag to check for duplicate samples.
#' @param use_igraph_mst Logical flag to use `igraph`'s MST 
#' implementation. Currently only `mst` from `igraph` is implemented.
#' @param BPPARAM BiocParallel params for multithreading (default none)
#' @param ... Ignored
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
#' @details This implementation will not fully reproduce the results of other 
#' existing implementations (e.g. \url{https://github.com/FelSiq/DBCV}) due to 
#' the different algorithms used for computing the Minimum Spanning Tree.
#' @export
#' @examples
#' data(noisy_moon)
#' data <- noisy_moon
#' dbcv(data[, c("x", "y")], data$kmeans_label)
#' dbcv(data[, c("x", "y")], data$hdbscan_label)
dbcv <- function(X, labels, distance = "euclidean", noise_id = -1, 
                 check_duplicates = FALSE,
                 use_igraph_mst = TRUE, BPPARAM = BiocParallel::SerialParam(), 
                 ...) {
  X <- as.matrix(X)
  labels <- as.integer(labels)
  n <- nrow(X)
  if (n != length(labels)) {
    stop("Mismatch in number of samples and cluster labels.")
  }
  
  labels <- .convert_singleton_clusters_to_noise(labels, noise_id = noise_id)
  non_noise_inds <- labels != noise_id
  X <- X[non_noise_inds, ]
  labels <- labels[non_noise_inds]
  
  if (length(labels) == 0) {
    return(list(vcs = numeric(0), dbcv = 0.0))
  }
  
  if (check_duplicates) {
    .check_duplicated_samples(X)
  }
  
  cluster_ids <- sort(unique(labels))
  n_clusters <- length(cluster_ids)
  cls_inds <- lapply(cluster_ids, function(id) which(labels == id))
  
  # Initialize storage
  dscs <- numeric(n_clusters)
  min_dspcs <- rep(Inf, n_clusters)
  internal_objects_per_cls <- vector("list", n_clusters)
  internal_core_dists_per_cls <- vector("list", n_clusters)
  
  # Process each cluster
  density_sparseness_results <- BiocParallel::bplapply(
    seq_along(cls_inds),
    function(i) {
      inds <- cls_inds[[i]]
      X_cls <- X[inds, , drop = FALSE]
      dists_intra <- .compute_pair_to_pair_dists(X_cls, distance = distance)
      mutual_reach <- .compute_mutual_reach_dists(dists_intra, d = ncol(X))
      internal_objects <- .get_internal_objects(
        mutual_reach$mutual_reach_dists, 
        use_igraph_mst = use_igraph_mst
      )
      dsc <- max(internal_objects$internal_edge_weights)
      internal_node_inds_global <- inds[internal_objects$internal_node_inds]
      internal_core_dists <- mutual_reach$core_dists[internal_objects$internal_node_inds]
      list(
        dsc = dsc,
        internal_node_inds_global = internal_node_inds_global,
        internal_core_dists = internal_core_dists
      )
    },
    BPPARAM = BPPARAM
  )
  
  # Extract results
  for (i in seq_along(density_sparseness_results)) {
    dscs[i] <- density_sparseness_results[[i]]$dsc
    internal_objects_per_cls[[i]] <- density_sparseness_results[[i]]$internal_node_inds_global
    internal_core_dists_per_cls[[i]] <- density_sparseness_results[[i]]$internal_core_dists
  }
  
  # Process cluster pairs
  if (n_clusters > 1) {
    pairs <- utils::combn(n_clusters, 2, simplify = FALSE)
    density_separation_results <- BiocParallel::bplapply(
      pairs,
      function(pair) {
        i <- pair[1]
        j <- pair[2]
        inds_i <- internal_objects_per_cls[[i]]
        inds_j <- internal_objects_per_cls[[j]]
        if (length(inds_i) == 0 || length(inds_j) == 0) {
          return(list(cls_i = i, cls_j = j, dspc_ij = Inf))
        }
        dists_inter <- .compute_cross_dists(
          X, inds_i, inds_j, distance
        )
        core_i <- internal_core_dists_per_cls[[i]]
        core_j <- internal_core_dists_per_cls[[j]]
        mat_core_i <- matrix(core_i, nrow = length(core_i), ncol = length(core_j), byrow = FALSE)
        mat_core_j <- matrix(core_j, nrow = length(core_i), ncol = length(core_j), byrow = TRUE)
        sep <- pmax(dists_inter, mat_core_i, mat_core_j)
        dspc_ij <- min(sep)
        list(cls_i = i, cls_j = j, dspc_ij = dspc_ij)
      },
      BPPARAM = BPPARAM
    )
    
    for (result in density_separation_results) {
      i <- result$cls_i
      j <- result$cls_j
      min_dspcs[i] <- min(min_dspcs[i], result$dspc_ij)
      min_dspcs[j] <- min(min_dspcs[j], result$dspc_ij)
    }
  }
  
  # Compute validity scores
  vcs <- (min_dspcs - dscs) / pmax(min_dspcs, dscs)
  vcs[is.nan(vcs)] <- 0.0
  dbcv <- sum(vcs * table(labels)) / nrow(X)
  return(list(vcs = vcs, dbcv = dbcv))
}