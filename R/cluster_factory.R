###########################################################################
# CLUSTER FACTORY
# cluster assembly and mapper graph statistics
###########################################################################


#' Subset a distance matrix
#'
#' @param patch A list of names of data points.
#' @param dists A distance matrix for data points in the patch, possibly including extra points.
subset_dists <- function(patch, dists) {
  patch_size = length(patch)
  if (patch_size == 0) {
    return(NA)
  } else if (patch_size == 1) {
    return(patch)
  } else {
    res = as.dist(as.matrix(dists)[patch, patch]) # this is how it's done in the usedist package
    return(res)
  }
}

# goblin clustering mines -------------------------------------------------

#' Ship data off to the clustering goblins
#'
#' This function tells the computer to look away for a second, so the goblins come and cluster your data while it's not watching.
#'
#' @param dist_mats A list of distance matrices of each bin that is to be clustered.
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix in a list.
#'
#' @return The output of `clusterer(dist_mats)`, which needs to be a list containing named vectors (one per bin), whose names are data point names and whose values are cluster labels (within each bin)
get_raw_clusters <- function(dist_mats, clusterer) {
  return(clusterer(dist_mats))
}



#' Perform the clustering step in mapper
#'
#' This function processes the binned data and global distance matrix to return freshly clustered data.
#'
#' @param bins A list containing "bins" of vectors of names of data points.
#' @param dists A distance matrix containing pairwise distances between named data points.
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix.
#'
#' @return The output of `clusterer` applied to a list of distance matrices, which should be a list containing named vectors (one per bin), whose names are data point names and whose values are cluster labels.
get_clusters <- function(bins, dists, clusterer) {
  # more than one bin, need more than one distance matrix
  if (is.list(bins)) {
    # subset the global distance matrix per bin
    dist_mats = mapply(subset_dists, bins, MoreArgs = list(dists = dists), SIMPLIFY = FALSE)

    # cluster the data
    clusters = get_raw_clusters(dist_mats, clusterer)

    # accurately total up clusters
    clusters_per_bin = sapply(clusters, max)
    offset = c(0, cumsum(clusters_per_bin))
    clusters = mapply(function(x, y)
      x + y, clusters, offset[-length(offset)], SIMPLIFY = FALSE)
    return(clusters)
  }
#
  # cluster the data
  clusters = get_raw_clusters(subset_dists(bins, dists), clusterer) # this fixed everything????

  return(clusters)
}

# node data --------------------------------------------------------------


#' Compute cluster sizes
#'
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster ids.
#'
#' @return A vector of integers representing the lengths of the clusters in the input data.
get_cluster_sizes <- function(binclust_data) {

  # no need to list by level set
  flattened_data = unlist(binclust_data)
  num_vertices = max(flattened_data)

  # find data points with a specific cluster id
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])

  # return counts of the different clusters
  return(sapply(clusters, length))
}

#' Recover bins
#'
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster ids.
#'
#' @return A vector of integers equal in length to the number of clusters, whose members identify which bin that cluster belongs to.
get_bin_vector <- function(binclust_data) {

  # find unique clusters in each bin, then count how many
  num_unique_clusters_per_bin = sapply(lapply(binclust_data, unique), length)

  # repeat the bin id for as many clusters belonging to that bin
  bin_by_clusters = unlist(mapply(
    function(x, y)
      rep(x, y),
    1:length(num_unique_clusters_per_bin),
    num_unique_clusters_per_bin, SIMPLIFY = FALSE))

  return(bin_by_clusters)
}

#' Compute dispersion of a single cluster
#'
#' @param dists A distance matrix for points in the cluster.
#' @param cluster A list containing named vectors, whose names are data point names and whose values are cluster labels
#'
#' @return A real number in \eqn{[0,1]} representing the mean distance to the medoid of the cluster.
#' @details This method finds the medoid of the input data set and returns the average distance to the medoid, i.e., \deqn{\tau(C) = \dfrac{1}{\left(|C|-1\right)}\displaystyle\sum_{i}\text{dist}(x_i, x_j)} where \deqn{x_j = \text{arg}\,\min\limits_{x_j\in C}\, \sum_{x_i \in C, i\neq j}\text{dist}(x_i, x_j)} A smaller value indicates a tighter cluster based on this metric.
compute_tightness <- function(dists, cluster) {

  # empty or singleton clusters have trivial tightness
  if (length(cluster) <= 1) {
    return(0)
  } else {
    # get the distances associated to points in this cluster
    cluster_names = names(cluster)
    these_dists = dists[cluster_names, cluster_names]

    # find minimum sum of distances
    sums = apply(these_dists, 1, sum)
    min_sum = min(sums)

    # use the maximum distance to the medoid to calculate tightness
    medoid_dists = sample(which(sums == min_sum), 1) # pick a medoid
    min_dists = these_dists[medoid_dists, ]
    closeness_factor = min_sum / (length(cluster) - 1)

    return(closeness_factor)
  }
}

#' Compute dispersion measures of a list of clusters
#'
#' @param dists A distance matrix for the data points inside all the input clusters
#' @param binclust_data A list of named vectors whose names are those of data points and whose values are cluster ids
#'
#' @return A vector of real numbers in \eqn{(0,\infty)} containing mean distances to the medoids of each cluster in `dists`.
get_cluster_tightness_vector <- function(dists, binclust_data) {

  # no need to list by level set
  flattened_data = unlist(binclust_data)
  num_vertices = max(flattened_data)

  # grab the data per cluster
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])

  # compute tightness of all clusters
  tightness_vector = sapply(clusters, function(x)
    compute_tightness(dists, x))

  return(tightness_vector)
}

#' Get data within a cluster
#'
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster ids
#'
#' @return A list of strings, each a comma separated list of the toString values of the data point names.
get_clustered_data <- function(binclust_data) {

  # no need to list by level set
  flattened_data = unlist(binclust_data)
  num_vertices = max(flattened_data)

  # grab the data per cluster
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])

  # make a big list of the names of the data in the cluster
  data_in_cluster = lapply(lapply(clusters, names), toString)

  return(data_in_cluster)
}

# edge data --------------------------------------------------------------

#' Calculate edge weights
#'
#' @param overlap_lengths A named vector of cluster overlap lengths, obtained by calling [length()] on the output from \code{[get_overlaps()]}.
#' @param cluster_sizes A vector of cluster sizes.
#' @param edges A 2D array of source and target nodes, representing an edge list. Should be ordered consistently with the `overlap_lengths` parameter.
#'
#' @return A vector of real numbers representing the Jaccard index of each overlap.
#' @details This value is calculated per edge by dividing the number of data points in the union of the two clusters by the number of data points in the intersection. Formally, \deqn{w(\{c_i, c_j\}) = \dfrac{|c_i \cap c_j|}{|c_i \cup c_j|} = \dfrac{|c_i \cap c_j|}{|c_i| + |c_j| - |c_i \cap c_j|}}
get_edge_weights <- function(overlap_lengths, cluster_sizes, edges) {

  # no edges? no weights
  if (length(edges) == 0) {
    return(NULL)
  }

  # grab source and target cluster sizes
  heads = edges[, 1]
  tails = edges[, 2]
  head_sizes = cluster_sizes[heads]
  tail_sizes = cluster_sizes[tails]

  jaccards = overlap_lengths / (head_sizes + tail_sizes - overlap_lengths)

  return(jaccards)
}

