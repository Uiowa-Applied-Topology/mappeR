###########################################################################
# STATS HUB
# cluster stats
###########################################################################


# node data --------------------------------------------------------------


#' Compute cluster sizes
#'
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster
#'
#' @return A vector of integers representing the lengths of the clusters in the input data.
get_cluster_sizes <- function(binclust_data) {
  flattened_data = unlist(binclust_data) # bins are not important here

  num_vertices = max(flattened_data)

  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])

  return(sapply(clusters, length))
}

#' Recover bins
#'
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster ids.
#'
#' @return A vector of integers equal in length to the number of clusters, whose members identify which bin that cluster belongs to.
get_bin_vector <- function(binclust_data) {
  if (!is.list(binclust_data)) {
    return(rep(1, max(binclust_data)))
  }
  num_unique_clusters_per_bin = sapply(lapply(binclust_data, unique), length)
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
#' @return A real number in \eqn{(0,\infty)} representing a measure of dispersion of a cluster. This method finds the medoid of the input data set, the point with the smallest sum of distances to all other points, and returns that sum divided by the largest distance from the medoid to another point. Formally, we say the tightness \eqn{\tau} of a cluster \eqn{C} is given by \deqn{\tau(C) = \dfrac{1}{\displaystyle\max_{x_i\in C, i\neq j}{\text{dist}(x_i, x_j)}}\displaystyle\sum_{i}\text{dist}(x_i, x_j)} where \deqn{x_j = \text{arg}\,\min\limits_{x_j\in C}\, \sum_{x_i \in C, i\neq j}\text{dist}(x_i, x_j)} A smaller value indicates a tighter cluster based on this metric.
compute_tightness <- function(dists, cluster) {
  if ((length(cluster) == 0) | (length(cluster) == 1)) {
    return(0)
  } else {
    cluster_names = names(cluster)
    these_dists = dists[cluster_names, cluster_names]
    sums = apply(these_dists, 1, sum)
    min_sum = min(sums) # sum of distances of the medoid
    min_dists = these_dists[which(sums == min_sum), ] # distance matrix subset of medoid
    closeness_factor = min_sum / max(min_dists)
    return(closeness_factor)
  }
}

#' Compute dispersion measures of a list of clusters
#'
#' @param dists A distance matrix for the data points inside the input clusters
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster ids
#'
#' @return A vector of real numbers in \eqn{(0,\infty)} representing a measure of dispersion of a cluster, calculated according to [compute_tightness()]
get_cluster_tightness_vector <- function(dists, binclust_data) {
  flattened_data = unlist(binclust_data) # we don't care about bins here
  num_vertices = max(flattened_data)

  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])
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
  flattened_data = unlist(binclust_data) # bins are not important here

  num_vertices = max(flattened_data)

  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x]) # sort by cluster

  # TODO: put spaces in between data names?
  data_in_cluster = lapply(lapply(clusters, names), toString)
}

## edge data --------------------------------------------------------------

#' Calculate edge weights
#'
#' @param overlap_lengths A named vector of cluster overlap lengths, obtained by calling [length()] on the output from \code{[get_overlaps()]}.
#' @param cluster_sizes A vector of cluster sizes.
#' @param edges A 2D array of source and target nodes, representing an edge list. Should be ordered consistently with the `overlap_lengths` parameter.
#'
#' @return A vector of real numbers representing cluster overlap strength. This is calculated per edge by dividing the number of data points in the overlap by the number of points in the cluster on either end, and taking the maximum value.
get_edge_weights <- function(overlap_lengths, cluster_sizes, edges) {

  if (length(edges) == 0) {
    return(NULL)
  }

  heads = edges[, 1]
  tails = edges[, 2]
  head_sizes = cluster_sizes[heads]
  tail_sizes = cluster_sizes[tails]
  total_size = head_sizes + tail_sizes

  if (nrow(edges) == 1) {
    return(max(length(overlap_lengths)/total_size, length(overlap_lengths)/total_size))
  }

  head_overlaps = overlap_lengths / total_size
  tail_overlaps = overlap_lengths / total_size
  edge_weights = mapply(max, head_overlaps, tail_overlaps, SIMPLIFY = FALSE)

  return(edge_weights)
}
