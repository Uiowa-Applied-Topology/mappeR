###########################################################################
# STATS HUB
# cluster stats
###########################################################################


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
#' @return A real number in \eqn{[0,1]} representing a measure of dispersion of a cluster.
#' @details This method computes a measure of cluster dispersion. It finds the medoid of the input data set and returns the sum of distances from the medoid scaled by the size of the data set. Formally, we say the tightness \eqn{\tau} of a cluster \eqn{C} is given by \deqn{\tau(C) = \dfrac{\displaystyle\sum_{i}\text{dist}(x_i, x_j)}{\left(|C| - 1\right)}} where \deqn{x_j = \text{arg}\,\min\limits_{x_j\in C}\, \sum_{x_i \in C, i\neq j}\text{dist}(x_i, x_j)} A smaller value indicates a tighter cluster based on this metric.
compute_tightness <- function(dists, cluster) {

  # empty, singleton, or doubleton clusters have trivial tightness
  if (length(cluster) <= 2) {
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
#' @return A vector of real numbers in \eqn{(0,\infty)} representing a measure of dispersion of a cluster, calculated according to [compute_tightness].
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
#' @return A vector of real numbers representing cluster overlap strength.
#' @details This value is calculated per edge by dividing the number of data points in the overlap by the number of points in the cluster on either end, and taking the maximum value. Formally, \deqn{w(\{c_i, c_j\}) = \displaystyle\max\left\{\dfrac{|c_i \cap c_j|}{|c_i|}, \dfrac{|c_i\cap c_j|}{|c_j|}\right\}}
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
  total_size = head_sizes + tail_sizes

  # compute edge weight as maximum relative overlap
  head_overlaps = overlap_lengths / total_size
  tail_overlaps = overlap_lengths / total_size
  edge_weights = mapply(max, head_overlaps, tail_overlaps)

  return(edge_weights)
}
