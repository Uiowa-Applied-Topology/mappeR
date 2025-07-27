###########################################################################
# CLUSTER FACTORY
# cluster assembly and mapper graph statistics
###########################################################################

#' Distance Matrix Splicer
#'
#' Subset a `dist` object.
#'
#' @param patch A list of names of data points.
#' @param dists A `dist` object for data points in the patch, possibly including extra points.
#' @noRd
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

#' Unvarnished Clustering
#'
#' This function tells the computer to look away for a second, so the goblins come and cluster your data while it's not watching.
#'
#' @param dist_mats A list of distance matrices of each bin that is to be clustered. Needs to be acceptable to `clusterer`.
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix;
#' that is, it should return a list of named vectors, whose names are the names of data points and whose values are cluster assignments (integers).
#'
#' @return The output of `clusterer(dist_mats)`, which needs to be a list containing named vectors (one per bin), whose names are data point names and whose values are cluster labels (within each bin)
#' @noRd
get_raw_clusters <- function(dist_mats, clusterer) {
  return(clusterer(dist_mats))
}

#' Varnished Clustering
#'
#' Process level sets of data and a global distance matrix to return fresh clusters.
#'
#' @param bins A `list` containing "bins" of vectors of names of data points.
#' @param dists A distance matrix containing pairwise distances between named data points. Needs to be acceptable to `clusterer`.
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix.
#'
#' @return The output of `clusterer` applied to a list of distance matrices, which should be a list containing named vectors (one per bin), whose names are data point names and whose values are cluster labels.
#' These labels are unique to the bin which has them.
#' @noRd
get_clusters <- function(bins, dists, clusterer) {
  # more than one bin, need more than one distance matrix
  if (is.list(bins)) {
    # subset the global distance matrix per bin
    dist_mats = lapply(bins, subset_dists, dists = dists)

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

#' Cluster Weigher
#'
#' Compute cluster sizes.
#'
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster IDs (integers).
#'
#' @return A vector of integers representing the lengths of the clusters in the input data.
#' @noRd
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

#' Patch Identifier
#'
#' Recover patch membership from cluster assignment list.
#'
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster ids.
#'
#' @return A vector of integers equal in length to the number of clusters, whose members identify which bin that cluster belongs to.
#' @noRd
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

#' Medoid Getter
#'
#' Get the medoid of a single cluster.
#'
#' @param dists A distance matrix for points in the cluster.
#' @param cluster A named vector whose names are data point names and whose values are cluster labels.
#'
#' @return The medoid of the `cluster`, which is the data point in that cluster with the minimum sum of distances to all other data points.
#' @noRd
get_cluster_medoid <- function(dists, cluster) {
  # singleton clusters are their own medoids
  if (length(cluster) == 0) {
    return(0)
  } else if (length(cluster) == 1) {
    # print(cluster)
    return(names(cluster))
  } else {
    # get the distances associated to points in this cluster
    cluster_names = names(cluster)
    these_dists = dists[cluster_names, cluster_names]

    # find minimum sum of distances
    sums = apply(as.matrix(these_dists), 1, sum)
    # names(sums) = cluster_names
    min_sum = min(sums)
    medoid = sums[sums == min_sum]
    if (length(medoid) > 1) {
      rand = sample(1:length(medoid), 1)
      medoid = medoid[rand]
    }
    return(names(medoid))
  }
}

#' Medoids Getter
#'
#' Get the medoids of a list of clusters.
#'
#' @param dists A distance matrix for the data points inside all the input clusters.
#' @param binclust_data A list of named vectors whose names are those of data points and whose values are cluster IDs (integers).
#'
#' @return A vector of names of the data points in each cluster with the minimum sum of distances to all other within-cluster data points (the medoid).
#' @noRd
get_cluster_medoids <- function(dists, binclust_data) {

  # no need to list by level set
  flattened_data = unlist(binclust_data)
  num_vertices = max(flattened_data)

  # grab the data per cluster
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])

  # compute tightness of all clusters
  medoids = sapply(clusters, function(x)
    get_cluster_medoid(dists, x))

  return(medoids)
}

#' Tightness Calculator
#'
#' Compute a measure of dispersion for a single cluster.
#'
#' @param dists A distance matrix for points in the cluster.
#' @param cluster A named vector whose names are data point names and whose values are cluster labels.
#'
#' @return A real number representing the mean distance to the medoid of the cluster, which is the data point with the smallest combined distance to every other point.
#' A smaller value indicates a tighter cluster based on this measure.
#' @noRd
compute_tightness <- function(dists, cluster) {

  # empty or singleton clusters have trivial tightness
  if (length(cluster) <= 1) {
    return(0)
  } else {
    # get the distances associated to points in this cluster
    cluster_names = names(cluster)
    these_dists = dists[cluster_names, cluster_names]
    medoid = get_cluster_medoid(dists, cluster)
    min_dists = these_dists[medoid, ]
    closeness_factor = sum(min_dists) / (length(cluster) - 1)

    return(closeness_factor)
  }
}

#' Tightnesses Calculator
#'
#' @param dists A distance matrix for the data points inside all the input clusters.
#' @param binclust_data A list of named vectors whose names are those of data points and whose values are cluster IDs (integers).
#'
#' @return A vector of real numbers containing mean distances to the medoids of each cluster in `dists`.
#' @noRd
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

#' Max Distance Calculator
#'
#' Compute the maximum distance from the medoid of a single cluster.
#'
#' @param dists A distance matrix for points in the cluster.
#' @param cluster A named vector whose names are data point names and whose values are cluster labels.
#'
#' @return A real number representing the maximum distance to the medoid within the input cluster.
#' @noRd
get_max_eccentricity <- function(dists, cluster) {
  # empty or singleton clusters have trivial tightness
  if (length(cluster) <= 1) {
    return(0)
  } else {
    # get the distances associated to points in this cluster
    cluster_names = names(cluster)
    these_dists = dists[cluster_names, cluster_names]
    medoid = get_cluster_medoid(dists, cluster)
    min_dists = these_dists[medoid, ]
    max_dist = max(min_dists)

    return(max_dist)
  }
}

#' Max Distances Calculator
#'
#' @param dists A distance matrix for the data points inside all the input clusters.
#' @param binclust_data A list of named vectors whose names are those of data points and whose values are cluster IDs (integers).
#'
#' @return A vector of real numbers in containing maximum distances to the medoids of each cluster in `dists`.
#' @noRd
get_max_eccentricities <- function(dists, binclust_data) {

  # no need to list by level set
  flattened_data = unlist(binclust_data)
  num_vertices = max(flattened_data)

  # grab the data per cluster
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])

  # compute tightness of all clusters
  max_vector = sapply(clusters, function(x)
    get_max_eccentricity(dists, x))

  return(max_vector)
}

#' WCSS Calculator
#'
#' Compute the WCSS for a single cluster.
#'
#' @param dists A distance matrix for points in the cluster.
#' @param cluster A named vector whose names are data point names and whose values are cluster labels.
#'
#' @return A real number representing the sum of squares of distances to the medoid of a single cluster.
#' @noRd
get_wcss <- function(dists, cluster) {

  # empty or singleton clusters have trivial tightness
  if (length(cluster) <= 1) {
    return(0)
  } else {
    # get the distances associated to points in this cluster
    cluster_names = names(cluster)
    these_dists = dists[cluster_names, cluster_names]
    medoid = get_cluster_medoid(dists, cluster)
    min_dists = these_dists[medoid, ]
    sum_of_squares = sum(min_dists^2)

    return(sum_of_squares)
  }
}

#' WCSSs Calculator
#'
#' @param dists A distance matrix for the data points inside all the input clusters.
#' @param binclust_data A list of named vectors whose names are those of data points and whose values are cluster IDs (integers).
#'
#' @return A vector of real numbers in containing sums of squares of distances to the medoids of each cluster in `dists`.
#' @noRd
get_all_wcss <- function(dists, binclust_data) {

  # no need to list by level set
  flattened_data = unlist(binclust_data)
  num_vertices = max(flattened_data)

  # grab the data per cluster
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])

  # compute tightness of all clusters
  wcss_vector = sapply(clusters, function(x)
    get_wcss(dists, x))

  return(wcss_vector)
}


#' Cluster Width Calculator
#'
#' Compute the maximum pairwise distance for a single cluster.
#'
#' @param dists A distance matrix for points in the cluster.
#' @param cluster A named vector whose names are data point names and whose values are cluster labels.
#'
#' @return A real number representing the maximum pairwise distance within a single cluster.
#' @noRd
get_width <- function(dists, cluster) {

  # empty or singleton clusters have trivial tightness
  if (length(cluster) <= 1) {
    return(0)
  } else {
    # get the distances associated to points in this cluster
    cluster_names = names(cluster)
    these_dists = dists[cluster_names, cluster_names]

    return(max(these_dists))
  }
}

#' Cluster Widths Calculator
#'
#' @param dists A distance matrix for the data points inside all the input clusters.
#' @param binclust_data A list of named vectors whose names are those of data points and whose values are cluster IDs (integers).
#'
#' @return A vector of real numbers in containing maximum pairwise distances within each cluster in `dists`.
#' @noRd
get_widths <- function(dists, binclust_data) {

  # no need to list by level set
  flattened_data = unlist(binclust_data)
  num_vertices = max(flattened_data)

  # grab the data per cluster
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])

  # compute tightness of all clusters
  maxes = sapply(clusters, function(x)
    get_width(dists, x))

  return(maxes)
}
#' Cluster Manifesto Logger
#'
#' Get the data names for data within a cluster.
#'
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster ids
#'
#' @return A `list` of strings, each a comma separated list of the `toString` values of the data point names.
#' @noRd
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

#' Edge Weigher
#'
#' Calculate Jaccard indices to represent edge weights.
#'
#' @param overlap_lengths A named vector of cluster overlap lengths, obtained by calling [length()] on the output from \code{[get_overlaps()]}.
#' @param cluster_sizes A vector of cluster sizes.
#' @param edges A 2D array of source and target nodes, representing an edge list. Should be ordered consistently with the `overlap_lengths` parameter.
#'
#' @return A vector of real numbers representing the Jaccard index of each overlap.
#' @details This value is calculated per edge by dividing the number of data points in the union of the two clusters by the number of data points in the intersection. Formally, \deqn{w(\{c_i, c_j\}) = \dfrac{|c_i \cap c_j|}{|c_i \cup c_j|} = \dfrac{|c_i \cap c_j|}{|c_i| + |c_j| - |c_i \cap c_j|}}
#' @noRd
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

