# goblin clustering mines -------------------------------------------------

#' Ship data off to the clustering goblins
#'
#' This function tells the computer to look away for a second, so the goblins come and cluster your data while it's not watching.
#'
#' @param dist_mats A list of distance matrices of each bin that is to be clustered.
#' @param method A string that suggests how the goblins will handle the data.
#'
#' @return A list containing named vectors (one per bin), whose names are datapoint names and whose values are cluster labels (within each bin)
run_cluster_machine <- function(dist_mats, method) {
  switch(tolower(method), "single" = return(get_single_linkage_clusters(dist_mats)))
}

#' Initate the clustering process
#'
#' This function processes the binned data and global distance matrix to return freshly clustered data.
#'
#' @param bins A list containing "bins" of vectors of names of datapoints.
#' @param dists A distance matrix containing pairwise distances between named datapoints.
#' @param method A clue!
#'
#' @return A list containing named vectors (one per bin), whose names are datapoint names and whose values are cluster labels
get_clusters <- function(bins, dists, method) {
  # subset the global distance matrix per bin
  dist_mats = sapply(1:length(bins), function(x)
    as.dist(as.matrix(dists)[bins[[x]], bins[[x]]]))

  # cluster the data
  clusters = run_cluster_machine(dist_mats, method)

  # accurately total up clusters
  clusters_per_bin = sapply(clusters, max)
  offset = c(0, cumsum(clusters_per_bin))
  clusters = mapply(function(x, y)
    x + y, clusters, offset[-length(offset)])
  return(clusters)
}

# TODO: add more clustering methods

## single linkage clustering -----------------------------------------------

# please don't ask
run_slink <- function(dist) {
  return(fastcluster::hclust(dist, "single"))
}

# performs single linkage clustering and outputs clusters based on a rough heuristic.
get_single_linkage_clusters <- function(dist_mats) {
  dends = lapply(dist_mats, run_slink)
  return(process_dendrograms(dends))
}

# TODO: add options for what clustering math is given to mapper graph

# goblin data processing center -------------------------------------------

# given our mappered data, returns the a vector with the number of datapoints per cluster.
get_cluster_sizes <- function(binclust_data, num_vertices) {
  size_vector = c()

  flattened_data = unlist(binclust_data) # bins are not important here
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])

  return(sapply(clusters, length))
}

# categorizes clusters by bin. returns a vector of bin numbers corresponding to each cluster.
get_bin_vector <- function(binclust_data) {
  num_unique_clusters_per_bin = sapply(lapply(binclust_data, unique), length)
  bin_by_clusters = unlist(mapply(
    function(x, y)
      rep(x, y),
    1:length(num_unique_clusters_per_bin),
    num_unique_clusters_per_bin
  ))

  return(bin_by_clusters)
}

# TODO: add rationale for this
compute_tightness <- function(dists, cluster) {
  if (length(cluster) == 0) {
    return(1)
  }
  if (length(cluster == 1)) {
    return(1)
  } else {
    cluster_names = names(cluster)
    these_dists = dists[cluster_names, cluster_names]
    sums = apply(these_dists, 1, sum)
    min_sum = min(sums) # sum of distances of the medoid
    min_dists = these_dists[which(sums == min_sum), ] # distance matrix subset of medoid
    min_dist = min(min_dists[min_dists > 0]) # distance of closest datapoint to medoid
    closeness_factor = (length(cluster) - 1) * (min(min_dist)) / min_sum # "normalized" so a score of 1 means the medoid is as close as possible to all the data in the cluster
    return(closeness_factor)
  }
}

# computes a measure of cluster "tightness" as a number between 0 and 1, based on the medoids of the clusters.
get_cluster_tightness_vector <- function(dists, binclust_data, num_vertices) {
  flattened_data = unlist(binclust_data) # we don't care about bins here

  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])
  tightness_vector = sapply(clusters, function(x)
    compute_tightness(dists, x))

  return(tightness_vector)
}

# calcuates edge weights. for each edge, its weight is the maximum of the relative percent overlaps between the clusters.
get_edge_weights <- function(overlap_lengths, cluster_sizes, edges) {
  heads = edges[, 1]
  tails = edges[, 2]
  head_sizes = cluster_sizes[heads]
  tail_sizes = cluster_sizes[tails]

  head_overlaps = overlap_lengths / head_sizes
  tail_overlaps = overlap_lengths / tail_sizes
  edge_weights = mapply(max, head_overlaps, tail_overlaps)

  return(edge_weights)
}

# gives which data is in each cluster
get_clustered_data <- function(binclust_data, num_vertices) {
  flattened_data = unlist(binclust_data) # bins are not important here
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x]) # sort by cluster
  data_in_cluster = lapply(lapply(clusters, names), toString)

  return(data_in_cluster)
}
