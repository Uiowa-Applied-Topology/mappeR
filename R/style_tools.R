# given our mappered data, returns the a vector with the number of datapoints per cluster.
get_cluster_sizes <- function(binclust_data, num_vertices) {
  size_vector = c()

  flattened_data = unlist(binclust_data) # bins are not important here
  clusters = lapply(1:num_vertices, function(x) flattened_data[flattened_data == x])

  return(sapply(clusters, length))
}

# categorizes clusters by bin. returns a vector of bin numbers corresponding to each cluster.
get_bin_vector <- function(binclust_data) {
  num_unique_clusters_per_bin = sapply(lapply(binclust_data, unique), length)
  bin_by_clusters = unlist(mapply(function(x, y) rep(x, y), 1:length(num_unique_clusters_per_bin), num_unique_clusters_per_bin))

  return(bin_by_clusters)
}

compute_tightness <- function(dists, cluster) {
  if (length(cluster) == 0) {
    return(1)
  } else {
    cluster_names = names(cluster)
    these_dists = dists[cluster_names, cluster_names]
    sums = apply(these_dists, 1, sum)
    min_sum = min(sums) # sum of distances of the medoid
    min_dists = these_dists[which(sums == min_sum),] # distance matrix subset of medoid
    min_dist = min(min_dists[min_dists > 0]) # distance of closest datapoint to medoid
    closeness_factor = (length(cluster)-1)*(min(min_dist))/min_sum # "normalized" so a score of 1 means the medoid is as close as possible to all the data in the cluster
    return(closeness_factor)
  }
}

# computes a measure of cluster "tightness" as a number between 0 and 1, based on the medoids of the clusters.
get_cluster_tightness_vector <- function(dists, binclust_data, num_vertices) {
  flattened_data = unlist(binclust_data) # we don't care about bins here

  clusters = lapply(1:num_vertices, function(x) flattened_data[flattened_data == x])
  tightness_vector = sapply(clusters, function(x) compute_tightness(dists, x))

  return(tightness_vector)
}

# calcuates edge weights. for each edge, its weight is the maximum of the relative percent overlaps between the clusters.
get_edge_weights <- function(overlap_lengths, cluster_sizes, edges) {
  heads = edges[,1]
  tails = edges[,2]
  head_sizes = cluster_sizes[heads]
  tail_sizes = cluster_sizes[tails]

  head_overlaps = overlap_lengths/head_sizes
  tail_overlaps = overlap_lengths/tail_sizes
  edge_weights = mapply(max, head_overlaps, tail_overlaps)

  return(edge_weights)
}

# gives which data is in each cluster
get_clustered_data <- function(binclust_data, num_vertices) {
  flattened_data = unlist(binclust_data) # bins are not important here
  clusters = lapply(1:num_vertices, function(x) flattened_data[flattened_data == x]) # sort by cluster
  data_in_cluster = lapply(lapply(clusters, names), toString)

  return(data_in_cluster)
}
