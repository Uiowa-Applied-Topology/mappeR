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

# computes a measure of cluster "tightness" as a number between 0 and 1, based on the medoids of the clusters.
get_cluster_tightness_vector <- function(dists, binclust_data, num_vertices) {
  flattened_data = unlist(binclust_data) # we don't care about bins here
  res = c()

  # there is some lapply way to do this but i actually don't know how to code
  for (i in 1:num_vertices) { # go cluster by cluster
    min = Inf
    min_name = NULL
    my_cluster = flattened_data[flattened_data == i] # gather points in this cluster
    if (length(my_cluster) == 1) { # trivially tight cluster
      res = append(res, 1)
    } else {
      for (datapoint in names(my_cluster)) { # go through all the data in the cluster
        these_dists = dists[datapoint,names(my_cluster)]
        dist_sum = sum(these_dists) # sum the distances from this datapoint to all others in the cluster
        if (dist_sum < min) { # we are looking for the smallest such sum
          min = dist_sum
          min_name = datapoint
        }
      }

      close_datapoint_dists = dists[min_name,names(my_cluster)] # distance set of the medoid
      min_dist = min(close_datapoint_dists[close_datapoint_dists > 0]) # find the closest datapoint to the medoid (assuming we are living in at least a semimetric)
      dist_sums = sum(close_datapoint_dists) # minimum sum
      closeness_factor = (length(my_cluster)-1)*min_dist/dist_sums # "normalized" so a score of 1 means the medoid is as close as possible to all the data in the cluster
      res = append(res, closeness_factor)
    }
  }
  return(res)
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
