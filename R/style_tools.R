# given our mappered data, returns the a vector with the number of datapoints per cluster.
get_cluster_sizes <- function(binclust_data, num_vertices) {
  size_vector = c()

  flattened_data = unlist(binclust_data) # bins are not important here

  for (i in 1:num_vertices) {
    my_cluster = flattened_data[flattened_data == i] # get just the data in this cluster
    size_vector = append(size_vector, length(my_cluster)) # record the number of datapoints
  }

  return(size_vector)
}

# categorizes clusters by bin. returns a vector of bin numbers corresponding to each cluster.
get_bin_vector <- function(binclust_data) {
  clusters_and_bins = c()
  for (i in 1:length(binclust_data)) { # go through each bin
    for (j in unique(binclust_data[[i]])) { # the bin contains clusters
      clusters_and_bins = append(clusters_and_bins, i) # tag each cluster with the bin number
    }
  }
  return(clusters_and_bins)
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
  edge_weights = c()
  for (i in 1:length(overlap_lengths)) {
    cluster1_id = edges[i, 1]
    cluster2_id = edges[i, 2]
    cluster1_size = cluster_sizes[cluster1_id]
    cluster2_size = cluster_sizes[cluster2_id]
    overlap = overlap_lengths[i]

    overlap_weight = max(overlap/cluster1_size, overlap/cluster2_size) # pick the bigger relative overlap
    edge_weights = append(edge_weights, overlap_weight)
  }
  return(edge_weights)
}

# gives which data is in each cluster
get_clustered_data <- function(binclust_data, num_vertices) {
  res = array()

  flattened_data = unlist(binclust_data) # bins are not important here

  for (i in 1:num_vertices) {
    my_cluster = flattened_data[flattened_data == i] # gather points in this cluster
    res = append(res, toString(names(my_cluster)))
  }

  return(res)
}
