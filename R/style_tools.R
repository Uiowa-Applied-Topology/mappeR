get_size_vector <- function(binclust_data, num_vertices) {
  size_vector = c()

  flattened_data = unlist(binclust_data)

  for (i in 1:num_vertices) {
    my_cluster = flattened_data[flattened_data == i]
    size_vector = append(size_vector, length(my_cluster))
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

get_cluster_tightness_vector <- function(dists, binclust_data, num_vertices) {
  flattened_data = unlist(binclust_data) # we don't care about bins here
  res = c()

  # there is some lapply way to do this but i actually don't know how to code
  for (i in 1:num_vertices) {
    min = Inf
    min_name = NULL
    my_cluster = flattened_data[flattened_data == i] # find all the points with the same cluster number
    if (length(my_cluster) == 1) { # trivially tight cluster
      res = append(res, 1)
    } else {
      # find the medoid of the dataset
      for (datapoint in names(my_cluster)) {
        these_dists = dists[datapoint,names(my_cluster)]
        dist_sum = sum(these_dists)
        if (dist_sum < min) {
          min = dist_sum
          min_name = datapoint
        }
      }
      close_datapoint_dists = dists[min_name,names(my_cluster)]
      n = length(close_datapoint_dists)
      min_dist = min(close_datapoint_dists[close_datapoint_dists > 0])
      dist_sums = sum(close_datapoint_dists)
      closeness_factor = (length(my_cluster)-1)*min_dist/dist_sums
      res = append(res, closeness_factor)
    }
  }
  return(res)
}

get_edge_weights <- function(overlap_lengths, cluster_sizes, edges) {
  edge_weights = c()
  for (i in 1:length(overlap_lengths)) {
    head_id = edges[i, 1]
    tail_id = edges[i, 2]
    node_size_1 = cluster_sizes[head_id]
    node_size_2 = cluster_sizes[tail_id]
    overlap = overlap_lengths[i]
    overlap_weight = max(overlap/node_size_1, overlap/node_size_2)
    edge_weights = append(edge_weights, overlap_weight)
  }
  return(edge_weights)
}
