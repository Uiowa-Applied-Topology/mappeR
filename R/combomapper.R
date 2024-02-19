# given binned data and the full data's distance matrix, clusters within each bin. keeps track of total clusters across bins.
# output is a list of named vectors; there is one named vector of data per bin containing a cluster number.
get_comboclusters <- function(bins, dists, method) {
  binclust_data = list()
  cluster_count = 0

  for (i in 1:length(bins)){
    if (length(bins[[i]]) == 0) {
      binclust_data[[i]] = list()
      next
    }
    if (length(bins[[i]]) == 1) {
      binclust_data[[i]] = setNames(1, names(bins[[i]])) + cluster_count
      cluster_count = cluster_count + 1
      next
    }
    d = dists[bins[[i]],bins[[i]]]
    # d = dist_subset(as.dist(dists), names(bins[[i]])) # assuming your data has identifiers! also I am not sure of this function's behavior on non-dissimilarity matrices.
    clust = switch(tolower(method), "single" = get_single_linkage_clusters(d), stop("please tell george to do more clustering methods"))
    binclust_data[[i]] = clust + cluster_count
    cluster_count = cluster_count + max(clust) # update the total
  }
  return(binclust_data)
}

construct_combomappergraph <- function(binclust_data, dists) {
  num_vertices = max(binclust_data[[length(binclust_data)]])

  node_ids = as.character(1:num_vertices)

  overlaps = get_overlaps(binclust_data)
  edges = get_edgelist_from_overlaps(overlaps, num_vertices)
  sources = as.character(edges[,1])
  targets = as.character(edges[,2])

  cluster_tightness = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_size = get_cluster_sizes(binclust_data, num_vertices)
  data_in_cluster = unlist(get_clustered_data(binclust_data, num_vertices))
  edge_weights = get_edge_weights(sapply(overlaps, length), cluster_size, edges)
  bins = get_bin_vector(binclust_data)

  nodes = data.frame(id=node_ids,
                     size=cluster_size,
                     tightness=cluster_tightness,
                     data=data_in_cluster,
                     bin=bins)

  edges = data.frame(source=sources,
                     target=targets,
                     weight=edge_weights)


  return(list(nodes, edges))
}

# runner function for combo mapper; outputs bins, clusters, and the mapper graph.
get_combomapper_data <- function(data, dist1, dist2, eps) {
  print("binning...")
  balls = create_balls(data, dist1, eps)

  print("clustering...")
  binclust_data = get_comboclusters(balls, dist2, "single")

  print("making mapper graph...")
  combomappergraph = construct_combomappergraph(binclust_data, dist2)

  return(combomappergraph)
}

cycombomapper <- function(data, dist1, dist2, eps) {
  visualize_mapper_data(get_combomapper_data(data, dist1, dist2, eps), FALSE)

  return(invisible(NULL))
}
