source("R/1D_mapper.R")
source("R/combomapper_viz.R")

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

visualize_mapper_data <- function(mapper_data, dists) {
  # get all of our data
  binclust_data = mapper_data[[1]]
  mapper_graph = mapper_data[[2]]
  edge_weights = mapper_data[[3]]

  num_vertices = gorder(mapper_graph)
  num_edges = gsize(mapper_graph)
  num_bins = length(binclust_data)
  bin_vector = get_bin_vector(binclust_data)
  tightness_vector = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_sizes = get_size_vector(binclust_data, num_vertices)
  bin_colors = palette(rainbow(num_bins))

  cygraph = set_vertex_attr(mapper_graph, "bin", value = bin_vector)
  cygraph = set_vertex_attr(cygraph, "cluster", value = 1:num_vertices)
  cygraph = set_vertex_attr(cygraph, "cluster_tightness", value = tightness_vector)
  cygraph = set_edge_attr(cygraph, "overlap", value = get_edge_weights(edge_weights, cluster_sizes, ends(mapper_graph, E(mapper_graph)))*10)
  cygraph = set_vertex_attr(cygraph, "cluster_size", value = 100*cluster_sizes/max(cluster_sizes))

  createNetworkFromIgraph(cygraph)

  style.name = paste("mapperstyle", runif(1))
  defaults <- list(NODE_SHAPE = "ellipse",
                   NODE_BORDER_WIDTH = 10,
                   EDGE_TRANSPARENCY = 255)

  nodeSizes <- mapVisualProperty('node size', 'cluster_size', 'p')
  edgeWidth <- mapVisualProperty('edge width', 'overlap', 'p')
  nodeBorderColors <- mapVisualProperty('node border color', 'bin', 'd', 1:num_bins, bin_colors)
  nodeFillColors <- mapVisualProperty('node fill color', 'cluster_tightness', 'c', c(0, mean(tightness_vector), 1), c("#ffffff", "#efefef", "#000000"))

  createVisualStyle(style.name, defaults, list(nodeSizes, edgeWidth, nodeBorderColors, nodeFillColors))

  setVisualStyle(style.name)
}

cymapper <- function(data, filtered_data, dists, num_bins, percent_overlap, clustering_method) {
  visualize_mapper_data(get_mapper_data(data, filtered_data, dists, num_bins, percent_overlap, clustering_method), dists)
  return(invisible(NULL))
}
