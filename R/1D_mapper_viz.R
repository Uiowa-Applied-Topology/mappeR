source("R/1D_mapper.R")

get_cluster_tightness_vector <- function(dists, binclust_data, num_vertices) {
  flattened_data = unlist(binclust_data) # we don't care about bins here
  res = c()

  # there is some lapply way to do this but i actually don't know how to code
  for (i in 1:num_vertices) {
    min = Inf
    min_name = NULL
    my_cluster = flattened_data[flattened_data == i] # find all the points with the same cluster number
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

# calculates how large nodes should be based on how many datapoints it contains
# returns a vector of cluster sizes, normalized and rescaled
get_size_vector <- function(binclust_data, num_vertices) {
  size_vector = c()

  flattened_data = unlist(binclust_data) # we don't care about bins here

  # there is some lapply way to do this but i actually don't know how to code
  for (i in 1:num_vertices) {
    my_cluster = flattened_data[flattened_data == i] # find all the points with the same cluster number
    size_vector = append(size_vector, length(my_cluster)) # record that cluster's length
  }

  scale_factor = 1000 # change according to your eyeballs/monitor

  size_vector = (size_vector/sqrt(sum(size_vector^2)))*scale_factor

  return(size_vector)
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

  cygraph = set_vertex_attr(mapper_graph, "bin", value = bin_vector)
  cygraph = set_vertex_attr(cygraph, "cluster", value = 1:num_vertices)
  cygraph = set_vertex_attr(cygraph, "cluster_tightness", value = tightness_vector)
  cygraph = set_edge_attr(cygraph, "overlap", value = (edge_weights/sqrt(sum(edge_weights^2)))*25)
  cygraph = set_vertex_attr(cygraph, "cluster_size", value = get_size_vector(binclust_data, num_vertices))

  createNetworkFromIgraph(cygraph)

  style.name = "mapperstyle"
  defaults <- list(NODE_SHAPE = "ellipse",
                   NODE_FILL_COLOR = "#FFFFFF",
                   NODE_BORDER_WIDTH = 10,
                   NODE_BORDER_TRANSPARENCY = 0,
                   EDGE_TRANSPARENCY = 120)

  nodeLabels <- mapVisualProperty('node label', 'cluster', 'p')
  nodeSizes <- mapVisualProperty('node size', 'cluster_size', 'p')
  edgeWidth <- mapVisualProperty('edge width', 'overlap', 'p')

  createVisualStyle(style.name, defaults, list(nodeLabels, nodeSizes, edgeWidth))

  setVisualStyle(style.name)

  setNodeBorderWidthDefault(10, style.name = style.name)
  setNodeBorderColorMapping("bin", c(1, num_bins/2, num_bins), c("#a50026", "#ffffbf", "#313695"), style.name = style.name)
  setNodeColorMapping("cluster_tightness", c(0,.5,1), c("#ffffff", "#8d8d8d", "#000000"), style.name = style.name)
}

cymapper <- function(data, filtered_data, dists, num_bins, percent_overlap, clustering_method) {
  visualize_mapper_data(get_mapper_data(data, filtered_data, dists, num_bins, percent_overlap, clustering_method), dists)
  return(invisible(NULL))
}
