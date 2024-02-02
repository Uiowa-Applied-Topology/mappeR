source("R/1D_mapper.R")

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

  scale_factor = 200 # change according to your eyeballs/monitor

  size_vector = (size_vector/sqrt(sum(size_vector^2)))*scale_factor

  return(size_vector)
}

visualize_mapper_data <- function(mapper_data) {
  # get all of our data
  binclust_data = mapper_data[[1]]
  mapper_graph = mapper_data[[2]]
  edge_weights = mapper_data[[3]]

  num_vertices = gorder(mapper_graph)
  num_edges = gsize(mapper_graph)
  num_bins = length(binclust_data)
  bin_vector = get_bin_vector(binclust_data)

  cygraph = set_vertex_attr(mapper_graph, "bin", value = bin_vector)
  cygraph = set_vertex_attr(cygraph, "cluster", value = 1:num_vertices)
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
  setNodeBorderColorMapping("bin", c(1, num_bins/2, num_bins), c("#0f62fe", "#0072c3", "#004144"), style.name = style.name)
  # setNodeBorderColorMapping("bin", colors = t(get_color_vector(bin_vector, num_vertices, num_bins)))
}

cymapper <- function(data, filtered_data, dists, num_bins, percent_overlap, clustering_method) {
  visualize_mapper_data(get_mapper_data(data, filtered_data, dists, num_bins, percent_overlap, clustering_method))
  return(invisible(NULL))
}
