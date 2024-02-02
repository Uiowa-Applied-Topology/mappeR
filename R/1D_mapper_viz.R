get_bin_vector <- function(clustered_data) {
  clusters_and_bins = c()
  for (i in 1:length(clustered_data)) {
    for (j in unique(clustered_data[[i]])) {
      clusters_and_bins = append(clusters_and_bins, i)
    }
  }
  return(clusters_and_bins)
}

get_size_vector <- function(clustered_data, num_vertices) {
  size_vector = c()

  flattened_data = unlist(clustered_data)

  for (i in 1:num_vertices) {
    my_cluster = flattened_data[flattened_data == i]
    size_vector = append(size_vector, length(my_cluster))
  }

  size_vector = (size_vector/sqrt(sum(size_vector^2)))*200

  return(size_vector)
}

get_color_vector <- function(bin_vector, num_vertices, num_bins) {
  colfun = colorRampPalette(c('#998ec3', '#f7f7f7', '#f1a340'))
  colors = colfun(num_bins)
  color_vector = c()

  for (i in 1:num_vertices) {
    color_vector = append(color_vector, colors[bin_vector[i]])
  }

  return(color_vector)
}

visualize_mapper_data <- function(mapper_data) {
  clustered_data = mapper_data[[1]]
  mapper_graph = mapper_data[[2]]
  edge_weights = mapper_data[[3]]

  num_vertices = gorder(mapper_graph)
  num_edges = gsize(mapper_graph)
  num_bins = length(clustered_data)
  bin_vector = get_bin_vector(clustered_data)

  cygraph = set_vertex_attr(mapper_graph, "bin", value = bin_vector)
  cygraph = set_vertex_attr(cygraph, "cluster", value = 1:num_vertices)
  cygraph = set_edge_attr(cygraph, "overlap", value = (edge_weights/sqrt(sum(edge_weights^2)))*25)
  cygraph = set_vertex_attr(cygraph, "cluster_size", value = get_size_vector(clustered_data, num_vertices))
  cygraph = set_vertex_attr(cygraph, "color", value = get_color_vector(bin_vector, num_vertices, num_bins))

  createNetworkFromIgraph(cygraph)

  style.name = "mapperstyle"
  defaults <- list(NODE_SHAPE = "ellipse",
                   EDGE_TRANSPARENCY = 120)

  nodeLabels <- mapVisualProperty('node label', 'cluster', 'p')
  nodeSizes <- mapVisualProperty('node size', 'cluster_size', 'p')
  edgeWidth <- mapVisualProperty('edge width', 'overlap', 'p')

  createVisualStyle(style.name, defaults, list(nodeLabels, nodeSizes, edgeWidth))

  setVisualStyle(style.name)

  setNodeColorMapping("bin", c(1, num_bins/2, num_bins), c("#998ec3", "#f7f7f7", "#f1a340"), style.name = style.name)
}

cymapper <- function(data, filtered_data, dists, num_bins, percent_overlap, clustering_method) {
  visualize_mapper_data(get_mapper_data(data, filtered_data, dists, num_bins, percent_overlap, clustering_method))
  return(invisible(NULL))
}
