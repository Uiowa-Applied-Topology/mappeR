source("R/ballmapper.R")
source("R/1D_mapper_viz.R")

get_size_vector <- function(binclust_data, num_vertices) {
  size_vector = c()

  flattened_data = unlist(binclust_data)

  for (i in 1:num_vertices) {
    my_cluster = flattened_data[flattened_data == i]
    size_vector = append(size_vector, length(my_cluster))
  }

  size_vector = (size_vector/sqrt(sum(size_vector^2)))*200

  return(size_vector)
}

visualize_ballmapper_data <- function(mapper_data, dists) {
  binclust_data = mapper_data[[1]]
  mapper_graph = mapper_data[[2]]
  edge_weights = mapper_data[[3]]

  num_vertices = gorder(mapper_graph)
  num_edges = gsize(mapper_graph)
  num_bins = length(binclust_data)
  tightness_vector = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)

  cygraph = set_vertex_attr(mapper_graph, "cluster", value = 1:num_vertices)
  cygraph = set_edge_attr(cygraph, "overlap", value = (edge_weights/sqrt(sum(edge_weights^2)))*25)
  cygraph = set_vertex_attr(cygraph, "cluster_size", value = get_size_vector(binclust_data, num_vertices))
  cygraph = set_vertex_attr(cygraph, "cluster_tightness", value = tightness_vector)

  createNetworkFromIgraph(cygraph)

  style.name = "mapperstyle"
  defaults <- list(NODE_SHAPE = "ellipse",
                   EDGE_TRANSPARENCY = 120)

  nodeLabels <- mapVisualProperty('node label', 'cluster', 'p')
  nodeSizes <- mapVisualProperty('node size', 'cluster_size', 'p')
  edgeWidth <- mapVisualProperty('edge width', 'overlap', 'p')

  createVisualStyle(style.name, defaults, list(nodeLabels, nodeSizes, edgeWidth))

  setVisualStyle(style.name)

  setNodeColorMapping("cluster_tightness", c(0,.5,1), c("#ffffff", "#8d8d8d", "#000000"), style.name = style.name)

}

cyballmapper <- function(data, dists, eps) {
  visualize_ballmapper_data(get_ballmapper_data(data, dists, eps), dists)
  return(invisible(NULL))
}
