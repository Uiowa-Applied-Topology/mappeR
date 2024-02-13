source("R/ballmapper.R")
source("R/1D_mapper_viz.R")

get_size_vector <- function(binclust_data, num_vertices) {
  size_vector = c()

  flattened_data = unlist(binclust_data)

  for (i in 1:num_vertices) {
    my_cluster = flattened_data[flattened_data == i]
    size_vector = append(size_vector, length(my_cluster))
  }

  return(size_vector)
}

visualize_ballmapper_data <- function(mapper_data, dists) {
  binclust_data = mapper_data[[1]]
  mapper_graph = mapper_data[[2]]
  overlap_data = mapper_data[[3]]

  num_vertices = gorder(mapper_graph)
  num_edges = gsize(mapper_graph)
  num_bins = length(binclust_data)
  tightness_vector = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_sizes = get_size_vector(binclust_data, num_vertices)
  edge_weights = get_edge_weights(overlap_data, cluster_sizes, ends(mapper_graph, E(mapper_graph)))

  cygraph = set_vertex_attr(mapper_graph, "cluster", value = 1:num_vertices)
  cygraph = set_edge_attr(cygraph, "overlap", value = edge_weights)
  cygraph = set_vertex_attr(cygraph, "cluster_size", value = 100*cluster_sizes/max(cluster_sizes))
  cygraph = set_vertex_attr(cygraph, "cluster_tightness", value = tightness_vector)
  cygraph = set.vertex.attribute(cygraph, "points in cluster", value=cluster_sizes)

  createNetworkFromIgraph(cygraph)

  style.name = paste("mapperstyle", runif(1))
  defaults <- list(NODE_SHAPE = "ellipse",
                   EDGE_TRANSPARENCY = 120)

  nodeSizes <- mapVisualProperty('node size', 'cluster_size', 'p')
  edgeWidth <- mapVisualProperty('edge width', 'overlap', 'p')
  nodeFillColor <- mapVisualProperty('node fill color', 'cluster_tightness', 'c', c(0, mean(tightness_vector), 1), c("#ffffff", "#fefefe", "#000000"))

  createVisualStyle(style.name, defaults, list(nodeSizes, edgeWidth, nodeFillColor))

  setVisualStyle(style.name)

}

cyballmapper <- function(data, dists, eps) {
  visualize_ballmapper_data(get_ballmapper_data(data, dists, eps), dists)
  return(invisible(NULL))
}
