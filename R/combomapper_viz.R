source("R/combomapper.R")
library(igraph)

visualize_combomapper_data <- function(mapper_data, dists) {
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

cycombomapper <- function(data, dist1, dist2, eps) {
  visualize_combomapper_data(get_combomapper_data(data, dist1, dist2, eps), dist2)

  return(invisible(NULL))
}
