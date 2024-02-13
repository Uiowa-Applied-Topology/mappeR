source("R/combomapper.R")
library(igraph)

get_edge_weights <- function(overlap_lengths, cluster_sizes, edges) {
  edge_weights = c()
  for (i in 1:length(overlap_lengths)) {
    head_id = edges[i, 1]
    tail_id = edges[i, 2]
    node_size_1 = cluster_sizes[head_id]
    node_size_2 = cluster_sizes[tail_id]
    overlap = overlap_lengths[i]
    overlap_weight = ((overlap/node_size_1) + (overlap/node_size_2))/2
    edge_weights = append(edge_weights, overlap_weight)
  }
  return(edge_weights)
}

visualize_combomapper_data <- function(mapper_data, dists) {
  # get all of our data
  binclust_data = mapper_data[[1]]
  mapper_graph = mapper_data[[2]]
  overlap_data = mapper_data[[3]]

  num_vertices = gorder(mapper_graph)
  num_edges = gsize(mapper_graph)
  num_bins = length(binclust_data)
  bin_vector = get_bin_vector(binclust_data)
  tightness_vector = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)

  cluster_sizes = get_size_vector(binclust_data, num_vertices)
  edge_weights = get_edge_weights(overlap_data, cluster_sizes, ends(mapper_graph, E(mapper_graph)))

  cygraph = set_vertex_attr(mapper_graph, "bin", value = bin_vector)
  cygraph = set_vertex_attr(cygraph, "cluster", value = 1:num_vertices)
  cygraph = set_vertex_attr(cygraph, "cluster_tightness", value = tightness_vector)
  cygraph = set_edge_attr(cygraph, "overlap", value = edge_weights*255)
  cygraph = set_vertex_attr(cygraph, "cluster_size", value = cluster_sizes/sqrt(sum(cluster_sizes^2))*200)

  createNetworkFromIgraph(cygraph)

  style.name = paste("mapperstyle", runif(1))
  defaults <- list(NODE_SHAPE = "ellipse",
                   NODE_BORDER_WIDTH = 10,
                   BORDER_TRANSPARENCY = 255)

  nodeSizes <- mapVisualProperty('node size', 'cluster_size', 'p')
  edgeTransparency <- mapVisualProperty('edge transparency', 'overlap', 'p')
  nodeBorderColors = mapVisualProperty('node border color', 'bin', 'c', c(1, num_bins/2, num_bins), c("#d9fbfb", "#008d89", "#081a1c"))
  nodeFillColors = mapVisualProperty('node fill color', 'cluster_tightness', 'c', c(0, .5, 1), c("#ffffff", "#efefef", "#000000"))

  createVisualStyle(style.name, defaults, list(nodeSizes, edgeTransparency, nodeBorderColors, nodeFillColors))

  setVisualStyle(style.name)
}

cycombomapper <- function(data, dist1, dist2, eps) {
  visualize_combomapper_data(get_combomapper_data(data, dist1, dist2, eps), dist2)

  return(invisible(NULL))
}
