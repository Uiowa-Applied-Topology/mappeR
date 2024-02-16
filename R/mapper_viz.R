source("R/style_tools.R")

visualize_mapper_data <- function(mapper_data, dists, is_ballmapper = TRUE) {
  # get all of our data
  binclust_data = mapper_data[[1]]
  mapper_graph = mapper_data[[2]]
  overlap_data = mapper_data[[3]]

  # graph basics
  num_vertices = gorder(mapper_graph)
  num_edges = gsize(mapper_graph)

  # both ballmapper and conventional mapper get these
  tightness_vector = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_sizes = get_cluster_sizes(binclust_data, num_vertices)
  clustered_data = get_clustered_data(binclust_data, num_vertices)
  edge_weights = get_edge_weights(overlap_data, cluster_sizes, ends(mapper_graph, E(mapper_graph)))

  # record data in igraph, will become tables in cytoscape
  cygraph = set_vertex_attr(mapper_graph, "cluster_id", value = 1:num_vertices)
  cygraph = set_vertex_attr(cygraph, "data_in_cluster", value=clustered_data[-1])
  cygraph = set_vertex_attr(cygraph, "num_points_in_cluster", value=cluster_sizes)
  cygraph = set_vertex_attr(cygraph, "cluster_tightness", value = tightness_vector)
  cygraph = set_edge_attr(cygraph, "overlap", value = edge_weights)

  # conventional mapper gets more style stuff
  if (!is_ballmapper) {
    bin_vector = get_bin_vector(binclust_data)
    cygraph = set_vertex_attr(cygraph, "bin", value = bin_vector)
  }

  createNetworkFromIgraph(cygraph)

  style.name = paste("mapperstyle", runif(1))
  defaults <- list(NODE_SHAPE = "ellipse",
                   NODE_BORDER_WIDTH = 10,
                   NODE_BORDER_PAINT = "#000",
                   EDGE_TRANSPARENCY = 255)

  nodeSizes <- mapVisualProperty('node size', 'cluster_id', 'd', 1:num_vertices, 100*cluster_sizes/max(cluster_sizes))
  edgeWidth <- mapVisualProperty('edge width', 'overlap', 'p')
  nodeFillColors <- mapVisualProperty('node fill color', 'cluster_tightness', 'c', c(0, mean(tightness_vector), 1), c("#ffffff", "#efefef", "#000000"))

  if (is_ballmapper) { # ballmapper needs no more styling
    createVisualStyle(style.name, defaults, list(nodeSizes, edgeWidth, nodeFillColors))
  } else { # conventional mapper needs bin coloring
    num_bins = length(binclust_data)
    colfunc <- colorRampPalette(c("blue", "purple", "red"))
    bin_colors = colfunc(num_bins)
    nodeBorderColors <- mapVisualProperty('node border color', 'bin', 'd', 1:num_bins, bin_colors)
    createVisualStyle(style.name, defaults, list(nodeSizes, edgeWidth, nodeBorderColors, nodeFillColors))
  }

  setVisualStyle(style.name)
}
