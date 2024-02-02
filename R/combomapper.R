source("R/ballmapper.R")
source("R/1Dmapper.R")
source("R/combomapper_viz.R")

# runner function for 1D mapper; outputs bins, clusters, and the mapper graph.
get_combomapper_data <- function(data, dist1, dist2) {
  print("binning...")
  balls = create_balls(data, dist1, eps)

  print("clustering...")
  binclust_data = get_clusters(balls, dist2, "single")

  print("making mapper graph...")
  graph_data = construct_graph(binclust_data)
  amat = graph_data[[1]]
  edge_overlaps = graph_data[[2]]
  mapper_graph = graph_from_adjacency_matrix(amat, mode="max")

  return(list(binclust_data, mapper_graph, edge_overlaps))
}

