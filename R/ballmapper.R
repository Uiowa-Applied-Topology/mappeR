# takes the output of the previous function and makes it suitable for the 1D mapper function
convert_balls <- function(balled_data) {
  ball_sizes = lapply(balled_data, length)
  ballball_data = unlist(mapply(function(x, y) rep(x, y), 1:length(ball_sizes), ball_sizes))
  names(ballball_data) = unlist(balled_data)

  return(ballball_data)
}

construct_ballmappergraph <- function(binclust_data, dists) {
  num_vertices = max(binclust_data[[length(binclust_data)]])

  node_ids = as.character(1:num_vertices)

  overlaps = get_overlaps(binclust_data)
  edges = get_edgelist_from_overlaps(overlaps, num_vertices)
  sources = as.character(edges[,1])
  targets = as.character(edges[,2])

  cluster_tightness = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_size = get_cluster_sizes(binclust_data, num_vertices)
  data_in_cluster = unlist(get_clustered_data(binclust_data, num_vertices))
  edge_weights = get_edge_weights(sapply(overlaps, length), cluster_size, edges)

  nodes = data.frame(id=node_ids,
                     size=cluster_size,
                     tightness=cluster_tightness,
                     data=data_in_cluster)

  edges = data.frame(source=sources,
                     target=targets,
                     weight=edge_weights)


  return(list(nodes, edges))
}

#' Runs ballmapper and returns two dataframes containing node and edge data.
#'
#' @param data Your input data. Ideally a dataframe.
#' @param dists A distance matrix for your data. Can be a `dist` object or 2D matrix.
#' @param eps A positive real number for your desired ball radius.
#' @returns A list of two dataframes, one with node data containing ball size,
#'  datapoints per ball, ball tightness, and one with edge data
#'  containing sources, targets, and weights representing overlap strength.
#' @export
get_ballmapper_data <- function(data, dists, eps) {
  balled_data = create_balls(data, dists, eps)
  formatted_balled_data = convert_balls(balled_data)

  # construct ballmapper graph
  ballmappergraph = construct_ballmappergraph(formatted_balled_data, dists)

  return(ballmappergraph)
}

#' Runs ballmapper and passes data to Cytoscape for visualization.
#'
#' @param data Your input data. Ideally a dataframe.
#' @param dists A distance matrix for your data. Can be a `dist` object or 2D matrix.
#' @param eps A positive real number for your desired ball radius.
#' @returns NULL
#' @export
cyballmapper <- function(data, dists, eps) {
  visualize_mapper_data(get_ballmapper_data(data, dists, eps))
  return(invisible(NULL))
}
