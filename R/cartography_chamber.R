
# 1D mapper ---------------------------------------------------------------

#' Run 1D mapper algorithm
#'
#' @param binclust_data A list containing named vectors whose names are datapoint names and whose values are cluster labels.
#' @param dists A distance matrix containing pairwise distances between named datapoints.
#'
#' @return A list containing:
#' \itemize{
#'    \item nodes - A dataframe containing vertex/cluster information:
#'      \itemize{
#'          \item id - The vertex label.
#'          \item size - How many datapoints are contained in that vertex.
#'          \item tightness - A measure of dispersion of the data inside the vertex.
#'          \item data - A string containing names of the datapoints within the vertex.
#'          \item bin - Which bin the vertex belongs to.
#'      }
#'    \item edges - A dataframe containing vertex/cluster overlap information:
#'      \itemize{
#'          \item source - One endpoint of an edge.
#'          \item target - The other endpoint of an edge.
#'          \item weight - A measure of the significance of the edge.
#'      }
#' }
construct_1Dmappergraph <- function(binclust_data, dists) {
  # grab basic graph characteristics
  num_vertices = max(binclust_data[[length(binclust_data)]])
  node_ids = as.character(1:num_vertices)
  overlaps = get_overlaps(binclust_data)
  edges = get_edgelist_from_overlaps(overlaps, num_vertices)
  sources = as.character(edges[, 1])
  targets = as.character(edges[, 2])

  # calculate some cluster stats
  cluster_tightness = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_size = get_cluster_sizes(binclust_data, num_vertices)
  data_in_cluster = unlist(get_clustered_data(binclust_data, num_vertices))
  edge_weights = get_edge_weights(sapply(overlaps, length), cluster_size, edges)
  bins = get_bin_vector(binclust_data)

  # assemble graph data
  nodes = data.frame(
    id = node_ids,
    cluster_size = cluster_size,
    tightness = cluster_tightness,
    data = data_in_cluster,
    bin = bins
  )

  edges = data.frame(source = sources,
                     target = targets,
                     weight = edge_weights)


  return(list(nodes, edges))
}

#' Runs 1D mapper returns a dataframe with node and edge data.
#'
#' @param data A dataframe.
#' @param filtered_data A single column of the input data.
#' @param dists A distance matrix containing pairwise relations among the input data. Can be a `dist` object or 2D matrix.
#' @param num_bins The number of "bins" to split the input data into based on the filter. A positive integer.
#' @param percent_overlap The percent overlap desired between each "bin." An integer between 0 and 100 (inclusive).
#' @param clustering_method Desired clustering method. A string from these options: "single" (single-linkage hierarchical)
#'
#' @return A list containing:
#' \itemize{
#'    \item nodes - A dataframe containing vertex/cluster information:
#'      \itemize{
#'          \item id - The vertex label.
#'          \item size - How many datapoints are contained in that vertex.
#'          \item tightness - A measure of dispersion of the data inside the vertex.
#'          \item data - A string containing names of the datapoints within the vertex.
#'          \item bin - Which bin the vertex belongs to.
#'      }
#'    \item edges - A dataframe containing vertex/cluster overlap information:
#'      \itemize{
#'          \item source - One endpoint of an edge.
#'          \item target - The other endpoint of an edge.
#'          \item weight - A measure of the significance of the edge.
#'      }
#' }
#' @export
get_1D_mapper_data <- function(data,
                               filtered_data,
                               dists,
                               num_bins,
                               percent_overlap,
                               clustering_method) {
  # bin data according to filter values
  bins = create_width_balanced_cover(min(filtered_data),
                                     max(filtered_data),
                                     num_bins,
                                     percent_overlap)
  binned_data = make_bins(data, filtered_data, bins)

  # cluster data
  binclust_data = get_clusters(binned_data, dists, clustering_method)

  # construct mapper graph
  mappergraph = construct_1Dmappergraph(binclust_data, dists)

  return(mappergraph)
}

#' Runs 1D mapper and passes data to Cytoscape for visualization.
#'
#' @param data A dataframe.
#' @param filtered_data A single column of the input data.
#' @param dists A distance matrix containing pairwise relations among the input data. Can be a `dist` object or 2D matrix.
#' @param num_bins The number of "bins" to split the input data into based on the filter. A positive integer.
#' @param percent_overlap The percent overlap desired between each "bin." An integer between 0 and 100 (inclusive).
#' @param clustering_method Desired clustering method. A string from these options: "single" (single-linkage hierarchical)
#'
#' @returns NULL
#'
#' @export
cymapper <- function(mapperobject) {

  # pass to visualizer for........visualizing...
  visualize_mapper_data(mapperobject, is_ballmapper = FALSE)

  # if this isn't here R will print something useless
  return(invisible(NULL))
}


# ballmapper --------------------------------------------------------------

# TODO: add rox docs and all that jazz

construct_ballmappergraph <- function(binclust_data, dists) {
  num_vertices = max(binclust_data[[length(binclust_data)]])

  node_ids = as.character(1:num_vertices)

  overlaps = get_overlaps(binclust_data)
  edges = get_edgelist_from_overlaps(overlaps, num_vertices)
  sources = as.character(edges[, 1])
  targets = as.character(edges[, 2])

  cluster_tightness = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_size = get_cluster_sizes(binclust_data, num_vertices)
  data_in_cluster = unlist(get_clustered_data(binclust_data, num_vertices))
  edge_weights = get_edge_weights(sapply(overlaps, length), cluster_size, edges)

  nodes = data.frame(
    id = node_ids,
    cluster_size = cluster_size,
    tightness = cluster_tightness,
    data = data_in_cluster
  )

  edges = data.frame(source = sources,
                     target = targets,
                     weight = edge_weights)


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
cyballmapper <- function(mapperobject) {
  visualize_mapper_data(mapperobject, FALSE)
  return(invisible(NULL))
}


# clusterball mapper ------------------------------------------------------

# TODO: test that this works lol

# runner function for combo mapper; outputs bins, clusters, and the mapper graph.
get_clusterballmapper_data <- function(data, dist1, dist2, eps, method) {
  balls = create_balls(data, dist1, eps)
  clusters = get_clusters(balls, dist2, method)

  clusterballmappergraph = construct_1Dmappergraph(clusters, dist2)

  return(clusterballmappergraph)
}

cyclusterballmapper <- function(mapperobject) {
  visualize_mapper_data(mapperobject, FALSE)

  return(invisible(NULL))
}

# graph construction ------------------------------------------------------

# TODO: add rox docs

get_overlaps <- function(binclust_data) {
  num_vertices = max(binclust_data[[length(binclust_data)]]) # id of last cluster in the last bin

  flattened_data = unlist(binclust_data)
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x]) # sort by cluster
  cluster_names = lapply(clusters, names) # it doesn't work if you don't do this

  pairs = combn(cluster_names, 2) # get all pairs of clusters
  raw_overlaps = apply(pairs, 2, function(x)
    intersect(x[[1]], x[[2]])) # get all intersections between clusters
  names(raw_overlaps) = 1:length(raw_overlaps)
  overlaps = Filter(length, raw_overlaps) # filter out the empty intersections

  return(overlaps)
}

next_triangular <- function(x) {
  next_triangle_indx = floor((1 + sqrt(1 + 8 * x)) / 2)
  prev_triangle_val = choose(next_triangle_indx, 2)
  if (prev_triangle_val == x) {
    return (next_triangle_indx - 1)
  } else {
    return (next_triangle_indx)
  }
}

get_edgelist_from_overlaps <- function(overlaps, num_vertices) {
  overlap_names = rev(-as.numeric(names(overlaps)) + choose(num_vertices, 2) + 1)
  sources = sapply(overlap_names, function(x)
    num_vertices - next_triangular(x))
  targets = sapply(overlap_names, function(x) {
    k = next_triangular(x)
    diff = k * (k + 1) / 2 - x
    num_vertices - k + diff + 1
  })
  edges = cbind(rev(sources), rev(targets))
  return(edges)
}

