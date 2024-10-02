
# mapper mapper -----------------------------------------------------------

#' Create a mapper object
#'
#' Run the mapper algorithm on a data frame.
#'
#' @param data A data frame.
#' @param dists A distance matrix for the data frame.
#' @param filtered_data The result of a function applied to the data frame; there should be one row per observation in the original data frame.
#' @param cover_element_tests A list of membership test functions for a list of cover elements. Each member of `cover_element_tests` should be able to identify (return `TRUE` or `FALSE`) if a single input data point is a member of the cover element it represents.
#' @param method The desired clustering method to use. e.g., "single"
#'
#' @return A list of two dataframes, one with node data containing bin membership,
#'  datapoints per cluster, and cluster dispersion, and one with edge data
#'  containing sources, targets, and weights representing overlap strength.
#' @export
#' @examples
#' \dontrun{
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#' projx = data$x
#'
#' num_bins = 10
#' percent_overlap = 25
#' xcover = create_width_balanced_cover(min(projx), max(projx), num_bins, percent_overlap)
#'
#' check_in_interval <- function(endpoints) {
#'  return(function(x) (endpoints[1] - x <= 0) & (endpoints[2] - x >= 0))
#' }
#'
#' # each of the "cover" elements will really be a function that checks if a data point lives in it
#' xcovercheck = apply(xcover, 1, check_in_interval)
#'
#' # build the mapper object
#' xmapper = create_mapper_object(
#'   data = data,
#'   dists = dist(data),
#'   filtered_data = projx,
#'   cover_element_tests = xcovercheck,
#'   method = "single"
#' )
#' }
create_mapper_object <- function(data, dists, filtered_data, cover_element_tests, method="none") {
  bins = create_bins(data, filtered_data, cover_element_tests)

  if (method == "none") {
    return(run_mapper(convert_to_clusters(bins), dists, binning = FALSE))
  } else {
    clusters = get_clusters(bins, dists, method)
    return(run_mapper(clusters, dists, binning = TRUE))
  }
}

#' Create a bin of data
#'
#' @param data A data frame.
#' @param filtered_data The result of a function applied to the data frame; there should be one row per observation in the original data frame.
#' @param cover_element_test A membership test function for a cover element. It should identify (return `TRUE` or `FALSE`) if a single input data point, is a member of the cover element it represents.
#'
#' @return A vector of names of points from the data frame, representing a bin.
create_single_bin <- function(data, filtered_data, cover_element_test) {
  in_bin = sapply(filtered_data, cover_element_test)
  bin_assignments = which(in_bin)
  if (length(bin_assignments) != 0) {
    return(rownames(data[bin_assignments, ])) # TODO: bother me about why I need the original dataset here, I think it's more safe but who knows!
  } else {
    return(vector()) # bin still exists, it's just empty
  }
}

#' Create bins of data
#'
#' @param data A data frame.
#' @param filtered_data The result of a function applied to the data frame; there should be one row per observation in the original data frame.
#' @param cover_element_tests A list of membership test functions for a list of cover elements. Each member of `cover_element_tests` should be able to identify (return `TRUE` or `FALSE`) if a single input data point is a member of the cover element it represents.
#'
#' @return A list of bins, each containing a vector of the names of the data inside it.
create_bins <- function(data, filtered_data, cover_element_tests) {
  return(mapply(create_single_bin, cover_element_test = cover_element_tests, MoreArgs = list(data = data, filtered_data = filtered_data)))
}

#' Construct mapper graph from data
#'
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster ids
#' @param dists A distance matrix for the data that has been binned and clustered.
#' @param binning Whether the output dataframe should sort vertices into "bins" or not. Should be true if using clustering, leave false otherwise
#'
#' @return A list of two dataframes, one with node data containing bin membership,
#'  datapoints per cluster, and cluster dispersion, and one with edge data
#'  containing sources, targets, and weights representing overlap strength.
run_mapper <- function(binclust_data, dists, binning=TRUE) {
  num_vertices = max(binclust_data[[length(binclust_data)]])
  node_ids = as.character(1:num_vertices)
  overlaps = get_overlaps(binclust_data)
  edgelist = get_edgelist_from_overlaps(overlaps, num_vertices)
  sources = as.character(edgelist[, 1])
  targets = as.character(edgelist[, 2])

  # calculate some cluster stats
  cluster_tightness = get_cluster_tightness_vector(as.matrix(dists), binclust_data)
  cluster_size = get_cluster_sizes(binclust_data)
  data_in_cluster = unlist(get_clustered_data(binclust_data))
  edge_weights = get_edge_weights(sapply(overlaps, length), cluster_size, edgelist)

  # if you care about bins
  if (binning) {
    nodes = data.frame(
      id = node_ids,
      cluster_size = cluster_size,
      tightness = cluster_tightness,
      data = data_in_cluster,
      bin = get_bin_vector(binclust_data)
    )

    edges = data.frame(source = sources,
                       target = targets,
                       weight = edge_weights)

    return(list(nodes, edges))

  # if you don't
  } else {
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
}


# 1D mapper ---------------------------------------------------------------
#
# a flavor of mapper based on projection to a single coordinate

#' Run 1D mapper
#'
#' Run mapper using a one-dimensional filter, a cover of intervals, and a clustering algorithm.
#'
#' @param data A data frame.
#' @param dists A distance matrix for the data frame.
#' @param filtered_data The result of a function applied to the data frame; there should be one row per observation in the original data frame.
#' @param cover A 2D array of interval left and right endpoints.
#' @param clustering_method Your favorite clustering algorithm.
#'
#' @return A list of two data frames, one with node data containing bin membership,
#'  data points per cluster, and cluster dispersion, and one with edge data
#'  containing sources, targets, and weights representing overlap strength.
#' @export
#' @examples
#' \dontrun{
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#' projx = data$x
#'
#' num_bins = 10
#' percent_overlap = 25
#'
#' cover = get_width_balanced_cover(min(projx), max(projx), num_bins, percent_overlap)
#'
#' create_1D_mapper_object(data, dist(data), projx, cover, "single")
#' }
create_1D_mapper_object <- function(data, dists, filtered_data, cover, clustering_method="single") {

  cover = apply(cover, 1, check_in_interval)

  return(create_mapper_object(data, dists, filtered_data, cover, clustering_method))
}

# ball mapper --------------------------------------------------------------
#
# a flavor of mapper all about the balls

#' Run mapper using a trivial filter, a cover of balls, and no clustering algorithm.
#'
#' Run mapper using an \eqn{\varepsilon}-net cover (greedily generated) and the 2D inclusion function as a filter.
#'
#' @param data A data frame.
#' @param dists A distance matrix for the data frame.
#' @param eps A positive real number for your desired ball radius.
#' @returns A list of two data frames, one with node data containing ball size,
#'  data points per ball, ball tightness, and one with edge data
#'  containing sources, targets, and weights representing overlap strength.
#' @export
#' @examples
#' \dontrun{
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#' eps = 1
#'
#' create_ball_mapper_object(data, dist(data), eps)
#' }
create_ball_mapper_object <- function(data, dists, eps) {
  balled_data = create_balls(data, dists, eps)

  ball_mapper_object = run_mapper(convert_to_clusters(balled_data), dists, binning = FALSE)

  return(ball_mapper_object)
}


# clusterball mapper ------------------------------------------------------
#
# a flavor of mapper that's just clustering in the balls of ball mapper

#' Run clusterball mapper
#'
#' Run ball mapper, but additionally cluster within the balls. Can use two different distance matrices to accomplish this.
#'
#' @param data A data frame.
#' @param dist1 A distance matrix for the data frame; this will be used to ball the data.
#' @param dist2 Another distance matrix for the data frame; this will be used to cluster the data after balling.
#' @param eps A positive real number for your desired ball radius.
#' @param clustering_method Your favorite clustering algorithm.
#'
#' @return A list of two dataframes, one with node data containing bin membership,
#'  datapoints per cluster, and cluster dispersion, and one with edge data
#'  containing sources, targets, and weights representing overlap strength.
#' @export
#' @examples
#' \dontrun{
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#' data.dists = dist(data)
#' eps = 1
#'
#' create_clusterball_mapper_object(data, data.dists, data.dists, eps, "single")
#' }
create_clusterball_mapper_object <- function(data, dist1, dist2, eps, clustering_method) {
  balls = create_balls(data, dist1, eps)

  is_in_ball <- function(ball) {
    return(function(x) x %in% ball)
  }

  return(create_mapper_object(data, dist2, rownames(data), lapply(balls, is_in_ball), clustering_method))
}

# graph construction ------------------------------------------------------

#' Find which triangular number you're on
#'
#' @param x A positive integer.
#'
#' @return The index of the next greatest or equal triangular number to \eqn{x}.
next_triangular <- function(x) {
  next_triangle_indx = floor((1 + sqrt(1 + 8 * x)) / 2)
  prev_triangle_val = choose(next_triangle_indx, 2)
  if (prev_triangle_val == x) {
    return (next_triangle_indx - 1)
  } else {
    return (next_triangle_indx)
  }
}

#' Get cluster overlaps
#'
#' @param binclust_data A list of bins, each containing named vectors whose names are those of data points and whose values are cluster ids.
#'
#' @return A named list of edges, whose elements contain the names of clusters in the overlap represented by that edge.
get_overlaps <- function(binclust_data) {
  num_vertices = max(binclust_data[[length(binclust_data)]]) # id of last cluster in the last bin
  if (num_vertices == 1) {
    return(0)
  }
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

#' Obtain edge list from cluster intersections
#'
#' @param overlaps A named list of edges, whose elements contain the names of clusters in the overlap represented by that edge; output of [get_overlaps()].
#' @param num_vertices The number of vertices in the graph.
#'
#' @return A 2D array representing the edge list of a graph.
get_edgelist_from_overlaps <- function(overlaps, num_vertices) {
  if (num_vertices == 2) {
    return(matrix(c(1,2), nrow = 1, ncol = 2))
  } else {
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
}

