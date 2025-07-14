###########################################################################
# CARTOGRAPHY CHAMBER
# mapper
###########################################################################

# mapper mapper -----------------------------------------------------------

#' Mapper
#'
#' Run the Mapper algorithm on a data frame.
#'
#' @param data A data frame.
#' @param dists A distance matrix for the data frame. Can be a `dist` object or `matrix`.
#' @param filtered_data The result of a function applied to the data frame; there should be one filter value per observation in the original data frame.
#' These values need to be named, and the names of these values must match the names of the original data set.
#' @param cover_element_tests A list of membership test functions for a set of cover elements. In other words, each element of `cover_element_tests` is a function that returns `TRUE` or `FALSE` when given a filter value.
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix;
#' that is, it should return a list of named vectors, whose names are the names of data points and whose values are cluster assignments (integers).
#' If this value is omitted, then trivial clustering will be done.
#' @return A `list` of two data frames, one with node data and one with edge data. The node data includes:
#'
#' - `id`: vertex ID
#' - `cluster_size`: number of data points in cluster
#' - `mean_dist_to_medoid`: mean distance to medoid of cluster
#' - `data`: names of data points in cluster
#' - `patch`: level set ID
#'
#' The edge data includes:
#'
#' - `source`: vertex ID of edge source
#' - `target`: vertex ID of edge target
#' - `weight`: Jaccard index of edge; intersection divided by union
#' - `overlap_data`: names of data points in overlap
#' - `overlap_size`: number of data points overlap
#'
#' @export
#' @examples
#' # Create noisy data around a circle
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#'
#' # Apply lens function to data
#' projx = data$x
#' names(projx) = row.names(data)
#'
#' # Build a width-balanced cover with 10 intervals and 25 percent overlap
#' num_bins = 10
#' percent_overlap = 25
#' xcover = create_width_balanced_cover(min(projx), max(projx), num_bins, percent_overlap)
#'
#' # Write a function which can check if a data point lives in an interval of our lens function
#' check_in_interval <- function(endpoints) {
#'  return(function(x) (endpoints[1] - x <= 0) & (endpoints[2] - x >= 0))
#' }
#'
#' # Each of the "cover" elements will really be a function that checks if a data point lives in it
#' xcovercheck = apply(xcover, 1, check_in_interval)
#'
#' # Build the mapper object
#' xmapper = create_mapper_object(
#'   data = data,
#'   dists = dist(data),
#'   filtered_data = projx,
#'   cover_element_tests = xcovercheck
#' )
create_mapper_object <- function(data,
                                 dists,
                                 filtered_data,
                                 cover_element_tests,
                                 clusterer = NULL) {
  if (!is.data.frame(data)) {
    stop("Input data needs to be a data frame.")
  } else if (!is.list(cover_element_tests)) {
    stop("Cover element tests need to be in a list.")
  } else if (!all(sapply(cover_element_tests, typeof) == "closure")) {
    stop("Cover element tests need to be boolean functions.")
  } else if (any(is.na(data))) {
    stop("Data cannot have NA values.")
  } else if (any(is.na(dists))) {
    stop("No distance value can be NA.")
  } else if (any(is.na(filtered_data))) {
    stop("Filtered data cannot have NA values.")
  } else if (any(is.na(cover_element_tests))) {
    stop("Cover element functions cannot be NA!")
  }

  if (length(data) == 0) {
    stop("Your data is missing!")
  } else if (length(dists) == 0) {
     stop("Your distance matrix is missing!")
  } else if (length(filtered_data) == 0) {
    stop("Your lens/filter is missing!")
  } else if (length(cover_element_tests) == 0) {
    stop("Your cover is missing!")
  }

  if (any(row.names(as.matrix(dists)) != row.names(data))) {
    stop("Names of points in distance matrix need to match names in data frame!")
  }


  if ((is.matrix(filtered_data))) {
    if (dim(filtered_data)[1] != nrow(data)) {
      stop("There should be as many filtered data points as there are data points.")
    }
    if (is.null(row.names(filtered_data)) | any(row.names(filtered_data) != row.names(data))) {
      stop("The names of the filtered data points should match the names of the original data points.")
    }
  } else if (is.data.frame(filtered_data)) {
    if (nrow(filtered_data) != nrow(data)) {
      stop("There should be as many filtered data points as there are data points.")
    } else if (is.null(row.names(filtered_data)) | any(row.names(filtered_data) != row.names(data))) {
      stop("The names of the filtered data points should match the names of the original data points.")
    }
  } else {
    if (length(filtered_data) != nrow(data)) {
      stop("There should be as many filtered data points as there are data points.")
    } else if (is.null(names(filtered_data)) | any(names(filtered_data) != row.names(data))) {
      stop("The names of the filtered data points should match the names of the original data points.")
    }
  }

  bins = create_bins(data, filtered_data, cover_element_tests)

  if (is.null(clusterer) | length(clusterer) == 0 | !is.function(clusterer)) {
    return(assemble_mapper_object(convert_to_clusters(bins), dists, binning = FALSE))
  } else {
    clusters = get_clusters(bins, dists, clusterer)
    return(assemble_mapper_object(clusters, dists, binning = TRUE))
  }
}

#' Level Set Maker
#'
#' @param data A data frame.
#' @param filtered_data The result of a function applied to the data frame; there should be one filter value per observation in the original data frame.
#' These values need to be named, and the names of these values must match the names of the original data set.
#' @param cover_element_test A membership test function for a cover element. It should return `TRUE` or `FALSE` when given a filtered data point.
#'
#' @return A vector of names of points from the data frame, representing a level set.
#' @noRd
create_single_bin <- function(data, filtered_data, cover_element_test) {

  # find which data points are part of the cover element
  in_bin = sapply(filtered_data, cover_element_test)
  bin_assignments = which(in_bin)

  # return the level set
  if (length(bin_assignments) != 0) {
    return(rownames(data[bin_assignments, ])) # TODO: bother me about why I need the original data set here, I think it's more safe but who knows!
  } else {
    return(vector()) # bin still exists, it's just empty
  }
}

#' Level Sets Maker
#'
#' @param data A data frame.
#' @param filtered_data The result of a function applied to the data frame; there should be one filter value per observation in the original data frame.
#' These values need to be named, and the names of these values must match the names of the original data set.
#' @param cover_element_tests A list of membership test functions for a set of cover elements. In other words, each element of `cover_element_tests` is a function that returns `TRUE` or `FALSE` when given a filter value.
#'
#' @return A `list` of vectors, where each one contains names of data points for which a specific cover element test was `TRUE`.
#' @noRd
create_bins <- function(data, filtered_data, cover_element_tests) {
  res = mapply(
    create_single_bin,
    cover_element_test = cover_element_tests,
    SIMPLIFY = FALSE,
    MoreArgs = list(data = data, filtered_data = filtered_data)
  )
  if (length(res) == 0) {
    stop("No filtered data is covered!")
  }
  return(res)
}

#' Mapper Object Assembler
#'
#' @param binclust_data A list of lists, each containing named vectors whose names are those of data points and whose values are cluster IDs (integers).
#' @param dists A distance matrix for the data that has been binned and clustered. Can be a `dist` object or a `matrix`.
#' @param binning Whether the output data frame should sort vertices into "patches" or not. Should be true if using clustering, leave `FALSE` otherwise (for, say Ball Mapper).
#'
#' @return A `list` of two data frames, one with node data and one with edge data. The node data includes:
#'
#' - `id`: vertex ID
#' - `cluster_size`: number of data points in cluster
#' - `mean_dist_to_medoid`: mean distance to medoid of cluster
#' - `data`: names of data points in cluster
#' - `patch`: level set ID (if `binning` was `TRUE`)
#'
#' The edge data includes:
#'
#' - `source`: vertex ID of edge source
#' - `target`: vertex ID of edge target
#' - `weight`: Jaccard index of edge; intersection divided by union
#' - `overlap_data`: names of data points in overlap
#' - `overlap_size`: number of data points overlap
#' @noRd
assemble_mapper_object <- function(binclust_data, dists, binning = TRUE) {

  # basic node info
  num_vertices = max(binclust_data[[length(binclust_data)]])
  node_ids = as.character(1:num_vertices)

  # basic edge info
  overlaps = get_overlaps(binclust_data)
  edgelist = get_edgelist_from_overlaps(overlaps, num_vertices)
  sources = as.character(edgelist[, 1])
  targets = as.character(edgelist[, 2])

  # calculate extra statistics
  cluster_tightness = get_cluster_tightness_vector(as.matrix(dists), binclust_data)
  cluster_size = get_cluster_sizes(binclust_data)
  data_in_cluster = unlist(get_clustered_data(binclust_data))
  edge_weights = get_edge_weights(sapply(overlaps, length), cluster_size, edgelist)

  # assemble edge dataframe
  if ((is.list(overlaps)) & length(overlaps) != 0) {
    data_in_overlap = sapply(overlaps, function(x)
      paste(x, collapse = ", "))
    edges = data.frame(
      source = sources,
      target = targets,
      weight = edge_weights,
      overlap_data = data_in_overlap,
      overlap_size = sapply(overlaps, length)
    )
  } else { # oops no edges
    edges = data.frame(source = "", target = "")
  }

  # assemble node data frame, depending on if we care about patch ID
  if (binning) {
    nodes = data.frame(
      id = node_ids,
      cluster_size = cluster_size,
      mean_dist_to_medoid = cluster_tightness,
      data = data_in_cluster,
      patch = get_bin_vector(binclust_data)
    )

    return(list(nodes, edges))

  } else {
    nodes = data.frame(
      id = node_ids,
      cluster_size = cluster_size,
      mean_dist_to_medoid = cluster_tightness,
      data = data_in_cluster
    )

    return(list(nodes, edges))
  }
}

# graph construction ------------------------------------------------------

#' Triangular Number Triangulator
#'
#' Given a number, find the cardinal order of the next closest triangular number.
#'
#' @param x A positive integer.
#'
#' @return The index of the next greatest or equal triangular number to \eqn{x}.
#' @noRd
next_triangular <- function(x) {
  # do a little algebra
  next_triangle_indx = floor((1 + sqrt(1 + 8 * x)) / 2)

  # reverse the algebra
  prev_triangle_val = choose(next_triangle_indx, 2)
  if (prev_triangle_val == x) {
    return(next_triangle_indx - 1)
  } else {
    return(next_triangle_indx)
  }
}

#' Overlap Finder
#'
#' Find single intersections between sets of data points.
#'
#' @param binclust_data A `list` of bins, each containing named vectors whose names are those of data points and whose values are cluster IDs (integers).
#'
#' @return A named `list` of edges, whose elements contain the names of clusters in the overlap represented by that edge.
#' @noRd
get_overlaps <- function(binclust_data) {

  # no need to list by level set
  num_vertices = max(binclust_data[[length(binclust_data)]])
  flattened_data = unlist(binclust_data)

  # find data points with a specific cluster id
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x])
  cluster_names = lapply(clusters, names)

  # we need at least 2 nodes to make an edge
  if (length(cluster_names) < 2) {
    return(0)
  }

  # get all pairs of data points between clusters and find intersections
  pairs = combn(cluster_names, 2, simplify = FALSE)
  raw_overlaps = lapply(pairs, function(x)
    intersect(x[[1]], x[[2]]))

  if (length(raw_overlaps) == 0) {
    return(0)
  } else {
    names(raw_overlaps) = 1:length(raw_overlaps)
    overlaps = Filter(length, raw_overlaps) # filter out the empty intersections
    return(overlaps)
  }
}

#' Edge List Maker
#'
#' Turn the edge labels given by [get_overlaps()] into a usable edge list.
#'
#' @param overlaps A named `list` of edges, whose elements contain the names of clusters in the overlap represented by that edge; i.e., the output of [get_overlaps()].
#' The names should be integers corresponding to their ordering among all possible pairs of vertices.
#' @param num_vertices The number of vertices in the graph.
#'
#' @return A 2D array representing the edge list of a graph.
#' @noRd
get_edgelist_from_overlaps <- function(overlaps, num_vertices) {

  # label all edges in order
  overlap_names = rev(-as.numeric(names(overlaps)) + choose(num_vertices, 2) + 1)

  # create source and target node list from edge labels
  sources = sapply(overlap_names, function(x)
    num_vertices - next_triangular(x))
  targets = sapply(overlap_names, function(x) {
    k = next_triangular(x)
    diff = k * (k + 1) / 2 - x
    num_vertices - k + diff + 1
  })

  # assemble edge list
  edges = cbind(rev(sources), rev(targets))
  return(edges)

}
