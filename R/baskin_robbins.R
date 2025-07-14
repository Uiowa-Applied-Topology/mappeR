###########################################################################
# BASKIN-ROBBINS
# Mapper flavors
###########################################################################


# 1D Mapper ---------------------------------------------------------------
#
# a flavor of Mapper based on projection to a single coordinate

#' One-Dimensional Mapper
#'
#' Run Mapper using a one-dimensional filter, a cover of the codomain of intervals, and a clusterer.
#'
#' @param data A data frame.
#' @param dists A distance matrix associated to the data frame. Can be a `dist` object or `matrix`.
#' @param filtered_data The result of a function applied to the data frame; there should be one filter value per observation in the original data frame.
#' These values need to be named, and the names of these values must match the names of the original data set.
#' @param cover An \eqn{n \times 2} `matrix` of interval left and right endpoints; rows should be intervals and columns left and right endpoints (in that order).
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix;
#' that is, it should return a list of named vectors, whose name are the names of data points and whose values are cluster assignments (integers).
#' If this value is omitted, then trivial clustering will be done.
#'
#' @return A `list` of two data frames, `nodes` and `edges`, which contain information about the Mapper graph constructed from the given parameters.
#'
#' The node data frame consists of:
#'
#' - `id`: vertex ID
#' - `cluster_size`: number of data points in vertex
#' - `mean_dist_to_medoid`: mean distance to medoid of vertex
#' - `data`: names of data points in cluster
#' - `patch`: level set ID
#'
#' The `edge` data frame contains consists of:
#'
#' - `source`: vertex ID of edge source
#' - `target`: vertex ID of edge target
#' - `weight`: Jaccard index of edge; this is the size of the intersection between the vertices divided by the union
#' - `overlap_data`: names of data points in overlap
#' - `overlap_size`: number of data points overlap
#'
#' @export
#' @examples
#' # Create noisy circle data
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#'
#' # Project to horizontal axis as lens
#' projx = data$x
#' names(projx) = row.names(data)
#'
#' # Create a one-dimensional cover
#' num_bins = 10
#' percent_overlap = 25
#' cover = create_width_balanced_cover(min(projx), max(projx), num_bins, percent_overlap)
#'
#' # Build Mapper object
#' create_1D_mapper_object(data, dist(data), projx, cover)
create_1D_mapper_object <- function(data,
                                    dists,
                                    filtered_data,
                                    cover,
                                    clusterer = global_hierarchical_clusterer("single", dists)) {
  if (!all(cover[, 1] - cover[, 2] <= 0)) {
    stop("Left endpoints in the cover must be less than or equal to right endpoints.")
  }

  cover = apply(cover, 1, check_in_interval)

  return(create_mapper_object(data, dists, filtered_data, cover, clusterer = clusterer))
}

# Ball Mapper --------------------------------------------------------------
#
# a flavor of Mapper all about the balls

#' Convert Balls to Clusters
#'
#' Perform trivial clustering on a set of balled data.
#'
#' @param bins A `list` of bins, each containing a named vector of data points.
#'
#' @return A named vector whose names are data point names and whose values are cluster labels (`integer`s).
#' @noRd
convert_to_clusters <- function(bins) {
  ball_sizes = lapply(bins, length)

  # repeat the cluster id for as many data points belonging to that bin
  ballball_data = unlist(mapply(function(x, y)
    rep(x, y), 1:length(ball_sizes), ball_sizes))

  # make sure names match up
  names(ballball_data) = unlist(bins)

  return(ballball_data)
}

#' Ball Mapper
#'
#' Run Mapper using the identity function as a lens and an \eqn{\varepsilon}-net cover, greedily generated using a distance matrix.
#'
#' @param data A data frame.
#' @param dists A distance matrix for the data frame. Can be a `dist` object or a `matrix`.
#' @param eps A positive real number for the desired ball radius.
#' @return A `list` of two data frames, `nodes` and `edges`, which contain information about the Mapper graph constructed from the given parameters.
#'
#' The node data frame consists of:
#'
#' - `id`: vertex ID
#' - `cluster_size`: number of data points in vertex
#' - `mean_dist_to_medoid`: mean distance to medoid of vertex
#' - `data`: names of data points in cluster
#'
#' The `edge` data frame contains consists of:
#'
#' - `source`: vertex ID of edge source
#' - `target`: vertex ID of edge target
#' - `weight`: Jaccard index of edge; this is the size of the intersection between the vertices divided by the union
#' - `overlap_data`: names of data points in overlap
#' - `overlap_size`: number of data points overlap
#'
#' @export
#' @examples
#' # Create noisy cirle data set
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#'
#' # Set ball radius
#' eps = .5
#'
#' # Create Mapper object
#' create_ball_mapper_object(data, dist(data), eps)
create_ball_mapper_object <- function(data, dists, eps) {
  if (!is.data.frame(data)) {
    stop("Input data needs to be a data frame.")
  } else if (any(is.na(data))) {
    stop("Data cannot have NA values.")
  } else if (any(is.na(dists))) {
    stop("No distance value can be NA.")
  } else if (!is.numeric(eps)) {
    stop("Epsilon parameter needs to be numeric.")
  } else if (eps <= 0) {
    stop("Epsilon parameter needs to be positive.")
  }

  if (nrow(data) != dim(as.matrix(dists))[1]) {
    stop("Your distance matrix dimensions are not correct for your data.")
  } else if (dim(as.matrix(dists))[1] != dim(as.matrix(dists))[2]) {
    stop("Your distance matrix is not square!")
  } else if (any(!is.numeric(dists))) {
    stop("Your distance matrix has non-numeric entries!")
  }

  if (length(data) == 0) {
    stop("Your data is missing!")
  } else if (length(dists) == 0) {
    stop("Your distance matrix is missing!")
  }

  if (any(row.names(as.matrix(dists)) != row.names(data))) {
    stop("Names of points in distance matrix need to match names in data frame!")
  }

  balled_data = create_balls(data, dists, eps)

  ball_mapper_object = assemble_mapper_object(convert_to_clusters(balled_data), dists, binning = FALSE)

  return(ball_mapper_object)
}


# clusterball Mapper ------------------------------------------------------
#
# a flavor of Mapper that's just clustering in the balls of ball Mapper

#' ClusterBall Mapper
#'
#' Run Ball Mapper, but non-trivially cluster within the balls. You can use two different distance matrices to for the balling and clustering.
#'
#' @param data A data frame.
#' @param dist1 A distance matrix for the data frame; this will be used to ball the data. It can be a `dist` object or a `matrix`.
#' @param dist2 Another distance matrix for the data frame; this will be used to cluster the data after balling. It can be a `dist` object or a `matrix`.
#' @param eps A positive real number for the desired ball radius.
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix;
#' that is, it should return a list of named vectors, whose name are the names of data points and whose values are cluster assignments (integers).
#' If this value is omitted, then single-linkage clustering will be done (and cutting heights will be decided for you).
#' @return A `list` of two data frames, `nodes` and `edges`, which contain information about the Mapper graph constructed from the given parameters.
#'
#' The node data frame consists of:
#'
#' - `id`: vertex ID
#' - `cluster_size`: number of data points in vertex
#' - `mean_dist_to_medoid`: mean distance to medoid of vertex
#' - `data`: names of data points in cluster
#' - `patch`: level set ID
#'
#' The `edge` data frame contains consists of:
#'
#' - `source`: vertex ID of edge source
#' - `target`: vertex ID of edge target
#' - `weight`: Jaccard index of edge; this is the size of the intersection between the vertices divided by the union
#' - `overlap_data`: names of data points in overlap
#' - `overlap_size`: number of data points overlap
#'
#' @export
#' @examples
#' # Create noisy circle data set
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#' data.dists = dist(data)
#'
#' # Set ball radius
#' eps = 1
#'
#' # Do single-linkage clustering in the balls to produce Mapper graph
#' create_clusterball_mapper_object(data, data.dists, data.dists, eps)
create_clusterball_mapper_object <- function(data, dist1, dist2, eps, clusterer = local_hierarchical_clusterer("single")) {
  if (!is.data.frame(data)) {
    stop("Input data needs to be a data frame.")
  } else if (any(is.na(data))) {
    stop("Data cannot have NA values.")
  } else if (any(is.na(dist1)) | any(is.na(dist2))) {
    stop("No distance value can be NA.")
  } else if (!is.numeric(eps)) {
    stop("Epsilon parameter needs to be numeric.")
  } else if (eps <= 0) {
    stop("Epsilon parameter needs to be positive.")
  }

  if (nrow(data) != dim(as.matrix(dists))[1]) {
    stop("Your distance matrix dimensions are not correct for your data.")
  } else if (dim(as.matrix(dists))[1] != dim(as.matrix(dists))[2]) {
    stop("Your distance matrix is not square!")
  } else if (any(!is.numeric(dists))) {
    stop("Your distance matrix has non-numeric entries!")
  }

  if (length(data) == 0) {
    stop("Your data is missing!")
  } else if (length(dist1) == 0 | length(dist2) == 0) {
    stop("Your distance matrix is missing!")
  }

  if (any(row.names(as.matrix(dist1)) != row.names(data)) | any(row.names(as.matrix(dist2)) != row.names(data))) {
    stop("Names of points in distance matrices need to match names in data frame!")
  }

  balls = create_balls(data, dist1, eps)

  projection = row.names(data)
  names(projection) = row.names(data) # label everything just trust me ok

  return(create_mapper_object(
    data,
    dist2,
    projection,
    lapply(balls, is_in_ball),
    clusterer = clusterer
  ))
}
