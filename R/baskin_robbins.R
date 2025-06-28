###########################################################################
# BASKIN-ROBBINS
# Mapper flavors
###########################################################################


# 1D Mapper ---------------------------------------------------------------
#
# a flavor of Mapper based on projection to a single coordinate

#' 1D Mapper
#'
#' Run Mapper using a one-dimensional filter, a cover of intervals, and a clusterer.
#'
#' @param data A data frame.
#' @param dists A distance matrix for the data frame.
#' @param filtered_data The result of a function applied to the data frame; there should be one filter value per observation in the original data frame.
#' @param cover An \eqn{n \times 2} matrix of interval left and right endpoints; rows should be intervals and columns left and right endpoints (in that order).
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix.
#'
#' @return A list of two data frames, `nodes` and `edges`, which contain information about the Mapper graph constructed from the given parameters.
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
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#' projx = data$x
#'
#' num_bins = 10
#' percent_overlap = 25
#'
#' cover = create_width_balanced_cover(min(projx), max(projx), num_bins, percent_overlap)
#'
#' create_1D_mapper_object(data, dist(data), projx, cover)
create_1D_mapper_object <- function(data,
                                    dists,
                                    filtered_data,
                                    cover,
                                    clusterer = global_hierarchical_clusterer("single", dists)) {
  if (!all(cover[, 1] - cover[, 2] <= 0)) {
    stop("left endpoints must be less than or equal to right endpoints")
  }

  cover = apply(cover, 1, check_in_interval)

  return(create_mapper_object(data, dists, filtered_data, cover, clusterer = clusterer))
}

# ball Mapper --------------------------------------------------------------
#
# a flavor of Mapper all about the balls

#' "Clustering" for ballmapper just means treating each bin as its own cluster.
#'
#' @param bins A list of bins, each containing names of data from some data frame.
#'
#' @return A named vector whose names are data point names and whose values are cluster labels
convert_to_clusters <- function(bins) {
  ball_sizes = lapply(bins, length)

  # repeat the cluster id for as many data points belonging to that bin
  ballball_data = unlist(mapply(function(x, y)
    rep(x, y), 1:length(ball_sizes), ball_sizes))

  # make sure names match up
  names(ballball_data) = unlist(bins)

  return(ballball_data)
}

#' Run Mapper using a trivial filter, a cover of balls, and no clustering algorithm.
#'
#' Run Mapper using an \eqn{\varepsilon}-net cover (greedily generated) and the 2D inclusion function as a filter.
#'
#' @param data A data frame.
#' @param dists A distance matrix for the data frame.
#' @param eps A positive real number for your desired ball radius.
#' @return A list of two data frames, `nodes` and `edges`, which contain information about the Mapper graph constructed from the given parameters.
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
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#' eps = .5
#'
#' create_ball_mapper_object(data, dist(data), eps)
create_ball_mapper_object <- function(data, dists, eps) {
  if (!is.data.frame(data)) {
    stop("input data needs to be a data frame.")
  } else if (!is.numeric(eps)) {
    stop("epsilon needs to be a number")
  } else if (eps <= 0) {
    stop("epsilon needs to be positive")
  }

  if (any(is.na(dists))) {
    stop("no distance value can be NA")
  }

  balled_data = create_balls(data, dists, eps)

  ball_mapper_object = assemble_mapper_object(convert_to_clusters(balled_data), dists, binning = FALSE)

  return(ball_mapper_object)
}


# clusterball Mapper ------------------------------------------------------
#
# a flavor of Mapper that's just clustering in the balls of ball Mapper

#' Run clusterball Mapper
#'
#' Run ball Mapper, but additionally cluster within the balls. Can use two different distance matrices to accomplish this.
#'
#' @param data A data frame.
#' @param dist1 A distance matrix for the data frame; this will be used to ball the data.
#' @param dist2 Another distance matrix for the data frame; this will be used to cluster the data after balling.
#' @param eps A positive real number for your desired ball radius.
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix.
#'
#' @return A list of two data frames, `nodes` and `edges`, which contain information about the Mapper graph constructed from the given parameters.
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
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#' data.dists = dist(data)
#' eps = 1
#'
#' create_clusterball_mapper_object(data, data.dists, data.dists, eps)
create_clusterball_mapper_object <- function(data, dist1, dist2, eps, clusterer = local_hierarchical_clusterer("single")) {
  if (!is.data.frame(data)) {
    stop("input data needs to be a data frame.")
  } else if (!is.numeric(eps)) {
    stop("epsilon needs to be a number")
  } else if (eps <= 0) {
    stop("epsilon needs to be positive")
  }

  if ((any(is.na(dist1))) | (any(is.na(dist2)))) {
    stop("no distance value can be NA")
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
