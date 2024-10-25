###########################################################################
# BASKIN-ROBBINS
# mapper flavors
###########################################################################


# 1D mapper ---------------------------------------------------------------
#
# a flavor of mapper based on projection to a single coordinate

#' Run 1D mapper
#'
#' Run mapper using a one-dimensional filter, a cover of intervals, and a clustering algorithm.
#'
#' @param data A data frame.
#' @param dists A distance matrix for the data frame.
#' @param filtered_data The result of a function applied to the data frame; there should be one filter value per observation in the original data frame.
#' @param cover A 2D array of interval left and right endpoints; rows should be intervals and columns left and right endpoints (in that order).
#' @param clustering_method A string to pass to [hclust] to determine clustering method.
#' @param global_clustering Whether you want clustering to happen in a global (all level visible) or local (only current level set visible) context.
#'
#' @return A list of two data frames, one with node data containing bin membership,
#'  data points per cluster, and cluster dispersion, and one with edge data
#'  containing sources, targets, and weights representing overlap strength.
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
#' create_1D_mapper_object(data, dist(data), projx, cover, "single")
create_1D_mapper_object <- function(data,
                                    dists,
                                    filtered_data,
                                    cover,
                                    clustering_method = "single",
                                    global_clustering = TRUE) {
  if (!all(cover[, 1] - cover[, 2] <= 0)) {
    stop("left endpoints must be less than or equal to right endpoints")
  }

  cover = apply(cover, 1, check_in_interval)

  return(create_mapper_object(data, dists, filtered_data, cover, clustering_method, global_clustering))
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
#' @param clustering_method A string to pass to [hclust] to determine clustering method.
#' @param global_clustering Whether you want clustering to happen in a global (all level visible) or local (only current level set visible) context.
#'
#' @return A list of two dataframes, one with node data containing bin membership,
#'  datapoints per cluster, and cluster dispersion, and one with edge data
#'  containing sources, targets, and weights representing overlap strength.
#' @export
#' @examples
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#' data.dists = dist(data)
#' eps = 1
#'
#' create_clusterball_mapper_object(data, data.dists, data.dists, eps, "single")
create_clusterball_mapper_object <- function(data, dist1, dist2, eps, clustering_method, global_clustering = TRUE) {
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

  return(create_mapper_object(
    data,
    dist2,
    rownames(data),
    lapply(balls, is_in_ball),
    clustering_method,
    global_clustering
  ))
}
