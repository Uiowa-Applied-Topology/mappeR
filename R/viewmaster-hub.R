#' Compute eccentricity of data points

#' @param dists A distance matrix associated to a data frame.
#'
#' @return A vector of centrality measures, calcuated per data point as the sum of its distances to every other data point, divided by the number of points.
#' @export
#'
#' @examples
#' num_points = 5000
#'
#' P.data = data.frame(
#'   x = sapply(1:num_points, function(x)
#'     sin(x) * 10) + rnorm(num_points, 0, 0.1),
#'   y = sapply(1:num_points, function(x)
#'     cos(x) ^ 2 * sin(x) * 10) + rnorm(num_points, 0, 0.1),
#'   z = sapply(1:num_points, function(x)
#'     10 * sin(x) ^ 2 * cos(x)) + rnorm(num_points, 0, 0.1)
#' )
#'
#' P.dist = dist(P.data)
#' eccentricity = eccentricity_filter(P.dist)
eccentricity_filter <- function(dists) {
  dists = as.matrix(dists)
  apply(dists, 1, sum, simplify = FALSE) / nrow(dists)
}
