###########################################################################
# VIEWMASTER LIBRARY
# lens/filter functions
###########################################################################

#' Eccentric Lens
#'
#' Compute eccentricity of data points.
#'
#' @param dists A distance matrix associated to a data frame. Can be a `dist` object or `matrix`.
#'
#' @return A vector of centrality measures, calculated per data point as the sum of its distances to every other data point, divided by the number of points.
#' @export
#'
#' @examples
#' # Generate some noisy data along a 2D curve
#' num_points = 100
#' P.data = data.frame(
#'   x = sapply(1:num_points, function(x)
#'     sin(x) * 10) + rnorm(num_points, 0, 0.1),
#'   y = sapply(1:num_points, function(x)
#'     cos(x) ^ 2 * sin(x) * 10) + rnorm(num_points, 0, 0.1)
#' )
#' P.dist = dist(P.data)
#'
#' # Apply eccentricity filter
#' eccentricity = eccentricity_filter(P.dist)
eccentricity_filter <- function(dists) {
  dists = as.matrix(dists)
  apply(dists, 1, sum) / nrow(dists)
}
