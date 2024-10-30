###########################################################################
# QUILTING HUT
# methods to make covers
###########################################################################

# 1D filtering --------------------------------------------------------------

#' Generate an overlapping cover of an interval
#'
#' This is a function that generates a cover of an interval \eqn{[a,b]} with
#' some number of (possibly) overlapping, evenly spaced, identical width subintervals.
#'
#' @param min_val The left endpoint \eqn{a}. A real number.
#' @param max_val The right endpoint \eqn{b}. A real number.
#' @param num_bins The number of cover intervals with which to cover the interval. A positive integer.
#' @param percent_overlap How much overlap desired between the cover intervals
#'  (the percent of the intersection of each interval with its immediate
#'   neighbor relative to its length, e.g., \eqn{[0,2]} and \eqn{[1,3]} would have \eqn{50\%} overlap).
#'   A real number between 0 and 100, inclusive.
#'
#' @return A 2D numeric array.
#' \itemize{
#'    \item left_ends - The left endpoints of the cover intervals.
#'    \item right_ends - The right endpoints of the cover intervals.
#' }
#' @export
#'
#' @examples
#' create_width_balanced_cover(min_val=0, max_val=100, num_bins=10, percent_overlap=15)
#' create_width_balanced_cover(-11.5, 10.33, 100, 2)
create_width_balanced_cover <- function(min_val,
                                        max_val,
                                        num_bins,
                                        percent_overlap) {

  if (any(is.na(c(min_val, max_val, num_bins, percent_overlap)))) {
    stop("cannot create a cover with an NA input!")
  }
  if ((num_bins <= 0) | (!(num_bins %% 1 == 0))) {
    stop("number of bins must be a positive integer")
  } else if ((percent_overlap < 0) | (percent_overlap > 100)) {
    stop("percent overlap must be between 0 and 100")
  }

  # calculate widths as if percent overlap is zero
  even_length = (max_val - min_val) / num_bins

  # calculate appropriate extension length to ends of zero overlap intervals
  nudge = (percent_overlap / 100) * even_length

  # construct intervals with the correct overlap but wrong width
  left_ends = min_val + (even_length - nudge) * (0:(num_bins - 1))
  right_ends = left_ends + even_length

  # scale intervals to correct width by translating to zero, scaling, then translating back
  scale_factor = (max_val - min_val) / (right_ends[num_bins] - min_val)
  bin_ends = cbind(left_ends, right_ends) # make bins
  bin_ends = scale_factor * (bin_ends - min_val) + min_val

  return(bin_ends)
}

#' Get a tester function for an interval.
#'
#' @param endpoints A vector of interval endpoints, namely a left and a right. Must be in order.
#'
#' @return A function that eats a data point and outputs TRUE if the datapoint is in the interval and FALSE if not.
check_in_interval <- function(endpoints) {
  return(function(x) (endpoints[1] - x <= 0) & (endpoints[2] - x >= 0))
}

# balls -------------------------------------------------------------------

# greedy epsilon net algorithm from Dłotko
# output is a list of bins, each containing names of datapoints.
#' Make a cover of balls
#'
#' @param data A data frame.
#' @param dists A distance matrix for the data frame.
#' @param eps A positive real number.
#'
#' @return A list of vectors of data point names, one list element per ball. The output is such that every data point is contained in a ball of radius \eqn{\varepsilon}, and no ball center is contained in more than one ball. The centers are datapoints themselves.
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
#' balls = create_balls(data = P.data, dists = P.dist, eps = .25)
create_balls <- function(data, dists, eps) {

  # safety first
  if (eps <= 0) {
    stop("epsilon needs to be positive")
  }

  dists = as.matrix(dists) # because I am stupid and usedists isn't working we use a symmetric matrix
  balls = list()
  marked = rep(FALSE, nrow(data)) # keep track of which points we've covered
  datanames = rownames(data) # actually keep track of the data

  names(marked) = datanames

  while (FALSE %in% marked) {
    # keep going until we have covered all the data
    current_ball_center = NULL

    # find a ball center
    if (length(which(marked)) == 0) {
      current_ball_center = sample(datanames, 1) # pick a random point if no points are marked
    } else {
      unmarked_points = datanames[which(!marked)]
      current_ball_center = sample(unmarked_points, 1) # otherwise pick from the set of unmarked points
    }
    all_dists = dists[current_ball_center, ] # get all distances away from ball center
    balled_data_names = datanames[which(all_dists < eps)] # restrict to within the (open???) ball
    marked[balled_data_names] = TRUE # mark points inside the ball as covered
    balls = append(balls, list(balled_data_names)) # add the ball to our big list of balls
  }
  return(balls)
}

#' Get a tester function for a ball.
#'
#' @param ball A list of data points.
#'
#' @return A function that eats a data point and returns TRUE or FALSE depending if the point is in the ball or not.
is_in_ball <- function(ball) {
  return(function(x) x %in% ball)
}
