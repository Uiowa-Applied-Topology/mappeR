# TODO: add n-D filtering
# TODO: add data-balanced covers
# TODO: add other types of filter functions

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
#'   An integer between 0 and 100, inclusive.
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
  even_length = (max_val - min_val) / num_bins # widths with zero percent overlap
  nudge = (percent_overlap / 100) * even_length

  left_ends = min_val + (even_length - nudge) * (0:(num_bins - 1)) # construct correctly overlapping bins
  right_ends = left_ends + even_length # we will scale everything after

  scale_factor = (max_val - min_val) / (right_ends[num_bins] - min_val) # scale by pretending min_val = 0
  bin_ends = cbind(left_ends, right_ends) # make bins

  bin_ends = scale_factor * (bin_ends - min_val) + min_val # translate to zero, scale, then translate back

  return(bin_ends)
}

#' Create a single bin of data
#'
#' This function identifies datapoints whose "filter" value falls within a numeric range.
#'
#' @param bin_left The minimum value to be filtered for.
#' @param bin_right The maximum value to be filtered for.
#' @param data A dataframe.
#' @param filtered_data A single column of the input dataframe.
#'
#' @return A vector of datapoint names which fall within the filter range.
make_one_bin <- function(bin_left, bin_right, data, filtered_data) {
  in_bin = sapply(filtered_data, function(x)
    (bin_left - x <= 0) & (bin_right - x >= 0))
  bin_assignments = which(in_bin) # tells us indices of binned datapoints
  if (length(bin_assignments) != 0) {
    return(rownames(data[bin_assignments, ])) # TODO: bother me about why I need the original dataset here, I think it's more safe but who knows!
  } else {
    return(vector()) # bin still exists, it's just empty
  }
}


#' Create bins of data from a 1D filter function
#'
#' This function identifies datapoints whose "filter" value falls within a specified set of numeric ranges.
#'
#' @param data A dataframe.
#' @param filtered_data A single column of the input dataframe.
#' @param bin_ends A 2D numeric array of input ranges, whose first column contains left endpoints and whose second column contains right endpoints.
#'
#' @return A list of "bins" containing vectors of names of datapoints.
make_bins <- function(data, filtered_data, bin_ends) {
  left_ends = bin_ends[, 1]
  right_ends = bin_ends[, 2]
  bins = mapply(
    make_one_bin,
    left_ends,
    right_ends,
    MoreArgs = list(data = data, filtered_data = filtered_data)
  )
  return(bins)
}

# balls -------------------------------------------------------------------

# greedy epsilon net algorithm from Dłotko
# output is a list of bins, each containing names of datapoints.
create_balls <- function(data, dists, eps) {
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

# takes the output of the previous function and makes it suitable for the 1D mapper function
convert_balls <- function(balled_data) {
  ball_sizes = lapply(balled_data, length)
  ballball_data = unlist(mapply(function(x, y)
    rep(x, y), 1:length(ball_sizes), ball_sizes))
  names(ballball_data) = unlist(balled_data)

  return(ballball_data)
}
