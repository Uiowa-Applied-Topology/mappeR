#=========================================================
# Darkroom: for lenscrafting
#=========================================================


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
#' create_width_balanced_cover(0, 100, 10, 15)
#' create_width_balanced_cover(-11.5, 10.33, 100, 2)
create_width_balanced_cover <- function(min_val, max_val, num_bins, percent_overlap) {

  even_length = (max_val - min_val)/num_bins # widths with zero percent overlap
  nudge = (percent_overlap/100)*even_length

  left_ends = min_val + (even_length-nudge)*(0:(num_bins-1)) # construct correctly overlapping bins
  right_ends = left_ends + even_length # we will scale everything after

  scale_factor = (max_val - min_val)/(right_ends[num_bins] - min_val) # scale by pretending min_val = 0
  bin_ends = cbind(left_ends, right_ends) # make bins

  bin_ends = scale_factor*(bin_ends - min_val) + min_val # translate to zero, scale, then translate back

  return(bin_ends)
}
