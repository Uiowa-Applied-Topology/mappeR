test_interval_math <- function(left_end, right_end, num_bins, percent_overlap) {
  bins = create_width_balanced_cover(left_end, right_end, num_bins, percent_overlap)

  expect_equal(nrow(bins), num_bins) # output should contain correct number of bins
  expect_equal(as.numeric(bins[1,1]), left_end) # leftmost endpoint should be correct
  expect_equal(as.numeric(bins[num_bins, 2]), right_end) # rightmost endpoint should be correct

  bin_length = as.numeric(abs(bins[1,1] - bins[1,2]))

  # TODO: probably vectorize this idk
  for (i in 2:(num_bins - 1)) {
    expect_equal(as.numeric(abs(bins[i, 1] - bins[i, 2])), bin_length) # bin lengths should be consistent

    left_overlap_percent = as.numeric(100*abs(bins[i, 1] - bins[i-1, 2])/bin_length)
    right_overlap_percent = as.numeric(100*abs(bins[i, 2] - bins[i+1, 1])/bin_length)

    expect_equal(left_overlap_percent, percent_overlap) # overlap on the left is correct
    expect_equal(right_overlap_percent, percent_overlap) # overlap on the right is correct
  }

  expect_equal(as.numeric(100*abs(bins[1, 2] - bins[2, 1])/bin_length), percent_overlap) # first bin overlap is correct
  expect_equal(as.numeric(100*abs(bins[nrow(bins), 1] - bins[nrow(bins)-1, 2])/bin_length), percent_overlap) # last bin overlap is correct
}

# TODO: more intervals?
# TODO: degenerate intervals, zero bins, zero overlap, forbidden inputs
test_that("interval math works", {
  # random case
  vars = c(sort(runif(2)))
  left_end  = vars[1]*100
  right_end = vars[2]*100
  num_bins = sample(1:1000, 1)
  percent_overlap = runif(1)*100
  test_interval_math(left_end, right_end, num_bins, percent_overlap)

  # case that broke ethan's code
  test_interval_math(0, 1, 5, 40)
})

# TODO: test empty bins, singletons, doubletons
# TODO: clean this up? why a 2D dataframe?
# idea from https://stackoverflow.com/questions/32345484/does-vector-exist-in-matrix
test_that("data binning works", {
  # 100 random 2D points, each coordinate with value from 0-100
  random_data = data.frame(x = runif(100)*100, y = runif(100)*100)

  min_x = min(random_data$x)
  max_x = max(random_data$x)

  # 20 bins (not evenly spaced), endpoints include data extremes so as to cover
  random_left_ends = runif(20)*100
  random_right_ends = c(random_left_ends + runif(20)*100, max_x)
  bin_ends = cbind(c(min_x, random_left_ends), random_right_ends)

  bins = make_bins(random_data, random_data$x, bin_ends)

  # TODO: vectorize this
  for (i in 1:20) {
    for (datapoint in bins[[i]]) {
      if ((random_data[datapoint,1] >= bin_ends[i, 1]) & (random_data[datapoint,1] <= bin_ends[i, 2])) {
        expect_true(any(bins[[i]] == datapoint)) # datapoint should correctly be in appropriate bin
      } else {
        if(length(bins[[i]]) != 0) {
          expect_false(any(bins[[i]] == datapoint)) # should also correctly not be in inappropriate bin
        }
      }
    }
  }
})
