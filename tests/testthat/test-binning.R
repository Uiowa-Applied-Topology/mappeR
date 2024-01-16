test_that("interval math works", {
  vars = c(0, 1, 5, 40)
  left_end  = vars[1]
  right_end = vars[2]
  num_bins = floor(vars[3]*100)
  percent_overlap = vars[4]*100

  bins = get_width_balanced_endpoints(left_end, right_end, num_bins, percent_overlap)

  expect_equal(nrow(bins), num_bins) # output should contain correct number of bins
  expect_equal(as.numeric(bins[1,1]), left_end) # leftmost endpoint should be correct
  expect_equal(as.numeric(bins[num_bins, 2]), right_end) # rightmost endpoint should be correct

  bin_length = as.numeric(abs(bins[1,1] - bins[1,2]))

  for (i in 2:(num_bins - 1)) {
    expect_equal(as.numeric(abs(bins[i, 1] - bins[i, 2])), bin_length) # bin lengths should be consistent

    left_overlap_percent = as.numeric(100*abs(bins[i, 1] - bins[i-1, 2])/bin_length)
    right_overlap_percent = as.numeric(100*abs(bins[i, 2] - bins[i+1, 1])/bin_length)

    expect_equal(left_overlap_percent, percent_overlap) # overlap on the left is correct
    expect_equal(right_overlap_percent, percent_overlap) # overlap on the right is correct
  }

  expect_equal(as.numeric(100*abs(bins[1, 2] - bins[2, 1])/bin_length), percent_overlap) # first bin overlap is correct
  expect_equal(as.numeric(100*abs(bins[nrow(bins), 1] - bins[nrow(bins)-1, 2])/bin_length), percent_overlap) # last bin overlap is correct
})

# idea from https://stackoverflow.com/questions/32345484/does-vector-exist-in-matrix
test_that("data binning works", {
  random_left_ends = runif(20)*100
  bin_ends = cbind(random_left_ends, random_left_ends + runif(20)*100)
  data = data.frame(x = runif(50)*100, y= runif(50)*100)
  bins = make_bins(data, data$x, bin_ends)
  for (data_idx in 1:50) {
    for (bin_idx in 1:20) {
      if ((data$x[data_idx] >= bin_ends[bin_idx, 1]) & (data$x[data_idx] <= bin_ends[bin_idx, 2])) {
        expect_true(any(bins[[bin_idx]] == data[data_idx,][col(bins[[bin_idx]])])) # datapoint should correctly be in appropriate bin
      } else {
        if(length(bins[[bin_idx]]) != 0) {
          expect_false(any(bins[[bin_idx]] == data[data_idx,][col(bins[[bin_idx]])])) # should also correctly not be in inappropriate bin
        }
      }
    }
  }
})
