
# 1D binning --------------------------------------------------------------

generate_random_intervals <- function(num_bins, min_val, max_val) {
  return(t(replicate(num_bins, sort(runif(2, min_val, max_val)))))
}

num_intervals = round(100*runif(1)) + 1 # to avoid zero case

# generate 1000 random data points between -2 and 2
data = data.frame(x = runif(1000, -2, 2), drop = FALSE)

test_that("data is binned properly into 1D cover elements", {
  intervals = generate_random_intervals(num_intervals, -3, 3)
  bins = create_bins(data, data$x, apply(intervals, 1, check_in_interval))

  sapply(1:length(bins), function(i) expect_true(all((intervals[i,2] - data[bins[[i]],1] >= 0) & (data[bins[[i]],1] - intervals[i,1] >= 0))))
})

test_that("data with endpoints outside the cover binned properly", {
  intervals = generate_random_intervals(num_intervals, -1, 1)
  bins = create_bins(data, data$x, apply(intervals, 1, check_in_interval))

  sapply(1:length(bins), function(i) expect_true(all((intervals[i,2] - data[bins[[i]],1] >= 0) & (data[bins[[i]],1] - intervals[i,1] >= 0))))
})
