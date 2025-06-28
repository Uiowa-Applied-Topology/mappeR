
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

# I dunno how necessary this test is
test_that("data with endpoints outside the cover (not) binned properly", {
  intervals = generate_random_intervals(num_intervals, -1, 1)
  bins = create_bins(data, data$x, apply(intervals, 1, check_in_interval))
  sapply(1:length(bins), function(i) expect_true(all((intervals[i,2] - data[bins[[i]],1] >= 0) & (data[bins[[i]],1] - intervals[i,1] >= 0))))
})

data = data.frame(x = runif(1000, -2, 2), drop = FALSE)

test_that("two bins does not cause error", {
  intervals = generate_random_intervals(2, -1, 1)
  expect_no_warning(create_bins(data, data$x, apply(intervals, 1, check_in_interval)))
})

test_that("one bin does not cause error", {
  interval = generate_random_intervals(1, -1, 1)
  expect_no_warning(create_bins(data, data$x, list(check_in_interval(interval))))
})

test_that("more bins than data points does not cause error", {
  intervals = generate_random_intervals(nrow(data) + 1, -1, 1)
  expect_no_warning(create_bins(data, data$x, apply(intervals, 1, check_in_interval)))
})

test_that("every bin empty quits", {
  expect_error(create_1D_mapper_object(data, dist(data), one_dim_projection, create_width_balanced_cover(max(one_dim_projection) + 1, max(one_dim_projection) + 2, 10, runif(1)*100)))
})

test_that("we can ball with equal sized bins", {
  data = data.frame(x = 1:100, drop = FALSE)
  cover = create_width_balanced_cover(1, 100, 10, 0)
  projection = data$x
  names(projection) = row.names(data)
  bins = create_bins(data, projection, apply(cover, 1, check_in_interval))
  expect_no_warning(assemble_mapper_object(convert_to_clusters(bins), dist(data), binning = FALSE))
})

test_that("we can clusterball with equal sized bins", {
  data = data.frame(x = 1:100, drop = FALSE)
  cover = create_width_balanced_cover(1, 100, 10, 0)
  projection = data$x
  names(projection) = row.names(data)
  bins = create_bins(data, data$x, apply(cover, 1, check_in_interval))
  expect_no_warning(create_mapper_object(data, dist(data), projection, lapply(bins, is_in_ball)))
})

test_that("we can ball with ballmapper", {
  expect_no_warning(create_ball_mapper_object(data, dist(data), .3))
})
