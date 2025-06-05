# # generate 1000 random data points between -2 and 2
data = data.frame(x = runif(1000, -2, 2), drop = FALSE)
cover = create_width_balanced_cover(min(data$x), max(data$x), 10, 25)

test_that("we can clusterball with equal sized bins", {
  data = data.frame(x = 1:100, drop = FALSE)
  cover = create_width_balanced_cover(1, 100, 10, 0)
  bins = create_bins(data, data$x, apply(cover, 1, check_in_interval))
  expect_no_warning(create_mapper_object(data, dist(data), data$x, lapply(bins, is_in_ball)))
})

test_that("we can clusterball with clusterball", {
  expect_no_warning(create_clusterball_mapper_object(data, dist(data), dist(data), .3))
})

test_that("mapper happens ok with distance matrix as a matrix", {
  expect_no_warning(create_1D_mapper_object(data, as.matrix(dist(data)), data$x, cover))
})

test_that("mapper happens ok with distance matrix as a dist", {
  expect_no_warning(create_1D_mapper_object(data, dist(data), data$x, cover))
})

test_that("mapper is ok with no clustering method", {
  expect_no_warning(create_1D_mapper_object(data, dist(data), data$x, cover))
})

test_that("mapper works with differently formatted filtered data", {
  expect_no_warning(create_1D_mapper_object(data, dist(data), as.data.frame(data$x), cover))
  expect_no_warning(create_1D_mapper_object(data, dist(data), as.list(data$x), cover))
  expect_no_warning(create_1D_mapper_object(data, dist(data), as.vector(data$x), cover))
  expect_no_warning(create_1D_mapper_object(data, dist(data), as.matrix(data$x), cover))
})

test_that("we can ball with equal sized bins", {
  data = data.frame(x = 1:100, drop = FALSE)
  cover = create_width_balanced_cover(1, 100, 10, 0)
  bins = create_bins(data, data$x, apply(cover, 1, check_in_interval))
  expect_no_warning(assemble_mapper_object(convert_to_clusters(bins), dist(data), binning = FALSE))
})

test_that("we can ball with ballmapper", {
  expect_no_warning(create_ball_mapper_object(data, dist(data), .3))
})

test_that("two bins does not cause error", {
  expect_no_warning(create_1D_mapper_object(data, dist(data), data$x, create_width_balanced_cover(min(data$x), max(data$x), 2, 0)))
})

test_that("one bin does not cause error", {
  expect_no_warning(create_1D_mapper_object(data, dist(data), data$x, create_width_balanced_cover(min(data$x), max(data$x), 1, runif(1)*100)))
})

test_that("bad bins quits", {
  expect_error(create_1D_mapper_object(data, dist(data), data$x, create_width_balanced_cover(max(data$x) + 1, max(data$x) + 2, 10, runif(1)*100)))
})

test_that("we can hierarchically cluster differently", {
  expect_no_warning(create_1D_mapper_object(data, dist(data), data$x, cover, clusterer = local_hierarchical_clusterer("mcquitty")))
})


