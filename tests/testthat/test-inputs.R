# # generate 1000 random data points between -2 and 2
data = data.frame(x = runif(1000, -2, 2), drop = FALSE)

names = runif(length(data$x), 0, 1)
row.names(data) = names # for testing
one_dim_projection = data$x
names(one_dim_projection) = names

cover = create_width_balanced_cover(min(one_dim_projection), max(one_dim_projection), sample(1:50, 1), sample(0:100, 1))

test_that("we can clusterball with clusterball", {
  expect_no_warning(create_clusterball_mapper_object(data, dist(data), dist(data), .3))
})

test_that("mapper happens ok with distance matrix as a matrix", {
  expect_no_warning(create_1D_mapper_object(data, as.matrix(dist(data)), one_dim_projection, cover))
})

test_that("mapper happens ok with distance matrix as a dist", {
  expect_no_warning(create_1D_mapper_object(data, dist(data), one_dim_projection, cover))
})

test_that("mapper is ok with no clustering method", {
  expect_no_warning(create_1D_mapper_object(data, dist(data), one_dim_projection, cover))
})

test_that("mapper works with differently formatted filtered data", {
  expect_no_warning(create_1D_mapper_object(data, dist(data), as.data.frame(one_dim_projection), cover))
  expect_no_warning(create_1D_mapper_object(data, dist(data), as.vector(one_dim_projection), cover)) # removes names
  expect_no_warning(create_1D_mapper_object(data, dist(data), as.list(one_dim_projection), cover))
  expect_no_warning(create_1D_mapper_object(data, dist(data), as.matrix(one_dim_projection), cover))
})

test_that("shuffling data does not affect output", {
  unshuffled = create_1D_mapper_object(data, dist(data), one_dim_projection, cover)
  shuffled = create_1D_mapper_object(sample(data), dist(sample(data)), sample(one_dim_projection), cover[sample(1:nrow(cover)), ])

  # recall the symmetric difference of two sets is empty if and only the sets are the same
  symdiff = function(a, b) setdiff(union(a, b), intersect(a, b))

  # get the names of the data points per vertex
  unshuffled_data = lapply(unshuffled[[1]]$data, function(x) unlist(strsplit(x, ", ")))
  shuffled_data = lapply(shuffled[[1]]$data, function(x) unlist(strsplit(x, ", ")))

  # make sure the data in each unshuffled vertex matches the data is some shuffled vertex
  lapply(unshuffled_data, function(vertex) expect_true(any(unlist(lapply(shuffled_data, function(shuffled_vertex) length(symdiff(vertex, shuffled_vertex)) == 0)))))
})




