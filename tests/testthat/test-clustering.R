data = data.frame(x = runif(1000, -2, 2), drop= FALSE)
cover = create_width_balanced_cover(min(data$x), max(data$x), 10, 25)
dists = dist(data)


test_that("global clustering", {
  expect_no_warning(create_1D_mapper_object(data, dists, data$x, cover, global_hierarchical_clusterer("single", dists)))
})

test_that("local clustering", {
  expect_no_warning(create_1D_mapper_object(data, dists, data$x, cover, local_hierarchical_clusterer("complete")))
})
