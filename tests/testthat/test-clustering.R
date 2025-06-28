data = data.frame(x = runif(1000, -2, 2), drop= FALSE)
lens = data$x
names(lens) = row.names(data)
cover = create_width_balanced_cover(min(lens), max(lens), 10, 25)
dists = dist(data)


test_that("global clustering", {
  expect_no_warning(create_1D_mapper_object(data, dists, lens, cover, global_hierarchical_clusterer("single", dists)))
})

test_that("local clustering", {
  expect_no_warning(create_1D_mapper_object(data, dists, lens, cover, local_hierarchical_clusterer("complete")))
})
