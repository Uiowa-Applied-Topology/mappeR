test_that("igraph works with no overlap", {
  data = data.frame(x=1:100, drop=FALSE)
  cover = create_width_balanced_cover(1, 100, 10, 0)
  mapperobj = create_1D_mapper_object(data, dist(data), data$x, cover, "single")
  expect_no_warning(mapper_object_to_igraph(mapperobj))
})
