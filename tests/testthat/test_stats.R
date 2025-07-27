# # generate 1000 random data points in the unit cube
data = data.frame(x = runif(1000), y = runif(1000), z = runif(1000))
dists = as.matrix(dist(data))
min_sum = min(apply(dists, 1, sum))
cluster = convert_to_clusters(row.names(data))
medoid = get_cluster_medoid(dists, cluster)

test_that("medoid has minimum sum", {
  expect_true(sum(dists[medoid, ]) == min_sum)
})

test_that("maximum pairwise distance is maximum", {
  expect_true(max(dists) == get_width(dists, cluster))
})

test_that("max distance to medoid is maximum", {
  expect_true(max(dists[medoid, ]) == get_max_eccentricity(dists, cluster))
})
