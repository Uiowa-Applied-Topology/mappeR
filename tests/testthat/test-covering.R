
# 1D covering -------------------------------------------------------------

# get two random numbers between -1 and 1 for the left and right endpoints
vars = sort(runif(2, -1, 1))
left_end  = vars[1]
right_end = vars[2]

# get a random number of bins between 1 and 100
num_bins = sample(1:100, 1)

# get a random percent overlap between 0 and 100
percent_overlap = runif(1)*100

# create cover with random parameters
cover = create_width_balanced_cover(left_end, right_end, num_bins, percent_overlap)

# length of the first interval
interval_length = as.numeric(abs(cover[1,1] - cover[1,2]))

test_that("1D width balanced cover endpoints are correct", {
  expect_true(left_end == cover[1,1])
  expect_true(right_end == cover[nrow(cover), 2])
})

test_that("1D width balanced cover intervals are equally sized", {
  # each interval in the cover should be the same length, so the same as the first interval
  sapply(1:(num_bins - 1), function(i) expect_equal(as.numeric(abs(cover[i, 1] - cover[i, 2])), interval_length))
})

test_that("1D width balanced cover overlaps are correct and consistent", {
  # check that bins have correct overlap with the bin behind them
  sapply(2:num_bins, function(i) expect_equal(as.numeric(100*abs(cover[i, 1] - cover[i - 1, 2])/interval_length), percent_overlap))

  # check that bins have correct overlap with the bin ahead of them
  sapply(1:(num_bins - 1), function(i) expect_equal(as.numeric(100*abs(cover[i, 2] - cover[i+1, 1])/interval_length), percent_overlap))

  # first and last bin overlaps should also be correct
  expect_equal(as.numeric(100*abs(cover[1, 2] - cover[2, 1])/interval_length), percent_overlap)
  expect_equal(as.numeric(100*abs(cover[nrow(cover), 1] - cover[nrow(cover)-1, 2])/interval_length), percent_overlap)
})
