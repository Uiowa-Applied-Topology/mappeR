
# 1D covering -------------------------------------------------------------

# get two random numbers between -1 and 1 for the left and right endpoints
vars = sort(runif(2, -1, 1))
left_end  = vars[1]
right_end = vars[2]

# ethan's breaking case
ethan_left = 0
ethan_right = 1

# get a random number of bins between 1 and 100
num_bins = 5
ethan_bins = 5

# get a random percent overlap between 0 and 100
percent_overlap = runif(1)*100
ethan_overlap = 40

# create cover with random parameters
no_overlap_cover = create_width_balanced_cover(left_end, right_end, num_bins, 0)
full_overlap_cover = create_width_balanced_cover(left_end, right_end, num_bins, 100)
cover = create_width_balanced_cover(left_end, right_end, num_bins, percent_overlap)
ethan_cover = create_width_balanced_cover(ethan_left, ethan_right, ethan_bins, ethan_overlap)

# length of the first interval
no_length = as.numeric(abs(no_overlap_cover[1,1] - no_overlap_cover[1,2]))
full_length = as.numeric(abs(full_overlap_cover[1,1] - full_overlap_cover[1,2]))
interval_length = as.numeric(abs(cover[1,1] - cover[1,2]))
ethan_length = as.numeric(abs(ethan_cover[1,1] - ethan_cover[1,2]))

test_that("1D width balanced cover endpoints are correct", {
  expect_true(abs(cover[1,1] - left_end) < 1e-15)
  expect_true(abs(right_end - cover[nrow(cover), 2]) < 1e-15)

  expect_true(abs(ethan_cover[1,1] - ethan_left) < 1e-15)
  expect_true(abs(ethan_right - ethan_cover[nrow(ethan_cover), 2]) < 1e-15)
})

test_that("1D width balanced cover intervals are equally sized", {
  # each interval in the cover should be the same length, so the same as the first interval
  if (num_bins > 1) {
    sapply(1:(num_bins - 1), function(i) expect_equal(as.numeric(abs(cover[i, 1] - cover[i, 2])), interval_length))
  }
  sapply(1:(ethan_bins - 1), function(i) expect_equal(as.numeric(abs(ethan_cover[i, 1] - ethan_cover[i, 2])), ethan_length))
})

test_that("1D width balanced cover overlaps are correct and consistent", {
  if (num_bins > 1) {
    # check that bins have correct overlap with the bin behind them
    sapply(2:num_bins, function(i) expect_equal(as.numeric(100*abs(cover[i, 1] - cover[i - 1, 2])/interval_length), percent_overlap))
    sapply(2:ethan_bins, function(i) expect_equal(as.numeric(100*abs(ethan_cover[i, 1] - ethan_cover[i - 1, 2])/ethan_length), ethan_overlap))

    # check that bins have correct overlap with the bin ahead of them
    sapply(1:(num_bins - 1), function(i) expect_equal(as.numeric(100*abs(cover[i, 2] - cover[i+1, 1])/interval_length), percent_overlap))
    sapply(1:(ethan_bins - 1), function(i) expect_equal(as.numeric(100*abs(ethan_cover[i, 2] - ethan_cover[i+1, 1])/ethan_length), ethan_overlap))

    # first and last bin overlaps should also be correct
    expect_equal(as.numeric(100*abs(cover[1, 2] - cover[2, 1])/interval_length), percent_overlap)
    expect_equal(as.numeric(100*abs(cover[nrow(cover), 1] - cover[nrow(cover)-1, 2])/interval_length), percent_overlap)

    expect_equal(as.numeric(100*abs(ethan_cover[1, 2] - ethan_cover[2, 1])/ethan_length), ethan_overlap)
    expect_equal(as.numeric(100*abs(ethan_cover[nrow(ethan_cover), 1] - ethan_cover[nrow(ethan_cover)-1, 2])/ethan_length), ethan_overlap)
  }
})

test_that("1D width balanced cover can handle no or full overlap", {
  expect_true(abs(no_overlap_cover[1,1] - left_end) < 1e-15)
  expect_true(abs(right_end - no_overlap_cover[nrow(no_overlap_cover), 2]) < 1e-15)

  expect_true(abs(full_overlap_cover[1,1] - left_end) < 1e-15)
  expect_true(abs(right_end - full_overlap_cover[nrow(full_overlap_cover), 2]) < 1e-15)

  if (num_bins > 1) {
    sapply(1:(num_bins - 1), function(i) expect_equal(as.numeric(abs(cover[i, 1] - cover[i, 2])), interval_length))
    sapply(1:(num_bins - 1), function(i) expect_equal(as.numeric(abs(no_overlap_cover[i, 1] - no_overlap_cover[i, 2])), no_length))
    sapply(1:(num_bins - 1), function(i) expect_equal(as.numeric(abs(full_overlap_cover[i, 1] - full_overlap_cover[i, 2])), full_length))
    sapply(2:num_bins, function(i) expect_equal(as.numeric(100*abs(no_overlap_cover[i, 1] - no_overlap_cover[i - 1, 2])/no_length), 0))
    sapply(2:num_bins, function(i) expect_equal(as.numeric(100*abs(full_overlap_cover[i, 1] - full_overlap_cover[i - 1, 2])/full_length), 100))
    sapply(1:(num_bins - 1), function(i) expect_equal(as.numeric(100*abs(no_overlap_cover[i, 2] - no_overlap_cover[i+1, 1])/no_length), 0))
    sapply(1:(num_bins - 1), function(i) expect_equal(as.numeric(100*abs(full_overlap_cover[i, 2] - full_overlap_cover[i+1, 1])/full_length), 100))
    expect_equal(as.numeric(100*abs(no_overlap_cover[1, 2] - no_overlap_cover[2, 1])/no_length), 0)
    expect_equal(as.numeric(100*abs(full_overlap_cover[1, 2] - full_overlap_cover[2, 1])/full_length), 100)
    expect_equal(as.numeric(100*abs(no_overlap_cover[nrow(no_overlap_cover), 1] - no_overlap_cover[nrow(no_overlap_cover)-1, 2])/no_length), 0)
    expect_equal(as.numeric(100*abs(full_overlap_cover[nrow(full_overlap_cover), 1] - full_overlap_cover[nrow(full_overlap_cover)-1, 2])/full_length), 100)
  }
})
