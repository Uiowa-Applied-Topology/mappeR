create_bins <- function(data, dists, eps) {
  dists = as.matrix(dists) # because I am stupid and usedists isn't working
  balls = list()
  marked = rep(FALSE, nrow(data))
  names(marked) = rownames(dists)
  while (FALSE %in% marked) {
    current_ball_center = NULL
    if (length(which(marked)) == 0) {
      current_ball_center = sample(rownames(data), 1) # pick a random point if no points are marked

    } else {
      datanames = rownames(data)
      unmarked_points = datanames[which(!marked)]
      current_ball_center = sample(unmarked_points, 1)
    }
    all_dists = dists[current_ball_center,] # all distances away from ball center
    balled_dists = all_dists[all_dists < eps] # now restrict to in the ball
    marked[names(balled_dists)] = TRUE # mark points inside the ball as covered
    balls = append(balls, list(balled_dists)) # add to our bin list
  }
  return(balls)
}
