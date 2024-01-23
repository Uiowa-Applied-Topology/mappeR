create_cover_vector <- function(data, dists, eps) {
  cover_vector = list()
  marked = rep(FALSE, nrow(data))
  while(FALSE %in% marked) { # while points are still uncovered
    ballbin = list()
    uncovered_idxs = which(marked)
    random_point = sample(data[uncovered_idxs], 1)$name


  }
}
