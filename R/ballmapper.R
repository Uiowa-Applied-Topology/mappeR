source("R/graph_constructor.R")

# greedy epsilon net algorithm from Dłotko
# output is a list of bins, each containing names of datapoints.
create_balls <- function(data, dists, eps) {
  dists = as.matrix(dists) # because I am stupid and usedists isn't working we use a symmetric matrix
  balls = list()
  marked = rep(FALSE, nrow(data)) # keep track of which points we've covered
  datanames = rownames(data) # actually keep track of the data

  names(marked) = datanames

  while (FALSE %in% marked) { # keep going until we have covered all the data
    current_ball_center = NULL

    # find a ball center
    if (length(which(marked)) == 0) {
      current_ball_center = sample(datanames, 1) # pick a random point if no points are marked
    } else {
      unmarked_points = datanames[which(!marked)]
      current_ball_center = sample(unmarked_points, 1) # otherwise pick from the set of unmarked points
    }
    all_dists = dists[current_ball_center,] # get all distances away from ball center
    balled_data_names = datanames[which(all_dists < eps)] # restrict to within the (open???) ball
    marked[balled_data_names] = TRUE # mark points inside the ball as covered
    balls = append(balls, list(balled_data_names)) # add the ball to our big list of balls
  }
  return(balls)
}

# takes the output of the previous function and makes it suitable for the 1D mapper function
convert_balls <- function(balled_data) {
  # we will stitch these together to get a named vector of ball/cluster membership
  res = c()
  names = c()

  for (i in 1:length(balled_data)) {
    current_ball = balled_data[[i]]
    names = append(names, current_ball) # collect data ids from this ball
    for (j in 1:length(current_ball)) {
      res = append(res, i) # each datapoint will get an appropriate cluster membership tag
    }
  }
  names(res) = names # sewing time
  return(res)
}

# gets the ballmapper data: ball membership, graph structure, and cluster overlap information
get_ballmapper_data <- function(data, dists, eps) {
  print("making balls...")
  balled_data = create_balls(data, dists, eps)
  formatted_balled_data = convert_balls(balled_data)

  # construct ballmapper graph
  print("constructing ballmapper graph...")
  graph_data = construct_graph(formatted_balled_data) # gets adjacency matrix and edge overlaps
  amat = graph_data[[1]]
  edge_overlaps = graph_data[[2]]
  mapper_graph = graph_from_adjacency_matrix(amat, mode="max") # so we only have to record the upper half

  return(list(formatted_balled_data, mapper_graph, edge_overlaps))
}
