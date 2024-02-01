# greedy epsilon net algorithm from Dłotko
create_balls <- function(data, dists, eps) {
  dists = as.matrix(dists) # because I am stupid and usedists isn't working we use a symmetric matrix
  balls = list()
  marked = rep(FALSE, nrow(data)) # keep track of which points we've covered
  names(marked) = rownames(dists) # actually keep track of the data
  datanames = rownames(data) # for easy grabbing of unmarked points

  while (FALSE %in% marked) { # keep going until we have covered all the data
    current_ball_center = NULL

    # find a ball center
    if (length(which(marked)) == 0) {
      current_ball_center = sample(rownames(data), 1) # pick a random point if no points are marked
    } else {
      unmarked_points = datanames[which(!marked)]
      current_ball_center = sample(unmarked_points, 1) # otherwise pick from the set of unmarked points
    }

    all_dists = dists[current_ball_center,] # get all distances away from ball center
    balled_dists = all_dists[all_dists < eps] # restrict to within the (open???) ball
    marked[names(balled_dists)] = TRUE # mark points inside the ball as covered
    balls = append(balls, list(balled_dists)) # add the ball to our big list of balls
  }
  return(balls)
}

# takes the output of the previous function and makes it suitable for the 1D mapper function
convert_balls <- function(balled_data) {
  # we will stitch these together to get a named vector of ball/cluster membership
  res = c()
  names = c()

  print(length(balled_data))

  for (i in 1:length(balled_data)) {
    print(i)
    names = append(names, names(balled_data[[i]])) # collect data ids from this ball
    for (j in 1:length(balled_data[[i]])) {
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

get_size_vector <- function(clustered_data, num_vertices) {
  size_vector = c()

  flattened_data = unlist(clustered_data)

  for (i in 1:num_vertices) {
    my_cluster = flattened_data[flattened_data == i]
    size_vector = append(size_vector, length(my_cluster))
  }

  size_vector = (size_vector/sqrt(sum(size_vector^2)))*200

  return(size_vector)
}

visualize_ballmapper_data <- function(mapper_data) {
  clustered_data = mapper_data[[1]]
  mapper_graph = mapper_data[[2]]
  edge_weights = mapper_data[[3]]

  num_vertices = gorder(mapper_graph)
  num_edges = gsize(mapper_graph)
  num_bins = length(clustered_data)

  cygraph = set_vertex_attr(mapper_graph, "cluster", value = 1:num_vertices)
  cygraph = set_edge_attr(cygraph, "overlap", value = (edge_weights/sqrt(sum(edge_weights^2)))*25)
  cygraph = set_vertex_attr(cygraph, "cluster_size", value = get_size_vector(clustered_data, num_vertices))

  createNetworkFromIgraph(cygraph)

  style.name = "mapperstyle"
  defaults <- list(NODE_SHAPE = "ellipse",
                   EDGE_TRANSPARENCY = 120)

  nodeLabels <- mapVisualProperty('node label', 'cluster', 'p')
  nodeSizes <- mapVisualProperty('node size', 'cluster_size', 'p')
  edgeWidth <- mapVisualProperty('edge width', 'overlap', 'p')

  createVisualStyle(style.name, defaults, list(nodeLabels, nodeSizes, edgeWidth))

  setVisualStyle(style.name)

}

cyballmapper <- function(data, dists, eps) {
  visualize_ballmapper_data(get_ballmapper_data(data, dists, eps))
  return(invisible(NULL))
}
