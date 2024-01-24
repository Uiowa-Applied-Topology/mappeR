create_balls <- function(data, dists, eps) {
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

# takes the output of the previous function and makes it suitable for the 1D mapper function
convert_balls <- function(balled_data) {
  res = c()
  names = c()
  count = 0
  for (i in 1:length(balled_data)) {
    names = append(names, names(balled_data[[i]]))
    for (j in 1:length(balled_data[[i]])) {
      res = append(res, i)
    }
  }
  names(res) = names
  return(res)
}

construct_graph <- function(clustered_data) {
  num_vertices = max(clustered_data[[length(clustered_data)]]) # I don't know why this works

  flattened_data = unlist(clustered_data)

  amat = matrix(, nrow = num_vertices, ncol = num_vertices)

  overlap_vector = c()

  for (i in 1:(num_vertices-1)) {
    for (j in i:num_vertices) {
      if (i == j) {
        amat[i, j] = 0
      } else {
        my_cluster = flattened_data[flattened_data == i] # get the datapoints in the ith cluster
        # my_cluster.length = length(my_cluster)
        compare_cluster = flattened_data[flattened_data == j] # get the datapoints in the jth cluster
        # compare_cluster.length = length(compare_cluster)
        overlap = intersect(names(my_cluster), names(compare_cluster))
        overlap.length = length(overlap)
        # avg_overlap = .5*(overlap.length*(my_cluster.length + compare_cluster.length)/(my_cluster.length * compare_cluster.length))
        if (length(overlap) != 0) {
          amat[i, j] = 1
          overlap_vector = append(overlap_vector, overlap.length)
        } else {
          amat[i, j] = 0
        }
      }
    }
  }
  return(list(amat, overlap_vector))
}

# runner function for 1D mapper; outputs bins, clusters, and the mapper graph.
get_mapper_data <- function(data, dists, eps) {
  # bin data according to cool epsilon net thing
  print("binning...")
  binned_data = create_balls(data, dists, eps)

  # s a n i t i z e
  print("clustering...")
  clustered_data = convert_balls(binned_data)

  # construct mapper graph
  print("making mapper graph...")
  graph_data = construct_graph(clustered_data)
  amat = graph_data[[1]]
  edge_overlaps = graph_data[[2]]
  mapper_graph = graph_from_adjacency_matrix(amat, mode="max")

  return(list(clustered_data, mapper_graph, edge_overlaps))
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

visualize_mapper_data <- function(mapper_data) {
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
  visualize_mapper_data(get_mapper_data(data, dists, eps))
  return(invisible(NULL))
}
