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
  ball_sizes = lapply(balled_data, length)
  ballball_data = unlist(mapply(function(x, y) rep(x, y), 1:length(ball_sizes), ball_sizes))
  names(ballball_data) = unlist(balled_data)

  return(ballball_data)
}

construct_ballmappergraph <- function(binclust_data, dists) {
  num_vertices = max(binclust_data[[length(binclust_data)]])

  node_ids = as.character(1:num_vertices)

  overlaps = get_overlaps(binclust_data)
  edges = get_edgelist_from_overlaps(overlaps, num_vertices)
  sources = as.character(edges[,1])
  targets = as.character(edges[,2])

  cluster_tightness = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_size = get_cluster_sizes(binclust_data, num_vertices)
  data_in_cluster = unlist(get_clustered_data(binclust_data, num_vertices))
  edge_weights = get_edge_weights(sapply(overlaps, length), cluster_size, edges)

  nodes = data.frame(id=node_ids,
                     size=cluster_size,
                     tightness=cluster_tightness,
                     data=data_in_cluster)

  edges = data.frame(source=sources,
                     target=targets,
                     weight=edge_weights)


  return(list(nodes, edges))
}

#' Runs ballmapper and returns two dataframes containing node and edge data.
#'
#' @param data Your input data. Ideally a dataframe.
#' @param dists A distance matrix for your data. Can be a `dist` object or 2D matrix.
#' @param eps A positive real number for your desired ball radius.
#' @returns A list of two dataframes, one with node data containing ball size,
#'  datapoints per ball, ball tightness, and one with edge data
#'  containing sources, targets, and weights representing overlap strength.
#' @examples
#' circle.data = data.frame( x= sapply(1:1000, function(x) cos(x)) + rnorm(100, 500, .03),
#'   y = sapply(1:1000, function(x) sin(x)) + rnorm(100, 0, 0.03))
#' circle.dist = dist(circle.data)
#'
#' show(get_ballmapper_data(circle.data, circle.dist, .3))
#' @export
get_ballmapper_data <- function(data, dists, eps) {
  balled_data = create_balls(data, dists, eps)
  formatted_balled_data = convert_balls(balled_data)

  # construct ballmapper graph
  ballmappergraph = construct_ballmappergraph(formatted_balled_data, dists)

  return(ballmappergraph)
}

#' Runs ballmapper and passes data to Cytoscape for visualization.
#'
#' @param data Your input data. Ideally a dataframe.
#' @param dists A distance matrix for your data. Can be a `dist` object or 2D matrix.
#' @param eps A positive real number for your desired ball radius.
#' @returns NULL
#' @examples
#' circle.data = data.frame( x= sapply(1:1000, function(x) cos(x)) + rnorm(100, 500, .03),
#'   y = sapply(1:1000, function(x) sin(x)) + rnorm(100, 0, 0.03))
#' circle.dist = dist(circle.data)
#'
#' # make sure Cytoscape is running in the background or this will not work
#' # cyballmapper(circle.data, circle.dist, .3)
#' @export
cyballmapper <- function(data, dists, eps) {
  visualize_mapper_data(get_ballmapper_data(data, dists, eps))
  return(invisible(NULL))
}
