
# mapper mapper -----------------------------------------------------------

create_mapper_object <- function(data, dists, filtered_data, cover_element_tests, method="none") {
  bins = create_bins(data, filtered_data, cover_element_tests)

  if (method == "none") {
    print("WE'RE NOT READY FOR THAT YET AAAAAAAA")
    quit()
  } else {
    clusters = get_clusters(bins, dists, method)
    return(run_mapper(clusters, dists, binning = TRUE))
  }
}

create_single_bin <- function(data, filtered_data, cover_element_test) {
  in_bin = sapply(filtered_data, cover_element_test)
  bin_assignments = which(in_bin)
  if (length(bin_assignments) != 0) {
    return(rownames(data[bin_assignments, ])) # TODO: bother me about why I need the original dataset here, I think it's more safe but who knows!
  } else {
    return(vector()) # bin still exists, it's just empty
  }
}

create_bins <- function(data, filtered_data, cover_element_tests) {
  return(mapply(create_single_bin, cover_element_test = cover_element_tests, MoreArgs = list(data = data, filtered_data = filtered_data)))
}

run_mapper <- function(binclust_data, dists, binning=TRUE) {
  num_vertices = max(binclust_data[[length(binclust_data)]])
  node_ids = as.character(1:num_vertices)
  overlaps = get_overlaps(binclust_data)
  edges = get_edgelist_from_overlaps(overlaps, num_vertices)
  sources = as.character(edges[, 1])
  targets = as.character(edges[, 2])

  # calculate some cluster stats
  cluster_tightness = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_size = get_cluster_sizes(binclust_data, num_vertices)
  data_in_cluster = unlist(get_clustered_data(binclust_data, num_vertices))
  edge_weights = get_edge_weights(sapply(overlaps, length), cluster_size, edges)

  # if you care about bins
  if (binning) {
    nodes = data.frame(
      id = node_ids,
      cluster_size = cluster_size,
      tightness = cluster_tightness,
      data = data_in_cluster,
      bin = get_bin_vector(binclust_data)
    )

    edges = data.frame(source = sources,
                       target = targets,
                       weight = edge_weights)

    return(list(nodes, edges))

  # if you don't
  } else {
    nodes = data.frame(
      id = node_ids,
      cluster_size = cluster_size,
      tightness = cluster_tightness,
      data = data_in_cluster
    )

    edges = data.frame(source = sources,
                       target = targets,
                       weight = edge_weights)


    return(list(nodes, edges))
  }

  # assemble graph data
  nodes = data.frame(
    id = node_ids,
    cluster_size = cluster_size,
    tightness = cluster_tightness,
    data = data_in_cluster,
    bin = bins
  )

  edges = data.frame(source = sources,
                     target = targets,
                     weight = edge_weights)


  return(list(nodes, edges))
}

cymapper <- function(mapperobject) {

  # pass to visualizer for........visualizing...
  visualize_mapper_data(mapperobject, is_ballmapper = FALSE)

  # if this isn't here R will print something useless
  return(invisible(NULL))
}


# 1D mapper ---------------------------------------------------------------

# a flavor of mapper based on projection to a single coordinate, and a width-balanced cover

create_1D_mapper_object <- function(data, dists, filtered_data, cover, clustering_method="single") {
  cover = apply(cover, 1, check_in_interval)

  return(create_mapper_object(data, dists, filtered_data, cover, clustering_method))
}

# ballmapper --------------------------------------------------------------

# TODO: add rox docs and all that jazz

construct_ballmappergraph <- function(binclust_data, dists) {
  num_vertices = max(binclust_data[[length(binclust_data)]])

  node_ids = as.character(1:num_vertices)

  overlaps = get_overlaps(binclust_data)
  edges = get_edgelist_from_overlaps(overlaps, num_vertices)
  sources = as.character(edges[, 1])
  targets = as.character(edges[, 2])

  cluster_tightness = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_size = get_cluster_sizes(binclust_data, num_vertices)
  data_in_cluster = unlist(get_clustered_data(binclust_data, num_vertices))
  edge_weights = get_edge_weights(sapply(overlaps, length), cluster_size, edges)

  nodes = data.frame(
    id = node_ids,
    cluster_size = cluster_size,
    tightness = cluster_tightness,
    data = data_in_cluster
  )

  edges = data.frame(source = sources,
                     target = targets,
                     weight = edge_weights)


  return(list(nodes, edges))
}

#' Run ballmapper
#'
#' Run mapper using an \eqn{\varepsilon}-net cover (greedily generated) and the 2D inclusion function as a filter.
#'
#' @param data Your input data. Ideally a dataframe.
#' @param dists A distance matrix for your data. Can be a `dist` object or 2D matrix.
#' @param eps A positive real number for your desired ball radius.
#' @returns A list of two dataframes, one with node data containing ball size,
#'  datapoints per ball, ball tightness, and one with edge data
#'  containing sources, targets, and weights representing overlap strength.
#' @export
get_ballmapper_data <- function(data, dists, eps) {
  balled_data = create_balls(data, dists, eps)
  formatted_balled_data = convert_balls(balled_data)

  # construct ballmapper graph
  ballmappergraph = construct_ballmappergraph(formatted_balled_data, dists)

  return(ballmappergraph)
}

cyballmapper <- function(mapperobject) {
  visualize_mapper_data(mapperobject, FALSE)
  return(invisible(NULL))
}


# clusterball mapper ------------------------------------------------------

# TODO: test that this works lol

# runner function for combo mapper; outputs bins, clusters, and the mapper graph.
get_clusterballmapper_data <- function(data, dist1, dist2, eps, method) {
  balls = create_balls(data, dist1, eps)
  clusters = get_clusters(balls, dist2, method)

  clusterballmappergraph = construct_1Dmappergraph(clusters, dist2)

  return(clusterballmappergraph)
}

cyclusterballmapper <- function(mapperobject) {
  visualize_mapper_data(mapperobject, FALSE)

  return(invisible(NULL))
}

# graph construction ------------------------------------------------------

# TODO: add rox docs

get_overlaps <- function(binclust_data) {
  num_vertices = max(binclust_data[[length(binclust_data)]]) # id of last cluster in the last bin

  flattened_data = unlist(binclust_data)
  clusters = lapply(1:num_vertices, function(x)
    flattened_data[flattened_data == x]) # sort by cluster
  cluster_names = lapply(clusters, names) # it doesn't work if you don't do this

  pairs = combn(cluster_names, 2) # get all pairs of clusters
  raw_overlaps = apply(pairs, 2, function(x)
    intersect(x[[1]], x[[2]])) # get all intersections between clusters
  names(raw_overlaps) = 1:length(raw_overlaps)
  overlaps = Filter(length, raw_overlaps) # filter out the empty intersections

  return(overlaps)
}

next_triangular <- function(x) {
  next_triangle_indx = floor((1 + sqrt(1 + 8 * x)) / 2)
  prev_triangle_val = choose(next_triangle_indx, 2)
  if (prev_triangle_val == x) {
    return (next_triangle_indx - 1)
  } else {
    return (next_triangle_indx)
  }
}

get_edgelist_from_overlaps <- function(overlaps, num_vertices) {
  overlap_names = rev(-as.numeric(names(overlaps)) + choose(num_vertices, 2) + 1)
  sources = sapply(overlap_names, function(x)
    num_vertices - next_triangular(x))
  targets = sapply(overlap_names, function(x) {
    k = next_triangular(x)
    diff = k * (k + 1) / 2 - x
    num_vertices - k + diff + 1
  })
  edges = cbind(rev(sources), rev(targets))
  return(edges)
}

