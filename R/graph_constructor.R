#' @importFrom utils combn
get_overlaps <- function(binclust_data) {
  num_vertices = max(binclust_data[[length(binclust_data)]]) # id of last cluster in the last bin

  flattened_data = unlist(binclust_data)
  clusters = lapply(1:num_vertices, function(x) flattened_data[flattened_data == x]) # sort by cluster
  cluster_names = lapply(clusters, names) # it doesn't work if you don't do this

  pairs = combn(cluster_names, 2) # get all pairs of clusters
  raw_overlaps = apply(pairs, 2, function(x) intersect(x[[1]], x[[2]])) # get all intersections between clusters
  names(raw_overlaps) = 1:length(raw_overlaps)
  overlaps = Filter(length, raw_overlaps) # filter out the empty intersections

  return(overlaps)
}

next_triangular <- function(x) {
  next_triangle_indx = floor((1 + sqrt(1 + 8*x))/2)
  prev_triangle_val = choose(next_triangle_indx, 2)
  if (prev_triangle_val == x) {
    return (next_triangle_indx - 1)
  } else {
    return (next_triangle_indx)
  }
}

get_edgelist_from_overlaps <- function(overlaps, num_vertices) {
  overlap_names = rev(-as.numeric(names(overlaps)) + choose(num_vertices, 2) + 1)
  sources = sapply(overlap_names, function(x) num_vertices - next_triangular(x))
  targets = sapply(overlap_names, function(x) {
    k = next_triangular(x)
    diff = k*(k+1)/2 - x
    num_vertices - k + diff + 1
  })
  edges = cbind(rev(sources), rev(targets))
  return(edges)
}
