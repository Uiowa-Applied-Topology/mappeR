# expects a list of named lists with names of data and values cluster numbers
# returns the mapper graph (as an adjacency matrix)
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
        compare_cluster = flattened_data[flattened_data == j] # get the datapoints in the jth cluster
        overlap = intersect(names(my_cluster), names(compare_cluster))
        overlap.length = length(overlap)
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
