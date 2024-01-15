# given an interval, calculates endpoints of a fixed number of evenly spaced, equal length, overlapping subintervals with a fixed percent overlap.
get_width_balanced_endpoints <- function(min_val, max_val, num_bins, percent_overlap) {

  even_length = (max_val - min_val)/num_bins # widths with zero percent overlap
  nudge = (percent_overlap/100)*even_length

  left_ends = c(min_val) # first left endpoint is unchanged
  right_ends = c(min_val + even_length) # first right endpoint is unchanged

  for (j in 1:(num_bins-1)) {
    leftcut = min_val + j*(even_length - nudge) # move the left endpoint so it overlaps with its neighbor appropriately
    rightcut = leftcut + even_length # don't change the length of the segment yet
    left_ends = append(left_ends, leftcut)
    right_ends = append(right_ends, rightcut)
  }

  scale_factor = (max_val - min_val)/(right_ends[num_bins] - min_val) # scale by pretending min_val = 0

  bin_ends = cbind(left_ends, right_ends) # make bins

  bin_ends = scale_factor*(bin_ends - min_val) + min_val # scale the bins after translating, then translate back

  return(bin_ends)
}

# makes a list of subsets of the input data, according to specified "bins."
make_bins <- function(data, filtered_data, bin_ends) {
  bins = list()
  num_bins = nrow(bin_ends)

  for (i in 1:num_bins) {
    bin_assignments = c()
    bin_left = as.numeric(bin_ends[i,1])
    bin_right = as.numeric(bin_ends[i,2])
    for (j in 1:length(filtered_data)) { # loop through filtered data so we can look at level sets
      if ((bin_left - filtered_data[j] <= 0) & (bin_right - filtered_data[j] >= 0)) {
        bin_assignments = append(bin_assignments, j) # don't need the whole datapoint, just the index
      }
    }
    if (length(bin_assignments) != 0) { # this means our bin is empty
      bins[[i]] = data[bin_assignments,] # get a subset of the original data based on the indices we collected
    } else {
      bins[[i]] = list() # bin still exists, it's just empty
    }
  }
  return(bins)
}

# performs single linkage clustering and outputs clusters based on a rough heuristic.
get_single_linkage_clusters <- function(dists) {
  hcl = hclust(dists, method = "single")

  cutval = 0
  maxdiff = 0
  heights = sort(unique(cophenetic(hcl))) # merge heights of dendrogram

  # idea is the most persistent number of clusters is the best one
  # this number is the tallest section of the dendrogram without merges
  for (i in 1:(length(heights)-1)){
    currentdiff = abs(heights[i]-heights[i+1])
    if (currentdiff > maxdiff) {
      maxdiff = currentdiff
      cutval = heights[i]
    }
  }

  # if the dendrogram is "very short" we might expect the best number of clusters is one
  if (max(heights) < 0.1) { # I made this number up
    return(cutree(hcl, k=1))
  } else {
    return(cutree(hcl, h=cutval))
  }
}

# given binned data and the full data's distance matrix, clusters within each bin. keeps track of total clusters across bins.
get_clusters <- function(bins, dists, method) {
  clustered_data = list()
  cluster_count = 0

  for (i in 1:length(bins)){
    d = dist_subset(dists, rownames(bins[[i]])) # assuming your data has identifiers! also I am not sure of this function's behavior on non-dissimilarity matrices.
    clust = switch(tolower(method), "single" = get_single_linkage_clusters(d), stop("please tell george to do more clustering methods"))
    clustered_data[[i]] = clust + cluster_count
    cluster_count = cluster_count + max(clust) # update the total
  }
  return(clustered_data)
}

get_bin_vector <- function(clustered_data) {
  clusters_and_bins = c()
  for (i in 1:length(clustered_data)) {
    for (j in unique(clustered_data[[i]])) {
      clusters_and_bins = append(clusters_and_bins, i)
    }
  }
  print(clusters_and_bins)
  return(clusters_and_bins)
}

# constructs an abstract graph with clusters as vertices and nonempty intersections between clusters as edges.
construct_graph <- function(clustered_data) {
  num_vertices = max(clustered_data[[length(clustered_data)]]) # I don't know why this works

  flattened_data = unlist(clustered_data)

  amat = matrix(, nrow = num_vertices, ncol = num_vertices)
  for (i in 1:(num_vertices-1)) {
    for (j in i:num_vertices) {
      if (i == j) {
        amat[i, j] = 0
      } else {
        my_cluster = flattened_data[flattened_data == i] # get the datapoints in the ith cluster
        compare_cluster = flattened_data[flattened_data == j] # get the datapoints in the jth cluster
        if (length((intersect(names(my_cluster), names(compare_cluster)))) != 0) {
          amat[i, j] = 1
        } else {
          amat[i, j] = 0
        }
      }
    }
  }
  return(amat)
}

# runner function for 1D mapper; outputs bins, clusters, and the mapper graph.
get_mapper_data <- function(data, filtered_data, dists, num_bins, percent_overlap, clustering_method) {
  # bin data according to filter values
  binned_data = make_bins(data, filtered_data, get_width_balanced_endpoints(min(filtered_data), max(filtered_data), num_bins, percent_overlap))

  # cluster data
  clustered_data = get_clusters(binned_data, dists, clustering_method)

  # construct mapper graph
  amat = construct_graph(clustered_data)

  return(list(clustered_data, amat))
}

visualize_mapper_data <- function(mapper_data) {
  clustered_data = mapper_data[[1]]
  amat = mapper_data[[2]]

  num_vertices = nrow(amat)
  bin_vector = get_bin_vector(clustered_data)

  sourcevec = c()
  targetvec = c()

  for (i in 1:(num_vertices - 1)) {
    for (j in i:num_vertices) {
      if (amat[i, j] == 1) {
        sourcevec = append(sourcevec, c(i,j))
        targetvec = append(targetvec, c(j,i))
      }
    }
  }

  print(sourcevec)
  print(targetvec)

  nodes <- data.frame(id = 1:num_vertices,
                      group = bin_vector)
  edges <- data.frame(source = sourcevec,
                      target = targetvec)

  createNetworkFromDataFrames(nodes, edges, title = "mapper graph", collection = "example")
}

circle.data = data.frame( x= sapply(1:1000, function(x) cos(x)) + rnorm(100, 500, .03), y = sapply(1:1000, function(x) sin(x)) + rnorm(100, 0, 0.03))
circle.dist = dist(circle.data)
visualize_mapper_data(get_mapper_data(circle.data, circle.data$x, circle.dist, 5, 10, "single"))






