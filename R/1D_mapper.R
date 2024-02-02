source("R/graph_constructor.R")

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

  if (length(heights) == 1) {
    cutval = heights[1]
  } else {
    # idea is the most persistent number of clusters is the best one
    # this number is the tallest section of the dendrogram without merges
    for (i in 1:(length(heights)-1)){
      currentdiff = abs(heights[i]-heights[i+1])
      if (currentdiff > maxdiff) {
        maxdiff = currentdiff
        cutval = heights[i]
      }
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
# output is a list of named vectors; there is one named vector of data per bin containing a cluster number.
get_clusters <- function(bins, dists, method) {
  binclust_data = list()
  cluster_count = 0

  for (i in 1:length(bins)){
    if (length(bins[[i]]) == 0) {
      binclust_data[[i]] = list()
      next
    }
    if (nrow(bins[[i]]) == 1) {
      binclust_data[[i]] = setNames(1, rownames(bins[[i]])) + cluster_count
      cluster_count = cluster_count + 1
      next
    }
    d = dist_subset(dists, rownames(bins[[i]])) # assuming your data has identifiers! also I am not sure of this function's behavior on non-dissimilarity matrices.
    clust = switch(tolower(method), "single" = get_single_linkage_clusters(d), stop("please tell george to do more clustering methods"))
    binclust_data[[i]] = clust + cluster_count
    cluster_count = cluster_count + max(clust) # update the total
  }
  return(binclust_data)
}

# runner function for 1D mapper; outputs bins, clusters, and the mapper graph.
get_mapper_data <- function(data, filtered_data, dists, num_bins, percent_overlap, clustering_method) {
  # bin data according to filter values
  print("binning...")
  binned_data = make_bins(data, filtered_data, get_width_balanced_endpoints(min(filtered_data), max(filtered_data), num_bins, percent_overlap))

  # cluster data
  print("clustering...")
  binclust_data = get_clusters(binned_data, dists, clustering_method)

  # construct mapper graph
  print("making mapper graph...")
  graph_data = construct_graph(binclust_data)
  amat = graph_data[[1]]
  edge_overlaps = graph_data[[2]]
  mapper_graph = graph_from_adjacency_matrix(amat, mode="max")

  return(list(binclust_data, mapper_graph, edge_overlaps))
}
