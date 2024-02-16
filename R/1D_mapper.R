source("R/graph_constructor.R")
source("R/mapper_viz.R")

# given an interval, calculates endpoints of a fixed number of evenly spaced, equal length, overlapping subintervals with a fixed percent overlap.
get_width_balanced_endpoints <- function(min_val, max_val, num_bins, percent_overlap) {

  even_length = (max_val - min_val)/num_bins # widths with zero percent overlap
  nudge = (percent_overlap/100)*even_length

  left_ends = min_val + (even_length-nudge)*(0:(num_bins-1)) # construct correctly overlapping bins
  right_ends = left_ends + even_length # we will scale everything after

  scale_factor = (max_val - min_val)/(right_ends[num_bins] - min_val) # scale by pretending min_val = 0
  bin_ends = cbind(left_ends, right_ends) # make bins

  bin_ends = scale_factor*(bin_ends - min_val) + min_val # translate to zero, scale, then translate back

  return(bin_ends)
}

# makes a list of subsets of the input data, according to specified "bins."
# assumes a real-valued filter function.
make_bins <- function(data, filtered_data, bin_ends) {
  bins = list()
  num_bins = nrow(bin_ends)

  for (i in 1:num_bins) {
    # left and right bin endpoints
    bin_left = as.numeric(bin_ends[i,1])
    bin_right = as.numeric(bin_ends[i,2])

    # check which filtered datapoints fall between the bin endpoints
    in_bin = sapply(filtered_data, function(x) (bin_left - x <= 0) & (bin_right - x >= 0))
    bin_assignments = which(in_bin)

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
  branch_heights = diff(heights)

  tallest_branch_height = max(branch_heights)
  tall_branch_id = which(branch_lengths == tallest_branch_height)

  cutval = (tall_branch_height + heights[tall_branch_id + 1])/2

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
    if (length(bins[[i]]) == 0) { # nothing in this bin, moving on
      binclust_data[[i]] = list()
      next
    }
    if (nrow(bins[[i]]) == 1) { # I am atoning for the sins of TDAmapper
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
  bins = get_width_balanced_endpoints(min(filtered_data), max(filtered_data), num_bins, percent_overlap)
  binned_data = make_bins(data, filtered_data, bin_endpoints)

  # cluster data
  print("clustering...")
  binclust_data = get_clusters(binned_data, dists, clustering_method)

  # construct mapper graph
  print("making mapper graph...")
  graph_data = construct_graph(binclust_data)
  amat = graph_data[[1]]
  edge_overlaps = graph_data[[2]]
  mapper_graph = graph_from_adjacency_matrix(amat, mode="max")

  # return the binned and clustered data, the mapper graph, and the edge overlap data
  return(list(binclust_data, mapper_graph, edge_overlaps, bins))
}

cymapper <- function(data, filtered_data, dists, num_bins, percent_overlap, clustering_method) {
  # generate mapper data
  mapper_data = get_mapper_data(data, filtered_data, dists, num_bins, percent_overlap, clustering_method)

  # pass to visualizer for........visualizing...
  visualize_mapper_data(mapper_data, dists, is_ballmapper = FALSE)

  # if this isn't here R will print something useless
  return(invisible(NULL))
}
