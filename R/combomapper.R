source("R/ballmapper.R")
source("R/1D_mapper.R")
library(usedist)
library(RCy3)

# performs single linkage clustering and outputs clusters based on a rough heuristic.
get_single_linkage_clusters <- function(dists) {
  hcl = hclust(as.dist(dists), method = "single")

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
get_comboclusters <- function(bins, dists, method) {
  binclust_data = list()
  cluster_count = 0

  for (i in 1:length(bins)){
    if (length(bins[[i]]) == 0) {
      binclust_data[[i]] = list()
      next
    }
    if (length(bins[[i]]) == 1) {
      binclust_data[[i]] = setNames(1, names(bins[[i]])) + cluster_count
      cluster_count = cluster_count + 1
      next
    }
    d = dists[bins[[i]],bins[[i]]]
    # d = dist_subset(as.dist(dists), names(bins[[i]])) # assuming your data has identifiers! also I am not sure of this function's behavior on non-dissimilarity matrices.
    clust = switch(tolower(method), "single" = get_single_linkage_clusters(d), stop("please tell george to do more clustering methods"))
    binclust_data[[i]] = clust + cluster_count
    cluster_count = cluster_count + max(clust) # update the total
  }
  return(binclust_data)
}

# runner function for 1D mapper; outputs bins, clusters, and the mapper graph.
get_combomapper_data <- function(data, dist1, dist2, eps) {
  print("binning...")
  balls = create_balls(data, dist1, eps)

  print("clustering...")
  binclust_data = get_comboclusters(balls, dist2, "single")

  print("making mapper graph...")
  graph_data = construct_graph(binclust_data)
  amat = graph_data[[1]]
  edge_overlaps = graph_data[[2]]
  mapper_graph = graph_from_adjacency_matrix(amat, mode="max")

  return(list(binclust_data, mapper_graph, edge_overlaps))
}

balls_to_bins <- function(balls) {

}

