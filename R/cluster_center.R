#' Ship data off to the clustering goblins
#'
#' This function tells the computer to look away for a second, so the goblins come and cluster your data while it's not watching.
#'
#' @param dist_mats A list of distance matrices of each bin that is to be clustered.
#' @param method A string that suggests how the goblins will handle the data.
#'
#' @return A list containing named vectors (one per bin), whose names are datapoint names and whose values are cluster labels (within each bin)
run_cluster_machine <- function(dist_mats, method) {
  switch(tolower(method), "single" = return(get_single_linkage_clusters(dist_mats)))
}

#' Initate the clustering process
#'
#' This function processes the binned data and global distance matrix to return freshly clustered data.
#'
#' @param bins A list containing "bins" of vectors of names of datapoints.
#' @param dists A distance matrix containing pairwise distances between named datapoints.
#' @param method A clue!
#'
#' @return A list containing named vectors (one per bin), whose names are datapoint names and whose values are cluster labels
#'
#' @importFrom stats as.dist
get_clusters <- function(bins, dists, method) {
  # subset the global distance matrix per bin
  dist_mats = sapply(1:length(bins), function(x)
    as.dist(as.matrix(dists)[bins[[x]], bins[[x]]]))

  # cluster the data
  clusters = run_cluster_machine(dist_mats, method)

  # accurately total up clusters
  clusters_per_bin = sapply(clusters, max)
  offset = c(0, cumsum(clusters_per_bin))
  clusters = mapply(function(x, y)
    x + y, clusters, offset[-length(offset)])
  return(clusters)
}
