###########################################################################
# CLUSTER FACTORY
# clustering
###########################################################################

# goblin clustering mines -------------------------------------------------

#' Ship data off to the clustering goblins
#'
#' This function tells the computer to look away for a second, so the goblins come and cluster your data while it's not watching.
#'
#' @param dist_mats A list of distance matrices of each bin that is to be clustered.
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix in a list.
#'
#' @return The output of `clusterer(dist_mats)`, which needs to be a list containing named vectors (one per bin), whose names are data point names and whose values are cluster labels (within each bin)
get_raw_clusters <- function(dist_mats, clusterer) {
  return(clusterer(dist_mats))
}



#' Perform the clustering step in mapper
#'
#' This function processes the binned data and global distance matrix to return freshly clustered data.
#'
#' @param bins A list containing "bins" of vectors of names of data points.
#' @param dists A distance matrix containing pairwise distances between named data points.
#' @param clusterer A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix.
#'
#' @return The output of `clusterer` applied to a list of distance matrices, which should be a list containing named vectors (one per bin), whose names are data point names and whose values are cluster labels.
get_clusters <- function(bins, dists, clusterer) {
  # more than one bin, need more than one distance matrix
  if (is.list(bins)) {
    # subset the global distance matrix per bin
    dist_mats = mapply(subset_dists, bins, MoreArgs = list(dists = dists), SIMPLIFY = FALSE)

    # cluster the data
    clusters = get_raw_clusters(dist_mats, clusterer)

    # accurately total up clusters
    clusters_per_bin = sapply(clusters, max)
    offset = c(0, cumsum(clusters_per_bin))
    clusters = mapply(function(x, y)
      x + y, clusters, offset[-length(offset)], SIMPLIFY = FALSE)
    return(clusters)
  }
#
  # cluster the data
  clusters = get_raw_clusters(subset_dists(bins, dists), clusterer) # this fixed everything????

  return(clusters)
}

