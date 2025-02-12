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
#' @param method A string to pass to `fastcluster` to determine clustering method.
#' @param global_clustering Whether you want clustering to happen in a global (all level visible) or local (only current level set visible) context
#'
#' @return A list containing named vectors (one per bin), whose names are data point names and whose values are cluster labels (within each bin)
run_cluster_machine <- function(dist_mats, method, global_clustering = TRUE) {
  if (method %in% c("single", "complete", "average", "mcquitty", "centroid", "median", "ward.D", "ward.D2")) {
    return(get_hierarchical_clusters(dist_mats, method, global_clustering))
  } else {
    stop("not a valid clustering method")
  }
}

#' Subset a distance matrix
#'
#' @param bin A list of names of data points.
#' @param dists A distance matrix for data points in the bin, possibly including extra points.
#'
#' @return A distance matrix for only the data points in the input bin.
subset_dists <- function(bin, dists) {
  bin_length = length(bin)
  if (bin_length == 0) {
    return(NA)
  } else if (bin_length == 1) {
    return(bin)
  } else {
    res = as.dist(as.matrix(dists)[bin, bin]) # this is how it's done in the usedist package
    return(res)
  }
}

#' Initate the clustering process
#'
#' This function processes the binned data and global distance matrix to return freshly clustered data.
#'
#' @param bins A list containing "bins" of vectors of names of data points.
#' @param dists A distance matrix containing pairwise distances between named data points.
#' @param method A string to pass to [hclust] to determine clustering method.
#' @param global_clustering Whether you want clustering to happen in a global (all level visible) or local (only current level set visible) context
#'
#' @return A list containing named vectors (one per bin), whose names are data point names and whose values are cluster labels
get_clusters <- function(bins, dists, method, global_clustering = TRUE) {
  # more than one bin, need more than one distance matrix
  if (is.list(bins)) {
    # subset the global distance matrix per bin
    dist_mats = mapply(subset_dists, bins, MoreArgs = list(dists = dists), SIMPLIFY = FALSE)

    # cluster the data
    clusters = run_cluster_machine(dist_mats, method, global_clustering)

    # accurately total up clusters
    clusters_per_bin = sapply(clusters, max)
    offset = c(0, cumsum(clusters_per_bin))
    clusters = mapply(function(x, y)
      x + y, clusters, offset[-length(offset)], SIMPLIFY = FALSE)
    return(clusters)
  }
#
  # cluster the data
  clusters = run_cluster_machine(subset_dists(bins, dists), method, global_clustering) # this fixed everything????

  return(clusters)
}

# caveman clustering ------------------------------------------------------

#' The easiest clustering method
#'
#' @param bins A list of bins, each containing names of data from some data frame.
#'
#' @return A named vector whose names are data point names and whose values are cluster labels
convert_to_clusters <- function(bins) {
  ball_sizes = lapply(bins, length)

  # repeat the cluster id for as many data points belonging to that bin
  ballball_data = unlist(mapply(function(x, y)
    rep(x, y), 1:length(ball_sizes), ball_sizes))

  # make sure names match up
  names(ballball_data) = unlist(bins)

  return(ballball_data)
}


# hierarchical clustering -------------------------------------------------


#' Perform hierarchical clustering and process dendrograms
#'
#' @param dist_mats A list of distance matrices to be used for clustering.
#' @param method A string to pass to [hclust] to determine clustering method.
#' @param global_clustering Whether you want clustering to happen in a global (all level visible) or local (only current level set visible) context
#'
#' @return A list containing named vectors (one per dendrogram), whose names are data point names and whose values are cluster labels
get_hierarchical_clusters <- function(dist_mats, method, global_clustering = TRUE) {
  dends = lapply(dist_mats, run_link, method = method)
  real_dends = dends[lapply(dends, length) > 1]
  imposter_dends = dends[lapply(dends, length) == 1]
  processed_dends = process_dendrograms(real_dends, global_clustering)
  if (length(imposter_dends) != 0) {
    return(append(processed_dends, sapply(imposter_dends, function(x)
      list(unlist(x))))) # LMAO what is this
  } else {
    return(processed_dends)
  }
}

#' Perform single linkage clustering
#'
#' @param dist A distance matrix.
#' @param method A string to pass to [hclust] to determine clustering method.
#'
#' @return A dendrogram generated by `fastcluster`.
run_link <- function(dist, method) {
  if (!(inherits(dist, "dist")) & (any(is.na(dist)))) {
    return(vector())
  } else if (!(inherits(dist, "dist"))) {
    res = list(1)
    names(res) = dist
    return(res)
  } else {
    return(fastcluster::hclust(dist, method))
  }
}

# dendrogram processing ---------------------------------------------------

#' Find the tallest branch of a dendrogram
#'
#' @param dend A single dendrogram.
#'
#' @return The height of the tallest branch (longest time between merge heights) of the input dendrogram.
get_tallest_branch <- function(dend) {
  heights = sort(unique(cophenetic(dend)))
  if (length(heights) <= 1) {
    return(max(heights))
  }
  branch_lengths = diff(heights)
  return(max(branch_lengths))
}

#' Cut a dendrogram
#'
#' @param dend A single dendrogram.
#' @param threshold A mininum tallest branch value.
#'
#' @return A named vector whose names are data point names and whose values are cluster labels.
#' @details The number of clusters is determined to be 1 if the tallest branch of the dendrogram is less than the threshold, or if the index of dispersion (standard deviation squared divided by mean) of the branch heights is below 0.015. Otherwise, we cut at the longest branch of the dendrogram to determine the number of clusters.
cut_dendrogram <- function(dend, threshold) {
  # TODO remove all the duplicate code lol
  heights = sort(unique(cophenetic(dend))) # merge heights of dendrogram

  if (length(heights) <= 2) {
    if (max(heights) < threshold) {
      return(cutree(dend, k = 1))
    } else {
      return(cutree(dend, k = 2))
    }
  }

  branch_lengths = diff(heights) # differences are branch lengths

  tallest_branch_height = max(branch_lengths)
  tallest_branch_id = which(branch_lengths == tallest_branch_height)

  cutval = (tallest_branch_height + heights[tallest_branch_id + 1]) / 2 # midpoint of tallest branch
  if (length(cutval) > 1) {
    cutval = sample(cutval, 1)
  }

  # one cluster condition: dendrogram has no sufficiently tall branches
  thresholdcondition = tallest_branch_height < threshold

  # one cluster condition: lengths of branches are not well-dispersed
  indexofdispersion = sd(heights) ^ 2 / mean(heights)
  dispersioncondition = indexofdispersion < .015

  # uncomment this to plot the dendrograms that come through here with their stats
  plot(dend, xlab=paste("index of dispersion: ", round(indexofdispersion, 3), " lower than 0.015? ", dispersioncondition, ", tallest branch: ", round(tallest_branch_height, 3), "threshold: ", threshold, ", lower than threshold? ", thresholdcondition))
  # add horiz (well, vertical) line:
  if (thresholdcondition | dispersioncondition) {
    abline(h = max(heights), lty = 2)
    return(cutree(dend, k = 1))
  } else {
    abline(h = cutval, lty = 2)
    return(cutree(dend, h = cutval))
  }
}

#' Cut many dendrograms
#'
#' @param dends A list of dendrograms to be cut.
#'
#' @return A list of named vectors (one per dendrogram) whose names are data point names and whose values are cluster labels.
#' @details This function uses a value of 10 percent of the tallest branch across dendrograms as a threshold for [cut_dendrogram].
#' @param global_clustering Whether you want clustering to happen in a global (all level visible) or local (only current level set visible) context.
process_dendrograms <- function(dends, global_clustering = TRUE) {
  if (inherits(dends, "hclust")) {
    return(cut_dendrogram(dends, 0))
  }

  tallest_branches = sapply(dends, get_tallest_branch)
  biggest_branch_length = max(tallest_branches)
  threshold = ifelse(global_clustering, biggest_branch_length * .1, 0)

  snipped_dends = mapply(cut_dendrogram,
                         dend = dends, SIMPLIFY = FALSE,
                         MoreArgs = list(threshold = threshold))
  return(snipped_dends)
}
