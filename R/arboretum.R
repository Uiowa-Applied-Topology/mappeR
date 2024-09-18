#' # performs single linkage clustering and outputs clusters based on a rough heuristic.
#' #' @importFrom stats as.dist
#' #' @importFrom stats hclust
#' #' @importFrom stats cophenetic
#' #' @importFrom stats cutree
#' get_single_linkage_clusters <- function(dists) {
#'   dists = as.dist(dists)
#'   hcl = hclust(dists, method = "single")
#'
#'   heights = sort(unique(cophenetic(hcl))) # merge heights of dendrogram
#'   branch_lengths = diff(heights) # differences are branch lengths
#'
#'   tallest_branch_height = max(branch_lengths)
#'   tallest_branch_id = which(branch_lengths == tallest_branch_height)
#'
#'   cutval = (tallest_branch_height + heights[tallest_branch_id + 1])/2 # cut at midpoint of tallest branch
#'
#'   # if the dendrogram is "very short" we might expect the best number of clusters is one
#'   if (max(heights) < 0.1) { # I made this number up
#'     return(cutree(hcl, k=1))
#'   } else {
#'     return(cutree(hcl, h=cutval))
#'   }
#' }

# get tallest branch of a single dendrogram
get_tallest_branch <- function(dend) {
  heights = sort(unique(cophenetic(dend)))
  branch_lengths = diff(heights)
  return(max(branch_lengths))
}

# cuts...a dendrogram...
cut_dendrogram <- function(dend, threshold) {
  # TODO remove all the duplicate code lol

  heights = sort(unique(cophenetic(dend))) # merge heights of dendrogram
  branch_lengths = diff(heights) # differences are branch lengths

  tallest_branch_height = max(branch_lengths)
  tallest_branch_id = which(branch_lengths == tallest_branch_height)

  cutval = (tallest_branch_height + heights[tallest_branch_id + 1])/2 # midpoint of tallest branch

  if (tallest_branch_height < threshold) {
    return(cutree(dend, k=1))
  } else {
    return(cutree(dend, h=cutval))
  }
}

# given a bunch of dendrograms, calculates a threshold a dendrogram's merge heights must clear to represent more than one cluster
process_dendrograms <- function(dends) {
  tallest_branches = sapply(dends, get_tallest_branch)
  biggest_branch_length = max(tallest_branches)
  threshold = biggest_branch_length*.05

  snipped_dends = mapply(cut_dendrogram, dend=dends, MoreArgs = list(threshold=threshold))

  return(snipped_dends)
}

# please don't ask
run_slink <- function(dist) {
  return(hclust(dist, "single"))
}

# performs single linkage clustering and outputs clusters based on a rough heuristic.
#' @importFrom stats as.dist
#' @importFrom stats hclust
#' @importFrom stats cophenetic
#' @importFrom stats cutree
get_single_linkage_clusters <- function(dist_mats) {
  dends = lapply(dist_mats, run_slink)
  return(process_dendrograms(dends))
}
