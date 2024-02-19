# performs single linkage clustering and outputs clusters based on a rough heuristic.
#' @importFrom stats as.dist
#' @importFrom stats hclust
#' @importFrom stats cophenetic
#' @importFrom stats cutree
get_single_linkage_clusters <- function(dists) {
  dists = as.dist(dists)
  hcl = hclust(dists, method = "single")

  heights = sort(unique(cophenetic(hcl))) # merge heights of dendrogram
  branch_lengths = diff(heights) # differences are branch lengths

  tallest_branch_height = max(branch_lengths)
  tallest_branch_id = which(branch_lengths == tallest_branch_height)

  cutval = (tallest_branch_height + heights[tallest_branch_id + 1])/2 # cut at midpoint of tallest branch

  # if the dendrogram is "very short" we might expect the best number of clusters is one
  if (max(heights) < 0.1) { # I made this number up
    return(cutree(hcl, k=1))
  } else {
    return(cutree(hcl, h=cutval))
  }
}
