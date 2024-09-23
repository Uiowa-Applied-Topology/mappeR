# dendrogram processing ---------------------------------------------------

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
