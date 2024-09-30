# dendrogram processing ---------------------------------------------------

# get tallest branch of a single dendrogram
get_tallest_branch <- function(dend) {
  heights = sort(unique(cophenetic(dend)))
  if (length(heights) <= 1) {
    return(max(heights))
  }
  branch_lengths = diff(heights)
  return(max(branch_lengths))
}

# cuts...a dendrogram...
cut_dendrogram <- function(dend, threshold) {
  # TODO remove all the duplicate code lol
  heights = sort(unique(cophenetic(dend))) # merge heights of dendrogram

  if (length(heights) <= 1) {
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
  thresholdcondition = tallest_branch_height < threshold
  indexofdispersion = sd(heights) ^ 2 / mean(heights)
  dispersioncondition = indexofdispersion < .015

  # plot(dend, xlab=paste("index of dispersion: ", round(indexofdispersion, 3), " too low? ", dispersioncondition, ", tallest branch: ", round(tallest_branch_height, 3), ", too short? ", thresholdcondition))

  if (thresholdcondition | dispersioncondition) {
    return(cutree(dend, k = 1))
  } else {
    return(cutree(dend, h = cutval))
  }
}

# given a bunch of dendrograms, calculates a threshold a dendrogram's merge heights must clear to represent more than one cluster
process_dendrograms <- function(dends) {
  tallest_branches = sapply(dends, get_tallest_branch)
  biggest_branch_length = max(tallest_branches)
  threshold = biggest_branch_length * .1

  snipped_dends = mapply(cut_dendrogram,
                         dend = dends,
                         MoreArgs = list(threshold = threshold))

  return(snipped_dends)
}
