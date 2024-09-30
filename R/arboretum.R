# dendrogram processing ---------------------------------------------------

#' Find the tallest branch of a dendrogram
#'
#' @param dend A single dendrogram.
#'
#' @return The height of the tallest branch (longest time between merge heights) of the input dendrogram.
#' @export
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
#' @return A named vector whose names are data point names and whose values are cluster labels. The number of clusters is determined to be 1 if the tallest branch of the dendrogram is less than the threshold, or if the index of dispersion (standard deviation squared divided by mean) of the branch heights is too low. Otherwise, we cut at the longest branch of the dendrogram to determine the number of clusters.
#' @export
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

  # uncomment this to plot the dendrograms that come through here with their stats
  # plot(dend, xlab=paste("index of dispersion: ", round(indexofdispersion, 3), " too low? ", dispersioncondition, ", tallest branch: ", round(tallest_branch_height, 3), ", too short? ", thresholdcondition))

  if (thresholdcondition | dispersioncondition) {
    return(cutree(dend, k = 1))
  } else {
    return(cutree(dend, h = cutval))
  }
}

#' Cut many dendrograms
#'
#' @param dends A list of dendrograms to be cut.
#'
#' @return A list of named vectors (one per dendrogram) whose names are data point names and whose values are cluster labels. This function determines a global minimum threshold based on the longest branches in all the input dendrograms, and uses that as a heuristic to gauge if the best number of clusters is 1, or the value obtained by cutting the longest branch.
#' @export
process_dendrograms <- function(dends) {
  tallest_branches = sapply(dends, get_tallest_branch)
  biggest_branch_length = max(tallest_branches)
  threshold = biggest_branch_length * .1

  snipped_dends = mapply(cut_dendrogram,
                         dend = dends,
                         MoreArgs = list(threshold = threshold))

  return(snipped_dends)
}
