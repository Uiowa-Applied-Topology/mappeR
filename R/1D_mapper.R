#' Create a single bin of data
#'
#' This function identifies datapoints whose "filter" value falls within a numeric range.
#'
#' @param bin_left The minimum value to be filtered for.
#' @param bin_right The maximum value to be filtered for.
#' @param data A dataframe.
#' @param filtered_data A single column of the input dataframe.
#'
#' @return A vector of datapoint names which fall within the filter range.
make_one_bin <- function(bin_left, bin_right, data, filtered_data) {
  in_bin = sapply(filtered_data, function(x)
    (bin_left - x <= 0) & (bin_right - x >= 0))
  bin_assignments = which(in_bin) # tells us indices of binned datapoints
  if (length(bin_assignments) != 0) {
    return(rownames(data[bin_assignments, ])) # why do we need to take the name from the original data? weird R thing that I am too unbothered to address feel free to annoy me if you want
  } else {
    return(vector()) # bin still exists, it's just empty
  }
}


#' Create bins of data from a 1D filter function
#'
#' This function identifies datapoints whose "filter" value falls within a specified set of numeric ranges.
#'
#' @param data A dataframe.
#' @param filtered_data A single column of the input dataframe.
#' @param bin_ends A 2D numeric array of input ranges, whose first column contains left endpoints and whose second column contains right endpoints.
#'
#' @return A list of "bins" containing vectors of names of datapoints.
make_bins <- function(data, filtered_data, bin_ends) {
  left_ends = bin_ends[, 1]
  right_ends = bin_ends[, 2]
  bins = mapply(
    make_one_bin,
    left_ends,
    right_ends,
    MoreArgs = list(data = data, filtered_data = filtered_data)
  )
  return(bins)
}


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

#' Run 1D mapper algorithm!
#'
#' @param binclust_data A list containing named vectors whose names are datapoint names and whose values are cluster labels.
#' @param dists A distance matrix containing pairwise distances between named datapoints.
#'
#' @return A list containing:
#' \itemize{
#'    \item nodes - A dataframe containing vertex/cluster information:
#'      \itemize{
#'          \item id - The vertex label.
#'          \item size - How many datapoints are contained in that vertex.
#'          \item tightness - A measure of dispersion of the data inside the vertex.
#'          \item data - A string containing names of the datapoints within the vertex.
#'          \item bin - Which bin the vertex belongs to.
#'      }
#'    \item edges - A dataframe containing vertex/cluster overlap information:
#'      \itemize{
#'          \item source - One endpoint of an edge.
#'          \item target - The other endpoint of an edge.
#'          \item weight - A measure of the significance of the edge.
#'      }
#' }
construct_1Dmappergraph <- function(binclust_data, dists) {
  # grab basic graph characteristics
  num_vertices = max(binclust_data[[length(binclust_data)]])
  node_ids = as.character(1:num_vertices)
  overlaps = get_overlaps(binclust_data)
  edges = get_edgelist_from_overlaps(overlaps, num_vertices)
  sources = as.character(edges[, 1])
  targets = as.character(edges[, 2])

  # calculate some cluster stats
  cluster_tightness = get_cluster_tightness_vector(as.matrix(dists), binclust_data, num_vertices)
  cluster_size = get_cluster_sizes(binclust_data, num_vertices)
  data_in_cluster = unlist(get_clustered_data(binclust_data, num_vertices))
  edge_weights = get_edge_weights(sapply(overlaps, length), cluster_size, edges)
  bins = get_bin_vector(binclust_data)

  # assemble graph data
  nodes = data.frame(
    id = node_ids,
    size = cluster_size,
    tightness = cluster_tightness,
    data = data_in_cluster,
    bin = bins
  )

  edges = data.frame(source = sources,
                     target = targets,
                     weight = edge_weights)


  return(list(nodes, edges))
}

#' Runs 1D mapper returns a dataframe with node and edge data.
#'
#' @param data A dataframe.
#' @param filtered_data A single column of the input data.
#' @param dists A distance matrix containing pairwise relations among the input data. Can be a `dist` object or 2D matrix.
#' @param num_bins The number of "bins" to split the input data into based on the filter. A positive integer.
#' @param percent_overlap The percent overlap desired between each "bin." An integer between 0 and 100 (inclusive).
#' @param clustering_method Desired clustering method. A string from these options: "single" (single-linkage hierarchical)
#'
#' @return A list containing:
#' \itemize{
#'    \item nodes - A dataframe containing vertex/cluster information:
#'      \itemize{
#'          \item id - The vertex label.
#'          \item size - How many datapoints are contained in that vertex.
#'          \item tightness - A measure of dispersion of the data inside the vertex.
#'          \item data - A string containing names of the datapoints within the vertex.
#'          \item bin - Which bin the vertex belongs to.
#'      }
#'    \item edges - A dataframe containing vertex/cluster overlap information:
#'      \itemize{
#'          \item source - One endpoint of an edge.
#'          \item target - The other endpoint of an edge.
#'          \item weight - A measure of the significance of the edge.
#'      }
#' }
#' @export
get_1D_mapper_data <- function(data,
                               filtered_data,
                               dists,
                               num_bins,
                               percent_overlap,
                               clustering_method) {
  # bin data according to filter values
  bins = create_width_balanced_cover(min(filtered_data),
                                     max(filtered_data),
                                     num_bins,
                                     percent_overlap)
  binned_data = make_bins(data, filtered_data, bins)

  # cluster data
  binclust_data = get_clusters(binned_data, dists, clustering_method)

  # construct mapper graph
  mappergraph = construct_1Dmappergraph(binclust_data, dists)

  return(mappergraph)
}

#' Runs 1D mapper and passes data to Cytoscape for visualization.
#'
#' @param data A dataframe.
#' @param filtered_data A single column of the input data.
#' @param dists A distance matrix containing pairwise relations among the input data. Can be a `dist` object or 2D matrix.
#' @param num_bins The number of "bins" to split the input data into based on the filter. A positive integer.
#' @param percent_overlap The percent overlap desired between each "bin." An integer between 0 and 100 (inclusive).
#' @param clustering_method Desired clustering method. A string from these options: "single" (single-linkage hierarchical)
#'
#' @returns NULL
#'
#' @export
cymapper <- function(data,
                     filtered_data,
                     dists,
                     num_bins,
                     percent_overlap,
                     clustering_method) {
  # generate mapper data
  mappergraph = get_1D_mapper_data(data,
                                   filtered_data,
                                   dists,
                                   num_bins,
                                   percent_overlap,
                                   clustering_method)

  # pass to visualizer for........visualizing...
  visualize_mapper_data(mappergraph, is_ballmapper = FALSE)

  # if this isn't here R will print something useless
  return(invisible(NULL))
}
