source("R/ballmapper_viz.R")
source("R/combomapper_viz.R")

edit_dist_data = as.matrix(read.csv("edit_dist_edges.csv"))
tree_dist_data = as.matrix(read.csv("tree_dist_edges.csv"))
node_data = as.matrix(read.csv("aptamer-data-raw.csv"))

num_datapoints = nrow(node_data)

edit_amat = matrix(, nrow = num_datapoints, ncol = num_datapoints)
tree_amat = matrix(, nrow = num_datapoints, ncol = num_datapoints)


for (i in 1:num_datapoints) {
  for (j in i:num_datapoints) {
    if (i == j) {
      edit_amat[i,j] = 0
    } else {
      dist = as.numeric(edit_dist_data[j + (i*(i-1)/2) - 1,3])
      edit_amat[i,j] = dist
      edit_amat[j,i] = dist
    }
  }
}

for (i in 1:num_datapoints) {
  for (j in i:num_datapoints) {
    if (i == j) {
      tree_amat[i,j] = 0
    } else {
      dist = as.numeric(tree_dist_data[j + (i*(i-1)/2) - 1,3])
      tree_amat[i,j] = dist
      tree_amat[j,i] = dist
    }
  }
}

row.names(edit_amat) = rev(node_data[,1])
colnames(edit_amat) = rev(node_data[,1])
row.names(tree_amat) = rev(node_data[,1])
colnames(tree_amat) = rev(node_data[,1])
row.names(node_data) = rev(node_data[,1])

cycombomapper(node_data, tree_amat, edit_amat, 15) # bin by tree distance, cluster by edit distance
