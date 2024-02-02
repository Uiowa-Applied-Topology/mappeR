library(stringdist)
source("R/ballmapper_viz.R")

m = as.matrix(read.csv("edit_dist_edges.csv"))
node_data = as.matrix(read.csv("aptamer-data-raw.csv"))

num_datapoints = nrow(node_data)

amat = matrix(, nrow = num_datapoints, ncol = num_datapoints)


for (i in 1:num_datapoints) {
  for (j in i:num_datapoints) {
    if (i == j) {
      amat[i,j] = 0
    } else {
      dist = as.numeric(data[j + (i*(i-1)/2) - 1,3])
      amat[i,j] = dist
      amat[j,i] = dist
    }
  }
}

row.names(amat) = rev(node_data[,1])
colnames(amat) = rev(node_data[,1])
row.names(node_data) = rev(node_data[,1])

cyballmapper(node_data, amat, 4)
