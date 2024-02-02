library(stringdist)
source("R/ballmapper.R")

data = as.matrix(read.csv("edit_dist_edges.csv"))
node_data = as.matrix(read.csv("aptamer-data-raw.csv"))

amat = matrix(, nrow = nrow(m), ncol = nrow(m))


for (i in 1:nrow(m)) {
  for (j in i:nrow(m)) {
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
row.names(m) = rev(node_data[,1])

cyballmapper(m, amat, 8)
