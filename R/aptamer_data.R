library(stringdist)
source("R/ballmapper.R")

m = as.matrix(read.csv("aptamer-data-raw.csv"))

sequences = m[,3]

amat = matrix(, nrow = nrow(m), ncol = nrow(m))
rownames(m) = m[,1]
rownames(amat) = m[,1]


for (i in 1:nrow(m)) {
  for (j in i:nrow(m)) {
    if (i == j) {
      amat[i,j] = 0
    } else {
      dist = stringdist(sequences[i], sequences[j], method = "lv")
      amat[i,j] = dist
      amat[j,i] = dist
    }
  }
}

dists = as.dist(amat)
cyballmapper(m, dists, 8)
