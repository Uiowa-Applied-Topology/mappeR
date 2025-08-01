% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cartography_chamber.R
\name{create_mapper_object}
\alias{create_mapper_object}
\title{Mapper}
\usage{
create_mapper_object(
  data,
  dists,
  filtered_data,
  cover_element_tests,
  clusterer = NULL
)
}
\arguments{
\item{data}{A data frame.}

\item{dists}{A distance matrix for the data frame. Can be a \code{dist} object or \code{matrix}.}

\item{filtered_data}{The result of a function applied to the data frame; there should be one filter value per observation in the original data frame.
These values need to be named, and the names of these values must match the names of the original data set.}

\item{cover_element_tests}{A list of membership test functions for a set of cover elements. In other words, each element of \code{cover_element_tests} is a function that returns \code{TRUE} or \code{FALSE} when given a filter value.}

\item{clusterer}{A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix;
that is, it should return a list of named vectors, whose names are the names of data points and whose values are cluster assignments (integers).
If this value is omitted, then trivial clustering will be done.}
}
\value{
A \code{list} of two data frames, one with node data and one with edge data. The node data includes:
\itemize{
\item \code{id}: vertex ID
\item \code{cluster_size}: number of data points in cluster
\item \code{medoid}: the name of the medoid of the vertex
\item \code{mean_dist_to_medoid}: mean distance to medoid of cluster
\item \code{max_dist_to_medoid}: max distance to medoid of cluster
\item \code{cluster_width}: maximum pairwise distance within cluster
\item \code{wcss}: sum of squares of distances to cluster medoid
\item \code{data}: names of data points in cluster
\item \code{patch}: level set ID
}

The edge data includes:
\itemize{
\item \code{source}: vertex ID of edge source
\item \code{target}: vertex ID of edge target
\item \code{weight}: Jaccard index of edge; intersection divided by union
\item \code{overlap_data}: names of data points in overlap
\item \code{overlap_size}: number of data points overlap
}
}
\description{
Run the Mapper algorithm on a data frame.
}
\examples{
# Create noisy data around a circle
data = data.frame(x = sapply(1:1000, function(x) cos(x)) + runif(1000, 0, .25),
 y = sapply(1:1000, function(x) sin(x)) + runif(1000, 0, .25))

# Apply lens function to data
projx = data$x
names(projx) = row.names(data)

# Build a width-balanced cover with 10 intervals and 25 percent overlap
num_bins = 10
percent_overlap = 25
xcover = create_width_balanced_cover(min(projx), max(projx), num_bins, percent_overlap)

# Write a function which can check if a data point lives in an interval of our lens function
check_in_interval <- function(endpoints) {
 return(function(x) (endpoints[1] - x <= 0) & (endpoints[2] - x >= 0))
}

# Each of the "cover" elements will really be a function that checks if a data point lives in it
xcovercheck = apply(xcover, 1, check_in_interval)

# Build the mapper object
xmapper = create_mapper_object(
  data = data,
  dists = dist(data),
  filtered_data = projx,
  cover_element_tests = xcovercheck
)
}
