% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baskin_robbins.R
\name{create_ball_mapper_object}
\alias{create_ball_mapper_object}
\title{Ball Mapper}
\usage{
create_ball_mapper_object(data, dists, eps)
}
\arguments{
\item{data}{A data frame.}

\item{dists}{A distance matrix for the data frame. Can be a \code{dist} object or a \code{matrix}.}

\item{eps}{A positive real number for the desired ball radius.}
}
\value{
A \code{list} of two data frames, \code{nodes} and \code{edges}, which contain information about the Mapper graph constructed from the given parameters.

The node data frame consists of:
\itemize{
\item \code{id}: vertex ID
\item \code{cluster_size}: number of data points in cluster
\item \code{medoid}: the name of the medoid of the vertex
\item \code{mean_dist_to_medoid}: mean distance to medoid of cluster
\item \code{max_dist_to_medoid}: max distance to medoid of cluster
\item \code{cluster_width}: maximum pairwise distance within cluster
\item \code{wcss}: sum of squares of distances to cluster medoid
\item \code{data}: names of data points in cluster
}

The \code{edge} data frame contains consists of:
\itemize{
\item \code{source}: vertex ID of edge source
\item \code{target}: vertex ID of edge target
\item \code{weight}: Jaccard index of edge; this is the size of the intersection between the vertices divided by the union
\item \code{overlap_data}: names of data points in overlap
\item \code{overlap_size}: number of data points overlap
}
}
\description{
Run Mapper using the identity function as a lens and an \eqn{\varepsilon}-net cover, greedily generated using a distance matrix.
}
\examples{
# Create noisy cirle data set
data = data.frame(x = sapply(1:1000, function(x) cos(x)) + runif(1000, 0, .25),
y = sapply(1:1000, function(x) sin(x)) + runif(1000, 0, .25))

# Set ball radius
eps = .25

# Create Mapper object
create_ball_mapper_object(data, dist(data), eps)
}
