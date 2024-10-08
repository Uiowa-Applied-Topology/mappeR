% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cartography_chamber.R
\name{get_edgelist_from_overlaps}
\alias{get_edgelist_from_overlaps}
\title{Obtain edge list from cluster intersections}
\usage{
get_edgelist_from_overlaps(overlaps, num_vertices)
}
\arguments{
\item{overlaps}{A named list of edges, whose elements contain the names of clusters in the overlap represented by that edge; output of \code{\link[=get_overlaps]{get_overlaps()}}.}

\item{num_vertices}{The number of vertices in the graph.}
}
\value{
A 2D array representing the edge list of a graph.
}
\description{
Obtain edge list from cluster intersections
}
