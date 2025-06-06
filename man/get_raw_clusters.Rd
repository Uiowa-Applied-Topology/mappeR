% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_factory.R
\name{get_raw_clusters}
\alias{get_raw_clusters}
\title{Ship data off to the clustering goblins}
\usage{
get_raw_clusters(dist_mats, clusterer)
}
\arguments{
\item{dist_mats}{A list of distance matrices of each bin that is to be clustered.}

\item{clusterer}{A function which accepts a list of distance matrices as input, and returns the results of clustering done on each distance matrix in a list.}
}
\value{
The output of \code{clusterer(dist_mats)}, which needs to be a list containing named vectors (one per bin), whose names are data point names and whose values are cluster labels (within each bin)
}
\description{
This function tells the computer to look away for a second, so the goblins come and cluster your data while it's not watching.
}
