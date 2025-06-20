% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cartography_chamber.R
\name{assemble_mapper_object}
\alias{assemble_mapper_object}
\title{Construct mapper graph from data}
\usage{
assemble_mapper_object(binclust_data, dists, binning = TRUE)
}
\arguments{
\item{binclust_data}{A list of bins, each containing named vectors whose names are those of data points and whose values are cluster ids}

\item{dists}{A distance matrix for the data that has been binned and clustered.}

\item{binning}{Whether the output dataframe should sort vertices into "bins" or not. Should be true if using clustering, leave false otherwise}
}
\value{
A list of two dataframes, one with node data containing bin membership,
datapoints per cluster, and mean distance to the medoid, and one with edge data
containing sources, targets, and weights representing overlap strength.
}
\description{
Construct mapper graph from data
}
