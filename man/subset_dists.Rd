% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterer_farm.R
\name{subset_dists}
\alias{subset_dists}
\title{Subset a distance matrix}
\usage{
subset_dists(bin, dists)
}
\arguments{
\item{bin}{A list of names of data points.}

\item{dists}{A distance matrix for data points in the bin, possibly including extra points.}
}
\value{
A distance matrix for only the data points in the input bin.
}
\description{
Subset a distance matrix
}
