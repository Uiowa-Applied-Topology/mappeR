% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/art_studio.R
\name{mapper_object_to_igraph}
\alias{mapper_object_to_igraph}
\title{make igraph}
\usage{
mapper_object_to_igraph(mapperobject)
}
\arguments{
\item{mapperobject}{mapper object generated by mappeR}
}
\value{
an igraph object
}
\description{
make igraph
}
\examples{
data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))

projy = data$y

cover = create_width_balanced_cover(min(projy), max(projy), 10, 25)

mapperobj = create_1D_mapper_object(data, dist(data), data$y, cover)

mapper_object_to_igraph(mapperobj)
}
