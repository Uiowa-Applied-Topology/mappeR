mappeR
================
2025-02-16

<!-- badges: start -->

[![R-CMD-check](https://github.com/Uiowa-Applied-Topology/mappeR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Uiowa-Applied-Topology/mappeR/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

This is an implementation of the mapper algorithm by Singh, Mémoli, and
Carlsson (2007).

## Setup

To install and use the most recent CRAN upload of this package, run the
following:

`install.packages("mappeR")`

`library(mappeR)`

To install the latest development version of this package from Github,
run the following commands:

`install.packages("devtools")`

`library(devtools)`

`devtools::install_github("https://github.com/Uiowa-Applied-Topology/mappeR/tree/dev", upgrade=FALSE)`

`library(mappeR)`

If you’re installing from Github, you might need to do some more stuff:

- **Windows:** install Rtools
  (<https://cran.r-project.org/bin/windows/Rtools/>)
- **OS X:** install Xcode (from the Mac App Store)
- **Linux:** run `apt-get install r-base-dev` (or similar).

## Mathy Overview

Mapper is a way to view a point cloud $P$ through a “lens” of our
choice.

Consider a function

``` math
f: P \to \mathbb{R}
```

Cover $\mathbb{R}$ in a set of intervals $`\{I_i\}_{i=1}^n`$, so that
every point in $P$ is contained in some level set $L_i = f^{-1}(I_i)$.
We may then construct a graph

``` math
G = (V,E)
```

based on this cover to reflect the original data, where

``` math
V = \{L_i \mid L_i \neq \varnothing\}
```

and

``` math
E = \{\{L_i, L_j\}\mid L_i\cap L_j \neq \varnothing,\ i\neq j\}
```

This is the basic idea of the mapper algorithm, with the addition that
each level set is first refined into clusters based on the intrinsic
pairwise distances of the data according to some clustering algorithm.
That is, we partition each level set $L_i$ into $k_i$ disjoint clusters

``` math
L_i = \bigsqcup_{j=1}^{k_i} C_j
```

and build a new graph $G' = (V', E')$ that is homomorphic to $G$ defined
by

``` math
V' = \bigsqcup_{i=1}^{n}\{C_j\}_{j=1}^{k_i}
```

and

``` math
E' = \{\{C_i, C_j\}\mid C_i\cap C_j \neq \varnothing\}
```

So in general, the ingredients to construct a mapper graph are

- A data set, along with their pairwise distances
- A *lens* function with the data as its domain (above this was
  real-valued, but it does not have to be)
- A cover of the codomain of the lens function
- A clustering algorithm

## Example 1: 1D Mapper

``` r
num_points = 5000

P.data = data.frame(
  x = sapply(1:num_points, function(x)
    sin(x) * 10) + rnorm(num_points, 0, 0.1),
  y = sapply(1:num_points, function(x)
    cos(x) ^ 2 * sin(x) * 10) + rnorm(num_points, 0, 0.1),
  z = sapply(1:num_points, function(x)
    10 * sin(x) ^ 2 * cos(x)) + rnorm(num_points, 0, 0.1)
)

P.dist = dist(P.data)
```

Here is a point cloud $P$ formed by adding a bit of uniform noise to
5000 points regularly sampled from the parametric curve

``` math
\gamma(t) = \begin{cases}x = 10\ \sin(t)\\ y=10\ \sin(t)\ \cos^2(t)\\ z=10\ \sin^2(t)\ \cos(t) \end{cases}
```

<img src="README_files/figure-gfm/fig8-1.png"/>

This seems to form a kind of figure-8 curve just based on this
projection. But as we can see from the 2D projections, the “shape” of
the data set we’re seeing really does depend on how we’re looking at it:

<img src="README_files/figure-gfm/plotting_the_curve-1.png" style="display: block; margin: auto;" />

We will build graphs using the outline of the mapper algorithm
described, with real-valued lens functions.

Parameters:

- Data: figure-8
- Lens: Projection to each factor, or eccentricity (a measure of
  centrality per data point)
- Cover: A cover of $\mathbb{R}$ (just up to the extremes of the
  function values) using 10 equally spaced intervals with 25% overlap
  between each consecutive interval
- Clustering method: Single-linkage hierarchical clustering

``` r
# lens functions
projx = P.data$x
projy = P.data$y
projz = P.data$z
# eccentricity = apply(as.matrix(P.dist), 1, sum) / num_points
eccentricity = eccentricity_filter(P.dist)

# cover parameters to generate a width-balanced cover
num_bins = 10
percent_overlap = 25

# generate the cover
xcover = create_width_balanced_cover(min(projx), max(projx), num_bins, percent_overlap)
ycover = create_width_balanced_cover(min(projy), max(projy), num_bins, percent_overlap)
zcover = create_width_balanced_cover(min(projz), max(projz), num_bins, percent_overlap)
eccentriccover = create_width_balanced_cover(min(eccentricity),
                                             max(eccentricity),
                                             num_bins,
                                             percent_overlap)

# bin tester machine machine
check_in_interval <- function(endpoints) {
  return(function(x) (endpoints[1] - x <= 0) & (endpoints[2] - x >= 0))
}

# each of the "cover" elements will really be a function that checks if a data point lives in it
xcovercheck = apply(xcover, 1, check_in_interval)

# build the mapper objects
xmapper = create_mapper_object(
  data = P.data,
  dists = P.dist,
  filtered_data = projx,
  cover_element_tests = xcovercheck,
  clusterer = hierarchical_clusterer("single") # built-in mappeR method
)

# mappeR also has a built-in 1D mapper function that will do the above for you; a single linkage clustering scheme is selected by default
ymapper = create_1D_mapper_object(P.data, P.dist, projy, ycover)
zmapper = create_1D_mapper_object(P.data, P.dist, projz, zcover)
eccentricmapper = create_1D_mapper_object(P.data, P.dist, eccentricity, eccentriccover)

# mappeR also has functions which will convert the mapper outputs into igraph format
ixmapper = mapper_object_to_igraph(xmapper)
iymapper = mapper_object_to_igraph(ymapper)
izmapper = mapper_object_to_igraph(zmapper)
ieccentricmapper = mapper_object_to_igraph(eccentricmapper)
```

The vertices in each output graph below are colored according to the
level set the cluster belongs to, and scaled by (the square root of) the
number of data points in the cluster.

<img src="README_files/figure-gfm/mapping_the_mapper-1.png" width="50%" /><img src="README_files/figure-gfm/mapping_the_mapper-2.png" width="50%" /><img src="README_files/figure-gfm/mapping_the_mapper-3.png" width="50%" /><img src="README_files/figure-gfm/mapping_the_mapper-4.png" width="50%" />

## Example 2: ball mapper

By toying with the general mapper parameters, we can obtain different
flavors of the algorithm. In the *ball mapper* flavor, we simply use the
inclusion into the ambient space of the data as our lens function, and
let the cover do the work. Specifically, we cover the ambient space with
$\varepsilon$-balls by creating a $\varepsilon$-net, which can be done
with a greedy algorithm.

Parameters:

- Data: figure-8
- Cover: set of $\varepsilon$-balls in $\mathbb{R^3}$
- Lens function: inclusion from $P\hookrightarrow\mathbb{R}^3$
- Clustering method: none (or, “any data set is one big cluster”-type
  clustering)

There’s a secret parameter here, which is $\varepsilon$. Below are
output graphs for varying values of $\varepsilon$; the sizing is as with
the 1D mapper, but no coloring is done as each vertex would have to
receive its own color in this flavor, which is redundant.

``` r
# creates a cover using a greedy algorithm
balls1 = create_balls(data = P.data, dists = P.dist, eps = .25)

# ball tester machine machine
is_in_ball <- function(ball) {
  return(function(x) x %in% ball)
}

# filtering is just giving back the data (row names because my balls are lists of data point names, so the filter should match)
ballmapper1 = create_mapper_object(P.data, P.dist, rownames(P.data), lapply(balls1, is_in_ball))

# mappeR has a built-in ball mapper function to do this for you
ballmapper2 = create_ball_mapper_object(P.data, P.dist, .5)
ballmapper3 = create_ball_mapper_object(P.data, P.dist, 1)
ballmapper4 = create_ball_mapper_object(P.data, P.dist, 2)
```

<img src="README_files/figure-gfm/ballmapper_time-1.png" width="50%" /><img src="README_files/figure-gfm/ballmapper_time-2.png" width="50%" /><img src="README_files/figure-gfm/ballmapper_time-3.png" width="50%" /><img src="README_files/figure-gfm/ballmapper_time-4.png" width="50%" />

## Built-ins

# Mapper Flavors

`mappeR` has built-in methods for:

**1D mapper**

`create_1D_mapper_object(data, dists, filtered_data, cover, clusterer)`

- Lens: $P \to \mathbb{R}$
- Cover: intervals
- Clustering: yes

**Ball mapper**

`create_ball_mapper_object(data, dists, eps)`

- Lens: $P \to P$ by identity
- Cover: $\varepsilon$-balls in ambient $P$-space
- Clustering: no

**Clusterball mapper**

`create_clusterball_mapper_object(data, ball_dists, clustering_dists, eps, clusterer)`

- Lens: $P \to P$ by identity
- Cover: $\varepsilon$-balls in ambient $P$-space
- Clustering: yes

# Clustering

`mappeR` has a built in heuristic to perform any mode of hierarchical
clustering from the `hclust` package; the above example used single
linkage clustering this way, and `mappeR` will default to this, but you
can replace the keyword `"single"` with any of the following:

- `"single"`: single linkage
- `"complete"`: complete linkage
- `"average"`: average linkage
- `"mcquitty"`: McQuitty linkage
- `"centroid"`: centroid linkage
- `"median"`: median linkage
- `"ward.D"`: Ward linkage (assumes squared distances)
- `"ward.D2"`: Ward linkage (assumes nonsquared distances)

You may also build your own function to use for clustering; this
function must be able to accept a list of distance matrices and return a
list of cluster information vectors for each matrix. A “cluster
information vector” means a named vector (a vector whose elements have
names) with datapoint IDs for names and cluster IDs for values. The
cluster IDs should be positive integers but do not need to be unique
across all of the information vectors; `mappeR` will handle this for
you.
