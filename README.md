mappeR
================
2024-09-27

This is an implementation of the
[mapper](https://research.math.osu.edu/tgda/mapperPBG.pdf) algorithm by
Singh, Mémoli, and Carlsson, and also the
[ballmapper](https://arxiv.org/pdf/1901.07410.pdf) algorithm from
Dlotko.

## Setup

To install the latest version of this package from Github, run the
following commands:

`install.packages("devtools")`

`library(devtools)`

`devtools::install_github("Uiowa-Applied-Topology/mappeR", upgrade=FALSE)`

`library("mappeR")`

If you’re installing from Github, you might need to do some more stuff:

- **Windows:** install Rtools
  (<http://cran.r-project.org/bin/windows/Rtools/>)
- **OS X:** install Xcode (from the Mac App Store)
- **Linux:** run `apt-get install r-base-dev` (or similar).

## Mathy Overview

Mapper is a way to view a point cloud $P$ through a “lens” of our
choice.

Consider a function $$f: P \to \mathbb{R}$$. Cover $\mathbb{R}$ in a set
of intervals $\{I_i\}_{i=1}^n$, so that every point in $P$ is contained
in some level set $L_i = f^{-1}(I_i)$. We may then construct a graph
$$G = (V,E)$$ based on this cover to reflect the original data, where
$$V = \{L_i \mid L_i \neq \varnothing\}$$, and
$$E = \{\{L_i, L_j\}\mid L_i\cap L_j \neq \varnothing,\ i\neq j\}$$

This is the basic idea of the mapper algorithm, with the addition that
each level set is first refined into clusters based on the intrinsic
pairwise distances of the data according to some clustering algorithm.
That is, we partition each level set $L_i$ into $k_i$ disjoint clusters
$\{C_j\}_{j=1}^{k_i}$ and build a new graph $G' = (V', E')$ that is
homomorphic to $G$ defined by
$$V' = \bigsqcup_{i=1}^{n}\{C_j\}_{j=1}^{k_i}$$ and
$$E' = \{\{C_i, C_j\}\mid C_i\cap C_j \neq \varnothing\}$$

So in general, the ingredients to construct a mapper graph are

- A data set, along with their pairwise distances
- A *lens* function with the data as its domain (above this was
  real-valued, but it does not have to be)
- A cover of the codomain of the lens function
- A clustering algorithm

## Example: 1D Mapper

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

$$\gamma(t) = \begin{cases}x = 10\sin(x)\\ y=10\sin(x)\cos^2(x)\\ z=10\sin^2(x)\cos(x) \end{cases}$$
<img src="README_files/figure-gfm/fig8-1.png" />

This seems to form a kind of figure-8 curve just based on this
projection. But as we can see from the 2D projections, the “shape” of
the data set we’re seeing really does depend on how we’re looking at it:

<img src="README_files/figure-gfm/plotting the curve-1.png" style="display: block; margin: auto;" />
We will build graphs using the outline of the mapper algorithm
described, with real-valued lens functions. The parameters used to
generate the graphs below were:

- Data: figure-8
- Lens: Projection to each factor, or eccentricity (a measure of
  centrality per data point)
- Cover: A cover of $\mathbb{R}$ (really, just up to the extremes of the
  function values) using 10 equally spaced intervals with 25% overlap
  between each consecutive interval
- Clustering method: Single-linkage hierarchical clustering

``` r
# lens functions
projx = P.data$x
projy = P.data$y
projz = P.data$z
eccentricity = apply(as.matrix(P.dist), 1, sum) / num_points

# cover parameters to generate a width-balanced cover
num_bins = 10
percent_overlap = 25

xmapper = get_1D_mapper_data(P.data, projx, P.dist, num_bins, percent_overlap, "single")
ymapper = get_1D_mapper_data(P.data, projy, P.dist, num_bins, percent_overlap, "single")
zmapper = get_1D_mapper_data(P.data, projz, P.dist, num_bins, percent_overlap, "single")
eccentricmapper = get_1D_mapper_data(P.data, eccentricity, P.dist, num_bins, percent_overlap, "single")
```

The vertices in each output graph below are colored according to the
level set the cluster belongs to, and scaled by (the square root of) the
number of data points in the cluster.

<img src="README_files/figure-gfm/mapping the mapper-1.png" width="50%" /><img src="README_files/figure-gfm/mapping the mapper-2.png" width="50%" /><img src="README_files/figure-gfm/mapping the mapper-3.png" width="50%" /><img src="README_files/figure-gfm/mapping the mapper-4.png" width="50%" />

## Example: ballmapper

By toying with the general mapper parameters, we can obtain different
flavors of the algorithm. In the *ballmapper* flavor, we simply use the
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
ballmapper1 = get_ballmapper_data(P.data, P.dist, .25)
ballmapper2 = get_ballmapper_data(P.data, P.dist, .5)
ballmapper3 = get_ballmapper_data(P.data, P.dist, 1)
ballmapper4 = get_ballmapper_data(P.data, P.dist, 2)
```

<img src="README_files/figure-gfm/ballmapper time-1.png" width="50%" /><img src="README_files/figure-gfm/ballmapper time-2.png" width="50%" /><img src="README_files/figure-gfm/ballmapper time-3.png" width="50%" /><img src="README_files/figure-gfm/ballmapper time-4.png" width="50%" />
