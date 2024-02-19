Cymapper
================
2024-02-19

This is an implementation of the
[mapper](https://research.math.osu.edu/tgda/mapperPBG.pdf) algorithm by
Singh, Mémoli, and Carlsson, and also the
[ballmapper](https://arxiv.org/pdf/1901.07410.pdf) algorithm from
Dlotko.

You will need to install the
[RCy3](https://www.bioconductor.org/packages/release/bioc/html/RCy3.html)
package. To do this, start R (or RStudio), and enter the following
commands:

`install.packages("BiocManager")`

`BiocManager::install("RCy3")`

You will also need to install [Cytoscape](https://cytoscape.org/) and
have it running in the background.

## Examples

Here we will run conventional 1D mapper on this noisy circle dataset,
using 10 bins with 15 percent overlap, and single linkage clustering.

``` r
circle.data = data.frame( x= sapply(1:1000, function(x) cos(x)) + rnorm(100, 500, .03), y = sapply(1:1000, function(x) sin(x)) + rnorm(100, 0, 0.03))
circle.dist = dist(circle.data)

# Cytoscape needs to be running
cymapper(circle.data, circle.data$x, circle.dist, 10, 15, "single")
```

    ## Loading data...

    ## Applying default style...

    ## Applying preferred layout...

![](man/figures/circle-1.png "noisy circle data")
![](man/figures/1dmappercircle.png "wow a cycle")

The size of the nodes are proportional to how many datapoints are in
that cluster. The border color identifies which bin the cluster belongs
to, while the interior color signals how “tight” the cluster is. The
edge opacity represents the relative overlap between clusters; a darker
edge indicates a stronger overlap.

Here we will run ballmapper on the same dataset, with
$\varepsilon = 0.3$.

``` r
cyballmapper(circle.data, circle.dist, .3)
```

    ## Loading data...

    ## Applying default style...

    ## Applying preferred layout...

![](man/figures/ballmappercircle.png "the balls have made a cycle")

The same visual styling was applied here, though it may be hard to see
the size differences. No colors are needed as the balls are also the
bins for ballmapper, so any coloring would be redundant.
