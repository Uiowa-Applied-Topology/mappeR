# mappeR 1.2.0

* Add other hierarchical clustering methods (those available from `fastcluster`)
* Fix issues with 100 percent overlap situtations caused by `mapply` simplifications

# mappeR 1.1.0

* Added safety checks for mapper input functions (no `NA` inputs, etc)
* Removed Cytoscape interfacing for the release version; will remain in a separate branch

# mappeR 1.0.0

* Balls/bins with equal amounts of data no longer cause wonkiness
* Zero overlap situations should work correctly with `igraph`
* Added eccentricity lens/filter fuction

# mappeR 0.1.5

* 1D tests will not fail due to a one interval cover

# mappeR 0.1.4

* Reduced edge widths for Cytoscape visualization.
* Rescaled tightness calculation down a bit.

# mappeR 0.1.3

* Edge overlap data is now present in output mapper objects.
* Cytoscape edges are styled by width instead of transparency.

# mappeR 0.1.2

* The algorithm now still runs even if there is only one level set to be considered.

# mappeR 0.1.1

# mappeR 0.1.0

* Initial CRAN submission.
