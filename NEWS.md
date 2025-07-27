# mappeR (development version)

* Add more input sanitization
* Add more cluster metrics to output, namely medoid name, maximum distance to medoid, and sum of squared distances to medoid
* Examples now use actual noisy circle data; better examples coming soon

# mappeR 2.2.1

* Relax some tests as to avert floating point shenanigans

# mappeR 2.2.0

* Force user to have the names of the filtered data match the names of the original data set. You may need to call `names(filtered_data) = row.names(data)` or similar before making a call to generate a mapper object. 
The Ball Mapper flavored versions should still work fine.
* Add the option to specify a global cut height when using the global hierarchical clusterer.
* Update existing documentation and remove .Rd pages written about unexported functions (documentation still exists in source).

# mappeR 2.1.0

* Add better hierarchical clusterers; there is one which cuts each dendrogram uniformly, and one which cuts per patch. The cut height is at the beginning of the longest unbroken branch in the dendrogram (note this is not always the "correct" number of clusters!).
* Remove `igraph` dependency; see `shinymappeR` for an example of how to use `igraph` with `mappeR`.
* Change output vertex dataframe column names from `tightness` to `mean_dist_to_medoid` and `bin` to `patch`.
* Change edge `weight` to use Jaccard index.

# mappeR 2.0.2

* Really, for real this time, make the default clusterer cut correctly.

# mappeR 2.0.1

* Default hierarchical clusterer now correctly cuts dendrograms at the midpoint of their tallest branches.

# mappeR 2.0.0

* Clustering is now handled by a `clusterer`, which is a function that can handle a list of distance matrices (one for each bin/level set) and output clustering results for each one. The hierarchical clustering included previously is now available as a clusterer called `hierarchical_clusterer` because I am very creative.
* User-defined clusterers are possible as long as they can handle inputs/outputs correctly. The basic idea is that the output of a `clusterer` should look like a list of calls to `cutree` from the `hclust` package. Look to the clusterer farm for more examples in the future.

# mappeR 1.3.0

* Adjust `compute_tightness` to no longer normalize by the maximum distance from the medoid (easier to see behavior in a single mapper graph, may add options in future)
* Add option to consider each level set locally when clustering; default is still to do it globally

# mappeR 1.2.0

* Add other hierarchical clustering methods (those available from `fastcluster`)
* Fix issues with 100 percent overlap situations caused by `mapply` simplifications

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
