## mappeR 1.0.0

I noticed some of the examples were breaking sporadically on my machine and did not seem to be caught during checks. This was due to the example being small enough to generate an edge case which was not dealt with properly. That has now been done. In a nutshell:

* mapper examples no longer break sporadically
* igraph example now always works (and function itself works on previous broken example)
* tests check for similar breaking cases (plus more)

And maybe less urgent:

* README now includes CRAN installation instructions

I don't anticipate releasing a new version to CRAN for a little while unless something else breaks; I snuck in some final code with this fix, hence the 1.0.0 versioning.

## R CMD check results

0 errors \| 0 warnings \| 0 notes
