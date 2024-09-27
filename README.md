mappeR
================
2024-09-27

This is an implementation of the
[mapper](https://research.math.osu.edu/tgda/mapperPBG.pdf) algorithm by
Singh, Mémoli, and Carlsson, and also the
[ballmapper](https://arxiv.org/pdf/1901.07410.pdf) algorithm from
Dlotko.

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

## Examples

    ## Loading required package: rgl

    ## Loading required package: webshot2

<img src="../../AppData/Local/Temp/RtmpmgRuLP/file228fc304d395.png" width="50%" /><img src="../../AppData/Local/Temp/RtmpmgRuLP/file228fc696c3cc1.png" width="50%" /><img src="../../AppData/Local/Temp/RtmpmgRuLP/file228fc147a6446.png" width="50%" /><img src="../../AppData/Local/Temp/RtmpmgRuLP/file228fc7954365a.png" width="50%" /><img src="../../AppData/Local/Temp/RtmpmgRuLP/file228fc65d314be.png" width="50%" /><img src="../../AppData/Local/Temp/RtmpmgRuLP/file228fc36613e86.png" width="50%" /><img src="../../AppData/Local/Temp/RtmpmgRuLP/file228fc16695ec2.png" width="50%" /><img src="README_files/figure-gfm/plotting the curve-8.png" width="50%" />

<img src="README_files/figure-gfm/mapping the mapper-1.png" width="50%" /><img src="README_files/figure-gfm/mapping the mapper-2.png" width="50%" /><img src="README_files/figure-gfm/mapping the mapper-3.png" width="50%" />
