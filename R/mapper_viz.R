# TODO: add some DOCS

visualize_mapper_data <- function(mapper_data, is_ballmapper = TRUE) {
  nodes = mapper_data[[1]]
  edges = mapper_data[[2]]

  createNetworkFromDataFrames(nodes, edges)

  style.name = paste("mapperstyle", runif(1))
  defaults <- list(
    NODE_SHAPE = "ellipse",
    NODE_BORDER_WIDTH = 10,
    NODE_BORDER_PAINT = "#000",
    EDGE_WIDTH = 10
  )

  nodeSizes <- mapVisualProperty('node size',
                                 'id',
                                 'd',
                                 1:nrow(nodes),
                                 100 * nodes$size / max(nodes$size))
  edgeWidth <- mapVisualProperty('edge transparency', 'weight', 'c', c(0, .5, 1), c(0, 127, 255))
  nodeFillColors <- mapVisualProperty(
    'node fill color',
    'tightness',
    'c',
    c(0, mean(nodes$tightness), 1),
    c("#ffffff", "#efefef", "#000000")
  )

  if (is_ballmapper) {
    # ballmapper needs no more styling
    createVisualStyle(style.name,
                      defaults,
                      list(nodeSizes, edgeWidth, nodeFillColors))
  } else {
    # conventional mapper needs bin coloring
    num_bins = length(unique(nodes$bin))
    colfunc <- colorRampPalette(c("blue", "purple", "red"))
    bin_colors = colfunc(num_bins)
    nodeBorderColors <- mapVisualProperty('node border color', 'bin', 'd', 1:num_bins, bin_colors)
    createVisualStyle(
      style.name,
      defaults,
      list(nodeSizes, edgeWidth, nodeBorderColors, nodeFillColors)
    )
  }

  setVisualStyle(style.name)
}
