# runner function for combo mapper; outputs bins, clusters, and the mapper graph.
get_combomapper_data <- function(data, dist1, dist2, eps) {
  balls = create_balls(data, dist1, eps)
  formatted_balled_data = convert_balls(balls)

  clusters = get_clusters(formatted_balled_data, dist2, "single")

  combomappergraph = construct_1Dmappergraph(clusters, dist2)

  return(combomappergraph)
}

cycombomapper <- function(data, dist1, dist2, eps) {
  visualize_mapper_data(get_combomapper_data(data, dist1, dist2, eps), FALSE)

  return(invisible(NULL))
}
