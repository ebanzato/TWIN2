#' Finding neighbors
#'
#' Finding neighbors of nodes in a graph.
#'
#' @param graph A symmetric adjacency matrix representing the undirected graph.
#'
#' @return a list. Each list element represents a graph node and contains its neighboring nodes.
#'
#' @export


find_neighbors = function(graph){

  if(!isSymmetric(graph)){
    stop('The adjacency matrix is not symmetric.')
  }
  g_nodes = colnames(graph)

  out = lapply(1:length(g_nodes), function(g){

    gcol = graph[, which(colnames(graph) == g_nodes[g])]
    gnei = g_nodes[which(gcol == 1)]
    gnei

  })

  names(out) = g_nodes

  return(out)
}
