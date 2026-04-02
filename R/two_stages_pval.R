#' Two-stage procedure
#'
#' 1. Screening stage
#'
#' 2. Confirmatory stage
#'
#' @param graph A symmetric adjacency matrix representing the undirected graph.
#'
#' @param test.pval A list returned by the 'test_known' function.
#'
#' @param method.FDR Method for adjusting p-values in the first stage to control the FDR.
#'
#' @param method.FWER Method for adjusting p-values in the first stage to control the FWER.
#'
#' @param alpha The alpha level for conducting the test.
#'
#' @return A matrix of dimension p x (p + 2), where p is the number of nodes, indicating differences in the network.
#' The columns are as follows:
#' 1. First screening stage: 0 for equality, 1 for nodes with a difference.
#' 2. Difference in the mean parameter of the node in the row: 0 for no difference, 1 for a difference, NA if the first column is 0.
#' 3-(p+2). Differences in the association parameters between the node in the row and nodes in the columns: 0 for no difference, 1 for a difference, NA if not tested.
#'
#' @export
#'


two_stages_pval = function(graph, test.pval, method.FDR, method.FWER, alpha){

  p = ncol(graph)

  # Stage 1: screening (FDR)
  stage1 = ifelse(stats::p.adjust(test.pval[,1], method=method.FDR) <= alpha, 1, 0)

  # Stage2: confirmatory (FWER)
  if(sum(stage1)>0){
    nRej = sum(stage1)
    alpha2 = (nRej/p)*alpha

    stage2 = t(sapply(1:p, function(r) {
      if(stage1[r] == 0){
        out = c(0, rep(NA, p+1))
      }else{
        out = c(1,ifelse(stats::p.adjust(test.pval[r,-1], method=method.FWER) <= alpha2, 1, 0))
      }
      out
    }))
  }else{
    stage2 = cbind(stage1, matrix(NA, ncol=p+1, nrow=p))

  }
  colnames(stage2) = colnames(test.pval)
  rownames(stage2) = rownames(test.pval)

  # out
  return(stage2)
}
