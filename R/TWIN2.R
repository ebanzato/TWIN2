
#' TWIN2 function
#'
#' Function to identify differences in networks between two conditions.
#'
#' 1. Screening stage
#'
#' 2. Confirmatory stage
#'
#' @param graph A symmetric adjacency matrix representing the undirected graph. Node names must match the variable names in the data.
#'
#' @param data A data frame with the observed variables. Variable names in the data must match the node names.
#'
#' @param group A vector of group labels for comparison, one entry per observation.
#'
#' @param glm.family Description of the distribution and link function to be used in the model.
#'
#' @param alpha The alpha level for conducting the test. Default is 0.05.
#'
#' @param method.FDR Method for adjusting p-values in the first stage to control the FDR. Default is 'BH'.
#'
#' @param method.FWER Method for adjusting p-values in the first stage to control the FWER. Default is 'holm'.
#'
#' @param Z .................
#'
#' @param progressbar Logical, whether to display the progress bar. Default is TRUE.
#'
#' @return a list with four elements.
#'
#' 1. diff: A matrix of dimension p x (p + 2), where p is the number of nodes, indicating differences in the network.
#' The columns are as follows:
#' 1st column: First screening stage: 0 for equality, 1 for nodes with a difference.
#' 2nd column: Difference in the mean parameter of the node in the row: 0 for no difference, 1 for a difference, NA if the first column is 0.
#' 3-(p+2)th columns: Differences in the association parameters between the node in the row and nodes in the columns: 0 for no difference, 1 for a difference, NA if not tested.
#'
#' 2. delta_pvalues: corrected p-values of the two stages.
#'
#' 3. delta_coef: coefficients of the difference parameters.
#'
#' 4. base_coef: coefficients of the reference condition.
#'
#' @export
#'


TWIN2 = function(graph, data, group, glm.family, alpha=0.05, method.FDR='BH', method.FWER='holm', Z=NULL, progressbar=TRUE){

  # Check if adjm names match
  if(!identical(colnames(graph), rownames(graph))){
    stop('No match between columns names and rows names of \'graph\'.')
  }

  # Check if 'graph' is symmetric
  if(!isSymmetric(graph)){
    warning('The adjacency matrix is not symmetric.')
  }

  # Check if data contains all the variables in graph
  if(!sum(colnames(graph) %in% colnames(data)) == length(colnames(graph))){
    stop('Not all variables in \'graph\' are in \'data\' OR no match between nodes and variables names')
  }

  # Check for variables that are all zeros
  if(sum(colSums(data) == 0) > 0){
    stop('One or more variables contain only zeros.')
  }

  # Check whether there are variables that contain only one unique value.
  if(sum(apply(data, 2, function(x) length(unique(x))) == 1) > 0){
    stop('One or more variables contain only one unique value.')
  }

  # Group as factors
  if(!is.null(group)){
    group = as.factor(group)
  }

  # Test
  if(glm.family == 'binomial' | glm.family=='poisson'){
    test.res = test_known(graph, data, group, glm.family, Z, progressbar)
  }


  # Two stages algorithm
  res = two_stages_pval(graph, test.pval = test.res$p.mat, method.FDR, method.FWER, alpha)

  # OUT
  out = list('g.diff' = res, 'delta_pvalues' = test.res$p.mat,
             'delta_coef' = test.res$c.mat, 'base_coef'= test.res$b.mat)
  return(out)
}

