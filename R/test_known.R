
#' Local regressions and tests
#'
#' Function to fit node-conditional models for comparison and calculate p-values for the two stages.
#'
#' @param graph A symmetric adjacency matrix representing the undirected graph. Node names must match the variable names in the data.
#'
#' @param data A data frame with the observed variables. Variable names in the data must match the node names.
#'
#' @param group A vector of group labels for comparison, one entry per observation.
#'
#' @param glm.family Description of the distribution and link function to be used in the model.
#'
#' @param progressbar Logical, whether to display the progress bar. Default is TRUE.
#'
#' @return A list containing three elements: a matrix of p-values, a matrix of coefficients for the baseline group, and a matrix of coefficients for the differences between the two conditions.
#'
#' @export
#'


test_known = function(graph, data, group, glm.family, Z, progressbar){

  test.aov = 'Chisq'

  neigh = find_neighbors(graph)
  n_nodes = colnames(graph)

  if(progressbar){
    pb = utils::txtProgressBar(min = 0, max = length(n_nodes), style = 3)
  }


  results = lapply(1:length(n_nodes), function(j) {

    n_neigh = neigh[[j]]

    # out
    res.p = rep(NA, length(n_nodes)+2) # pvalues
    names(res.p) = c('LRT','group',n_nodes)
    res.c = res.p[-1] # coef
    res.b = res.c     # coef baseline

    if(length(n_neigh)==0){

      if(is.null(Z)){
        # define formula
        formula_h0 = paste0(n_nodes[j], ' ~ 1')
        formula_h1 = paste0(n_nodes[j], ' ~ group')
        # models
        mod_h0 = stats::glm(formula_h0, data.frame(data,'group'=group), family=glm.family)  # model H0
        mod_h1 = stats::glm(formula_h1, data.frame(data,'group'=group), family=glm.family)  # model H1
        z = 0
      }else{
        # define formula
        formula_h0 = paste0(n_nodes[j], ' ~ Z')
        formula_h1 = paste0(n_nodes[j], ' ~ Z + group')
        # models
        mod_h0 = stats::glm(formula_h0, data.frame(data,'group'=group,'Z'=Z), family=glm.family)  # model H0
        mod_h1 = stats::glm(formula_h1, data.frame(data,'group'=group,'Z'=Z), family=glm.family)  # model H1
        z = length(levels(as.factor(Z)))-1
      }

    }else{

      if(is.null(Z)){
        # define formula
        formula_h0 = paste0(n_nodes[j], ' ~ ', paste(neigh[[j]], collapse='+', sep=''))
        formula_h1 = paste0(n_nodes[j], ' ~ (', paste(neigh[[j]], collapse='+', sep=''),')*group')
        # models
        mod_h0 = stats::glm(formula_h0, data.frame(data,'group'=group), family=glm.family)  # model H0
        mod_h1 = stats::glm(formula_h1, data.frame(data,'group'=group), family=glm.family)  # model H1
        z = 0
      }else{
        # define formula
        formula_h0 = paste0(n_nodes[j], ' ~ Z +', paste(neigh[[j]], collapse='+', sep=''))
        formula_h1 = paste0(n_nodes[j], ' ~ Z + (', paste(neigh[[j]], collapse='+', sep=''),')*group')
        # models
        mod_h0 = stats::glm(formula_h0, data.frame(data,'group'=group,'Z'=Z), family=glm.family)  # model H0
        mod_h1 = stats::glm(formula_h1, data.frame(data,'group'=group,'Z'=Z), family=glm.family)  # model H1
        z = length(levels(as.factor(Z)))-1
      }
    }

    # Test
    test = stats::anova(mod_h0, mod_h1, test=test.aov)
    if(test.aov=='Chisq'){
      res.p[1] = test$`Pr(>Chi)`[2]
    }

    pvalbeta = summary(mod_h1)$coefficients[-c(1:(1+z+length(n_neigh))),4]
    # out
    res.p[2] = pvalbeta[1]
    res.p[-c(1:2)][which(n_nodes %in% neigh[[j]])] = pvalbeta[-1]

    # Coef delta
    coefdelta = summary(mod_h1)$coefficients[-c(1:(1+z+length(n_neigh))),1]
    # out
    res.c[1] = coefdelta[1]
    res.c[-1][which(n_nodes %in% neigh[[j]])] = coefdelta[-1]

    # Coef baseline (ref)
    coefbase = summary(mod_h1)$coefficients[c(1,(2+z):(1+z+length(n_neigh))),1]
    # out
    res.b[1] = coefbase[1]
    res.b[-1][which(n_nodes %in% neigh[[j]])] = coefbase[-1]



    if(progressbar){
      utils::setTxtProgressBar(pb, j)
    }

    # OUT
    res = list('p.mat'=res.p,'c.mat'=res.c, 'b.mat'=res.b)
    return(res)
  }
  )

  out.p = t(sapply(results, function(x) c(x$p.mat)))
  rownames(out.p) = n_nodes

  out.c = t(sapply(results, function(x) c(x$c.mat)))
  rownames(out.c) = n_nodes

  out.b = t(sapply(results, function(x) c(x$b.mat)))
  rownames(out.b) = n_nodes

  # out
  out = list('p.mat' = out.p, 'c.mat' = out.c, 'b.mat'=out.b)
  return(out)

}
