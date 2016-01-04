#does what the name suggests
#samples random values and stores the best
#' @export
random_search = function(fun, maxFE = NULL) {
  if (!"smoof" %in% rownames(installed.packages())) install.packages("smoof")
  require(smoof)
  #get box constraints
  
  #hardcode constraints since the smoof functions are buggy
  ub = 5
  lb = -5
  dimensions = length(getLowerBoxConstraints(fun))
  result = "Start"
  
  #if no stopping criterion specified, stop after 10000 function evaluations
  if (is.null(maxFE)) maxFE = 10000
  bestFitness = Inf
  for (i in 1:maxFE) {
    solution = runif(dimensions, lb, ub)
    fitness = fun(solution)
    if (fitness < bestFitness) {
      bestFitness = fitness
    }
    Fopt = getGlobalOptimum(fun)$value
    if ((i %% 10) == 0) result = c(result, paste(i, i, (bestFitness - Fopt)))
  }
  result = c(result, "End")
  return(result)
}