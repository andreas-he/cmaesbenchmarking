### DELETE THIS SECTION WHEN PACKAGE IS CREATED ###
#if (!"devtools" %in% rownames(installed.packages())) install.packages("devtools")
#if (!"snow" %in% rownames(installed.packages())) install.packages("snow")
#if (!"parallel" %in% rownames(installed.packages())) install.packages("parallel")
#if (!"smoof" %in% rownames(installed.packages())) install.packages("smoof")
#if (!"BBmisc" %in% rownames(installed.packages())) install.packages("BBmisc")
#require(devtools)
#install_github(repo = "MarcusCramer91/cmaesr")
#require(cmaesr)
#require(BBmisc)
#require(snow)
#require(parallel)
#require(smoof)
#source("./cmaes/cmaes.R")
#source("rbga.R")
### END ###

#' @name bbo_benchmarking
#' @title Customized BBo-Benchmarking
#' @aliases bbob_custom_parallel
#' @aliases bbob_custom
#' @description
#' provides an interface for benchmarking optimization algorithms over all 24 noiseless bbob functions.
#' Use either \code{bbob_custom} or \code{bbob_custom_parallel} for parallel execution.
#'
#' @details
#' bbob_custom provides a benchmarking function as a highlevel interface for running experiments
#' for an optimizer over all 24 noiseless functions of the bbob testset. bbob_custom_parallel does the same in an 
#' parallel manner in order to save execution time. The input parameters need to be specified as defined
#' above. Once the execution of an experiment has finished, bbob_custom created a folder according to the datapath specified and  
#' writes the results in a number of log files to this folder. The folder's name includes a time stamp as well as the specified
#' data directory (e.g. 2015-12-31_CMAES_default_stopping_conditions). The caption of the log files in that folder is a sequence of the defined
#' \code{algorithm.id, function_ids} and \code{dimensions} (e.g. CMAES_output_1_2.txt).
#' @references
#' TODO
#' @keywords optimize
#' @param optimizer 
#'   The first argument passed to customized_bbob is an optimization wrapper, i.e. the particular optimizer under test. 
#'   The optimizer has to adhere to the following function signature:
#'   \code{function(dimension, instance, function_id, maxit, maxFE, stopFitness, path, debug.logging, max_restarts, 
#'   restart_multiplier, restart_triggers, OCD, varLimit, nPreGen, maxGen, fitnessValue, dispersion, evolutionPath)}.
#'   However, the user does not need to specify every argument in an optimizer but it must be able to
#'   deal with functions of the type \code{smoof_function}, that gets \code{dimension, instance, function_id} as 
#'   input parameters.
#' @param algorithm.id
#'   The \code{algorithm.id} is a short descriptive name for the optimization algorithm under test.
#'   The \code{algorithm.id} will be part of the caption of the output log files. See detail section.
#' @param data_directory
#'   The \code{data_directory} specifies the location of the output log files. See detail section for more information.
#' @param dimensions
#'   The \code{dimensions} parameter determines the dimensions to be passed to the optimizer under test. Possible values are between
#'   2 and 40.
#'   Default is \code{c(2, 3, 5, 10, 20, 40)}.
#' @param instances
#'   Every \code{smoof_function} can be instantiated in different ways. Could be either a single intergerish value 
#'   greater of equals 1 or a vector with such values. Default is \code{c(1:15)}
#' @param function_ids
#'   The \code{function_ids} are the unique identifier of the 24 noiseless bbob function. 
#'   Default is \code{c(1:24)}.
#' @param maxit
#'   If \code{maxFE} is not passed, the number of iterations, \code{maxit}, serves as an upper execution limit for optimizing one function.
#'   Default is \code{NULL}.
#' @param stopFitness
#'   If \code{stopFitness} is specified, the optimization is terminated if the gap between the current function
#'   value and the global optimum of the function is below the value of \code{stopFitness}. 
#'   Default is \code{NULL}.
#' @param maxFE
#'   \code{maxFE} is another upper execution limit for a function optimization. The limit \code{maxFE} is favored over
#'   maxit/stopfitness, if both are not \code{NULL}.
#'   Default is \code{NULL}.
#' @param debug.logging
#'   \code{debug.logging} ...
#'   Default is \code{FALSE}.
#' @param max.restarts
#'   If \code{max.restarts} the optimizer restarts the optimization if one of the passed \code{restart_triggers} has fired.
#'   In order to use this feature, the optimizer under test must be able to restart the optimization.
#'   Default is \code{FALSE}.
#' @param restart.triggers
#'   If \code{restart.triggers} are passed, the optimizer restarts the optimization if one of those triggers has fired.
#'   In order to use this feature, the optimizer under test must be able to restart the optimization.
#'   Default is \code{FALSE}.
#' @param OCD
#'   \code{OCD} indicates if Online Convergence Detection (OCD) should be used as a stopping condition of the
#'   optimizer under test. OCD can only be activated if the optimizer has implemented the functionality necessary for OCD.
#'   See \code{\link{stopOnOCD}} for further necessary parameters.
#'   Default is \code{FALSE}.
#' @param varLimit
#'   \code{varLimit}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @param nPreGen
#'   \code{nPreGen}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @param maxGen
#'   \code{maxGen}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @param fitnessValue
#'   \code{fitnessValue}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @param dispersion
#'   \code{disperion}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @param evolutionPath
#'   \code{evolutionPath}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @return bbob_custom does not return anything but writes the results of the experiment to log files, to be
#' processed with \code{\link{readOutput}}
#' @examples
#' suppressWarnings(bbob_custom_parallel(optimizer = optimizerCMAES, algorithm_id = "CMAES_OCD", 
#' data_directory = "CMAES_OCD_no_restarts", 
#' dimensions = c(2, 5, 10, 20), instances = 1:15, function_ids = 1:24, maxit = NULL, 
#' stopFitness = 1e-08, maxFE = 100000, max_restarts = 0, 
#' OCD = TRUE, varLimit = 0.0001, nPreGen = 100, fitnessValue = TRUE, 
#' dispersion = FALSE,  evolutionPath = FALSE, restart_multiplier = 2, 
#' restart_triggers = "OCD"))
#' @rdname bbo-benchmarking
#' @importFrom BBmisc makeProgressBar
#' @export
#only non-noisy functions
bbob_custom = function(optimizer, algorithm_id, data_directory, dimensions = c(2, 3, 5, 10, 20, 40), 
                       instances = c(1:5, 41:50), function_ids = NULL, maxit = NULL, stopFitness = NULL, 
                       maxFE = NULL, debug.logging = FALSE, max_restarts = 0, 
                       restart_multiplier = 1, restart_triggers = character(0), OCD = FALSE, varLimit = NULL,
                       nPreGen = NULL, maxGen = NULL, fitnessValue = FALSE, dispersion = FALSE, evolutionPath = FALSE) {
  write(paste("Functions:", function_ids), file = "bbob_calls.txt", append = TRUE)
  write(paste("Dimensions:", dimensions), file = "bbob_calls.txt", append = TRUE)
  write(paste("Instances:", instances), file = "bbob_calls.txt", append = TRUE)
  write("================================", file = "bbob_calls.txt", append = TRUE)
  
  data_directory = paste(Sys.Date(), data_directory, sep = "_")
  dir.create(data_directory, showWarnings = FALSE)
  dimensions = sort(dimensions, decreasing = FALSE)
  if (is.null(function_ids)) {
    function_ids = 1:24
  }
  
  #some sanity checks
  if (is.null(c(maxit, maxFE)) && !is.null(stopFitness)) stop("To ensure termination, stopFitness must be combined with either maxit or maxFE")
  if (OCD == TRUE & (is.null(varLimit) | is.null(nPreGen) | (isFALSE(fitnessValue) & isFALSE(dispersion) & isFALSE(evolutionPath))))
    stop("If OCD is enabled, a value for varLimit and nPreGen must be passed and at least one performance indicator must be enabled.")
  
  nruns = length(function_ids)*length(dimensions)*length(instances)
  currentRun = 1
  pbar = makeProgressBar(min = 1, max = nruns)
  for (i in 1:length(function_ids)) {
    for (j in 1:length(dimensions)) {
      for (k in 1:length(instances)) {
        print(paste("Function:", function_ids[i], ",instance:", instances[k], ",dimensions:", dimensions[j]))
        result = optimizer(dimension = dimensions[j], instance = instances[k], function_id = function_ids[i], 
                           maxit = maxit, maxFE = maxFE, stopFitness = stopFitness, path = data_directory, 
                           debug.logging = debug.logging, max_restarts, restart_multiplier,
                           restart_triggers, OCD = OCD, varLimit = varLimit, nPreGen = nPreGen, maxGen = maxGen,
                           fitnessValue = fitnessValue, dispersion = dispersion, evolutionPath = evolutionPath)
        pbar$set(currentRun)
        currentRun = currentRun + 1
        outputFile = file.path(data_directory, paste(algorithm_id, "_output_", function_ids[i], "_", dimensions[j], 
                                                     ".txt", sep = ""))
        write(result, file = outputFile, append = TRUE)
      }
    }
  }
}

#' @name bbo_benchmarking_optimizers
#' @title Optimizer for BBo_Benchmarking
#' @aliases optimizerCMAES
#' @aliases optimizerCMAESWithoutDef
#' @aliases optimizerRS
#' @aliases optimizerGA
#' @param dimensions [\code{integer}]\cr
#'   The \code{dimensions} parameter determines the dimensions to be passed to the optimizer under test. Possible values are between
#'   2 and 40.
#' @param instances [\code{integer}]\cr
#'   Every \code{smoof_function} can be instantiated in different ways. Could be a single intergerish value 
#'   greater of equals 1.
#' @param function_ids [\code{integer(1)}]\cr
#'   The \code{function_ids} are the unique identifier of the 24 noiseless bbob function. 
#' @param maxit [\code{integer}]\cr
#'   If \code{maxFE} is not passed, the number of iterations, \code{maxit}, serves as an upper execution limit for optimizing one function.
#' @param stopFitness [\code{numeric}]\cr
#'   If \code{stopFitness} is specified, the optimization is terminated if the gap between the current function
#'   value and the global optimum of the function is below the value of \code{stopFitness}. 
#' @param maxFE [\code{interger}]\cr
#'   \code{maxFE} is another upper execution limit for an optimizer, i.e. if the specified number of function evaluations is reached,
#'   the optimizer terminates the optimization.
#' @param debug.logging [\code{logical}]\cr
#'   \code{debug.logging} ...
#'   Default is \code{FALSE}.
#' @param max.restarts [\code{integer}]\cr
#'   If \code{max.restarts} the optimizer restarts the optimization if one of the passed \code{restart_triggers} has fired.
#'   In order to use this feature, the optimizer under test must be able to restart the optimization.
#'   Default is \code{FALSE}.
#' @param restart.triggers [\code{string}]\cr
#'   If \code{restart.triggers} are passed, the optimizer restarts the optimization if one of those triggers has fired.
#'   In order to use this feature, the optimizer under test must be able to restart the optimization.
#'   Default is \code{FALSE}.
#' @param OCD [\code{logical}]\cr
#'   \code{OCD} indicates if Online Convergence Detection (OCD) should be used as a stopping condition of the
#'   optimizer under test. OCD can only be activated if the optimizer has implemented the functionality necessary for OCD.
#'   See \code{\link{stopOnOCD}} for further necessary parameters.
#'   Default is \code{FALSE}.
#' @param varLimit [\code{logical}]\cr
#'   \code{varLimit}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @param nPreGen [\code{logical}]\cr
#'   \code{nPreGen}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @param maxGen [\code{logical}]\cr
#'   \code{maxGen}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @param fitnessValue [\code{logical}]\cr
#'   \code{fitnessValue}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @param dispersion [\code{logical}]\cr
#'   \code{disperion}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @param evolutionPath [\code{logical}]\cr
#'   \code{evolutionPath}: OCD Parameter (See \code{\link{stopOnOCD}})
#' @description
#' Optimizers ready for usage in bbob_custom or bbob_custom parallel
#' @details
#' The optimizers are ready to be set up in a benchmarking experiment with \code{\link{bbob_custom}} or
#' \code{\link{bbob_custom_parallel}}. Every input parameter is required as it passed by
#' \code{\link{bbob_custom}} or \code{\link{bbob_custom_parallel}} to the particular optimizer.
#' @export
optimizerCMAES = function(dimension, instance, function_id, maxit, maxFE, stopFitness, path, 
                          debug.logging = FALSE, max_restarts, restart_multiplier, restart_triggers, OCD = FALSE,
                          varLimit = NULL, nPreGen = NULL, maxGen = NULL, fitnessValue= FALSE, dispersion = FALSE,
                          evolutionPath = FALSE) {
  if(!grepl(":/", path, fixed = TRUE)) path = file.path(getwd(), path)
  fun = makeBBOBFunction(dimension = dimension, fid = function_id, iid = instance)
  #create .txt creating monitor
  Fopt = getGlobalOptimum(fun)$value
  monitor = makeTXTMonitor(max.params = 4L, path, Fopt, function_id, dimension, instance)
  #use maxFE before maxit/stopfitness if both are not null
  condition1 = NULL
  condition2 = NULL
  #use maxit before stopfitness
  if (!is.null(maxFE)) condition1 = stopOnMaxEvals(maxFE)
  else if (!is.null(maxit)) condition1 = stopOnMaxIters(maxit)
  #stopFitness can only be used in combination with either maxFE or maxit (caught error)
  result = NULL
  if (OCD == TRUE) {
    OCDcond = stopOnOCD(varLimit = varLimit, nPreGen = nPreGen, maxGen = maxGen, 
                                       fitnessValue = fitnessValue, dispersion = dispersion, evolutionPath = evolutionPath)
    #also add the default stopping conditions to the restart triggers
    #these need to be enabled at the same time as they perform some sanity checks
    restart_triggers = c(restart_triggers, "tolX", "noEffectAxis", "noEffectCoord",
                                             "conditionCov", "indefCovMat")
  }
  if (!is.null(stopFitness)) {
    optValue = getGlobalOptimum(fun)$value
    condition2 = stopOnOptValue(optValue, stopFitness)
    if (OCD == FALSE) {
      result = cmaes_custom(fun, monitor = monitor, debug.logging = debug.logging,
                   control = list (stop.ons = c(list(condition1, condition2),
                                             getDefaultStoppingConditions()), 
                                             max.restarts = max_restarts,
                                             restart.triggers = restart_triggers,
                                             restart.multiplier = restart_multiplier))
    }
    else {
      result = cmaes_custom(fun, monitor = monitor, debug.logging = debug.logging,
                            control = list (stop.ons = c(list(condition1, condition2, OCDcond), getDefaultStoppingConditions()), 
                                            max.restarts = max_restarts,
                                            restart.triggers = restart_triggers,
                                            restart.multiplier = restart_multiplier))
    }
  }
  #if stop fitness is null
  else if (!is.null(condition1)) {
    if (OCD == FALSE) {
     result = cmaes_custom(fun, monitor = monitor, debug.logging = debug.logging, 
                            control = list (stop.ons = c(list(condition1), 
                                                         getDefaultStoppingConditions()),
                                            max.restarts = max_restarts,
                                            restart.triggers = restart_triggers,
                                            restart.multiplier = restart_multiplier))
    }
    else { 
      result = cmaes_custom(fun, monitor = monitor, debug.logging = debug.logging, 
                            control = list (stop.ons = c(list(condition1, OCDcond), getDefaultStoppingConditions()),
                                            max.restarts = max_restarts,
                                            restart.triggers = restart_triggers,
                                            restart.multiplier = restart_multiplier))    
    }
  }
  #use default if no stopping criterion is defined
  else result = cmaes_custom(fun, monitor = monitor, debug.logging = debug.logging)
  return(result)
}

#' @rdname bbo_benchmarking_optimizers
#reduced optimizer that does not apply default stopping conditions
#purely for showcasing purposes for the presentation, should never be used in practise
optimizerCMAESWithoutDef = function(dimension, instance, function_id, maxit, maxFE, stopFitness, path, 
                          debug.logging = FALSE, max_restarts, restart_multiplier, restart_triggers, OCD = FALSE,
                          varLimit = NULL, nPreGen = NULL, maxGen = NULL, fitnessValue= FALSE, dispersion = FALSE,
                          evolutionPath = FALSE) {
  if(!grepl(":/", path, fixed = TRUE)) path = file.path(getwd(), path)
  fun = makeBBOBFunction(dimension = dimension, fid = function_id, iid = instance)
  #create .txt creating monitor
  Fopt = getGlobalOptimum(fun)$value
  monitor = makeTXTMonitor(max.params = 4L, path, Fopt, function_id, dimension, instance)
  #use maxFE before maxit/stopfitness if both are not null
  condition1 = NULL
  condition2 = NULL
  #use maxit before stopfitness
  if (!is.null(maxFE)) condition1 = stopOnMaxEvals(maxFE)
  else if (!is.null(maxit)) condition1 = stopOnMaxIters(maxit)
  #stopFitness can only be used in combination with either maxFE or maxit (caught error)
  result = NULL
  if (OCD == TRUE) {
    OCDcond = stopOnOCD(varLimit = varLimit, nPreGen = nPreGen, maxGen = maxGen, 
                        fitnessValue = fitnessValue, dispersion = dispersion, evolutionPath = evolutionPath)
    #also add the default stopping conditions to the restart triggers
    #these need to be enabled at the same time as they perform some sanity checks
    restart_triggers = c(restart_triggers, "tolX", "noEffectAxis", "noEffectCoord",
                         "conditionCov", "indefCovMat")
  }
  if (!is.null(stopFitness)) {
    optValue = getGlobalOptimum(fun)$value
    condition2 = stopOnOptValue(optValue, stopFitness)
    if (OCD == FALSE) {
      result = cmaes_custom(fun, monitor = monitor, debug.logging = debug.logging,
                            control = list (stop.ons = c(list(condition1, condition2)), 
                                            max.restarts = max_restarts,
                                            restart.triggers = restart_triggers,
                                            restart.multiplier = restart_multiplier))
    }
    else {
      result = cmaes_custom(fun, monitor = monitor, debug.logging = debug.logging,
                            control = list (stop.ons = c(list(condition1, condition2, OCDcond)), 
                                            max.restarts = max_restarts,
                                            restart.triggers = restart_triggers,
                                            restart.multiplier = restart_multiplier))
    }
  }
  #if stop fitness is null
  else if (!is.null(condition1)) {
    if (OCD == FALSE) {
      result = cmaes_custom(fun, monitor = monitor, debug.logging = debug.logging, 
                            control = list (stop.ons = c(list(condition1)),
                                            max.restarts = max_restarts,
                                            restart.triggers = restart_triggers,
                                            restart.multiplier = restart_multiplier))
    }
    else { 
      result = cmaes_custom(fun, monitor = monitor, debug.logging = debug.logging, 
                            control = list (stop.ons = c(list(condition1, OCDcond)),
                                            max.restarts = max_restarts,
                                            restart.triggers = restart_triggers,
                                            restart.multiplier = restart_multiplier))    
    }
  }
  #use default if no stopping criterion is defined
  else result = cmaes_custom(fun, monitor = monitor, debug.logging = debug.logging)
  return(result)
}

#' @rdname bbo_benchmarking_optimizers
optimizerRS = function(dimension, instance, function_id, maxit, maxFE, stopFitness, path, OCD = FALSE,
                       debug.logging = FALSE, max_restarts = 0, 
                       restart_multiplier = 1, restart_triggers = character(0), ...) {
  source("random_search.R")
  fun = makeBBOBFunction(dimension = dimension, fid = function_id, iid = instance)
  result = random_search(fun, maxFE)
  return(result)
}

#' @rdname bbo_benchmarking_optimizers
optimizerGA = function(dimension, instance, function_id, maxit, maxFE, stopFitness, path, 
                       debug.logging = FALSE, max_restarts, restart_multiplier, restart_triggers, OCD = FALSE,
                       varLimit = NULL, nPreGen = NULL, maxGen = NULL, fitnessValue= FALSE, dispersion = FALSE,
                       evolutionPath = FALSE) {
  if(!grepl(":/", path, fixed = TRUE)) path = file.path(getwd(), path)
  fun = makeBBOBFunction(dimension = dimension, fid = function_id, iid = instance)
  #create .txt creating monitor
  npop = dimension * 10
  # there must be a value for maxit. In order to guarantee the upper bound of 100000 function evaluations, set maxit to maxFE
  if (!is.null(maxFE) && is.null(maxit)) maxit = maxFE
  Fopt = getGlobalOptimum(fun)$value
  monitor = makeTXTMonitor(max.params = 4L, path, Fopt, function_id, dimension, instance)
  #use maxFE before maxit/stopfitness if both are not null
  condition1 = NULL
  condition2 = NULL
  if(!is.null(maxFE)) condition1 = stopOnMaxEvals(maxFE)
  #stopFitness can only be used in combination with either maxFE or maxit (caught error)
  result = NULL
  if (OCD == TRUE) {
    OCDcond = stopOnOCD(varLimit = varLimit, nPreGen = nPreGen, maxGen = maxGen, 
                        fitnessValue = fitnessValue, dispersion = dispersion, evolutionPath = evolutionPath)
    #also add the default stopping conditions to the restart triggers
    #these need to be enabled at the same time as they perform some sanity checks
    restart_triggers = restart_triggers
  }
  if (!is.null(stopFitness) && !is.null(condition1)) {
    optValue = getGlobalOptimum(fun)$value
    condition2 = stopOnOptValue(optValue, stopFitness)
    if (OCD == FALSE) {
      result = rbga(popSize = dimension * 10, instance = instance, iters = maxit, evalFunc = fun, stringMin = rep(-5, dimension), 
                    stringMax = rep(5, dimension), control = list (stop.ons = c(list(condition1, condition2)), 
                                                                   max.restarts = max_restarts,
                                                                   restart.triggers = restart_triggers,
                                                                   restart.multiplier = restart_multiplier))
      
    }
    else {
      result = rbga(popSize = dimension * 10, instance = instance, iters = maxit, evalFunc = fun, stringMin = rep(-5, dimension), 
                    stringMax = rep(5, dimension),control = list (stop.ons = c(list(condition1, condition2, OCDcond)), 
                                                                  max.restarts = max_restarts,
                                                                  restart.triggers = restart_triggers,
                                                                  restart.multiplier = restart_multiplier))
    }
  }
  #if stop fitness is null
  else if (is.null(stopFitness) && !is.null(condition1)){
    if (OCD == FALSE) {
      result = rbga(popSize = dimension * 10, instance = instance, iters = maxit, evalFunc = fun, stringMin = rep(-5, dimension), 
                    stringMax = rep(5, dimension),control = list (stop.ons = list(condition1), 
                                                                  max.restarts = max_restarts,
                                                                  restart.triggers = restart_triggers,
                                                                  restart.multiplier = restart_multiplier))
    }
    else { 
      result = rbga(popSize = dimension * 10, iters = maxit, instance = instance, evalFunc = fun, stringMin = rep(-5, dimension), 
                    stringMax = rep(5, dimension),control = list (stop.ons = c(list(condition1, OCDcond)), 
                                                                  max.restarts = max_restarts,
                                                                  restart.triggers = restart_triggers,
                                                                  restart.multiplier = restart_multiplier))
    }
  }
  #use default if no stopping criterion is defined
  else stop("Enter at least one stopping condition!")
  return(result)
}

#' @rdname bbo-benchmarking
#' @importFrom snow makeCluster clusterCall clusterExport clusterApply stopCluster
#' @importFrom parallel detectCores
#' @import smoof
#' @export
#wrapper function for bbob_custom that parallelizes it
#disables the progressbar so check the output files to see how far the algorithm has gottenhh
bbob_custom_parallel = function(optimizer, algorithm_id, data_directory, dimensions = c(2, 3, 5, 10, 20, 40), 
                                instances = c(1:5, 41:50), function_ids = NULL, maxit = NULL, stopFitness = NULL, 
                                maxFE = NULL, debug.logging = FALSE, max_restarts = 0, 
                                restart_multiplier = 1, restart_triggers = character(0), OCD = FALSE, varLimit = NULL,
                                nPreGen = NULL, maxGen = NULL, fitnessValue = FALSE, dispersion = FALSE, evolutionPath = FALSE) {
  nCores = detectCores()
  cluster = makeCluster(nCores, type = "SOCK")
  #export relevant libraries + functions to the clusters
  clusterCall(cluster, function() require(cmaesr))
  clusterCall(cluster, function() require(bbob))
  #export all environment functions
  ex = Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
  clusterExport(cluster, ex)
  clusterApply(cl = cluster, x = function_ids, fun = function(x) bbob_custom(optimizer = optimizer, 
                                                                  algorithm_id = algorithm_id, 
                                                                  data_directory = data_directory, 
                                                                  dimensions = dimensions,
                                                                  instances = instances,
                                                                  function_ids = x,
                                                                  maxit = maxit,
                                                                  stopFitness = stopFitness,
                                                                  maxFE = maxFE,
                                                                  OCD = OCD, 
                                                                  debug.logging = debug.logging,
                                                                  max_restarts = max_restarts, 
                                                                  restart_multiplier = restart_multiplier, 
                                                                  restart_triggers = restart_triggers,
                                                                  varLimit = varLimit,
                                                                  nPreGen = nPreGen,
                                                                  maxGen = maxGen,
                                                                  fitnessValue = fitnessValue,
                                                                  dispersion = dispersion, 
                                                                  evolutionPath = evolutionPath))
  stopCluster(cluster)
}
