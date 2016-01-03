#source("./cmaes/helpers.R")
#source("./cmaes/makeMonitor.R")
#source("./cmaes/makeStoppingCondition.R")
#source("./cmaes/stoppingConditions.R")




#' @title Covariance-Matrix-Adaption
#'
#' @description
#' Performs non-linear, non-convex optimization by means of the Covariance
#' Matrix Adaption - Evolution Strategy (CMA-ES).
#'
#' @details
#' This a pure R implementation of the popular CMA-ES optimizer for numeric
#' black box optimization [2, 3]. It features a flexible system of stopping conditions
#' and enables restarts [1], which can be triggered by arbitrary stopping conditions.
#'
#' You may pass additional parameters to the CMA-ES via the \code{control} argument.
#' This argument must be a named list. The following control elements will be considered
#' by the CMA-ES implementation:
#' \describe{
#'   \item{lambda [\code{integer(1)}]}{Number of offspring generaded in each generation.}
#'   \item{mu [\code{integer(1)}]}{Number of individuals in each population. Defaults to \eqn{\lfloor \lambda / 2\rfloor}.}
#'   \item{weights [\code{numeric}]}{Numeric vector of positive weights.}
#'   \item{sigma [\code{numeric(1)}]}{Initial step-size.}
#'   \item{restart.triggers [\code{character}]}{List of stopping condition codes / short names (see
#'   \code{\link{makeStoppingCondition}}). All stopping conditions which are placed in this vector do trigger a restart
#'   instead of leaving the main loop. Default is the empty character vector, i.e., restart is not triggered.}
#'   \item{max.restarts [\code{integer(1)}]}{Maximal number of restarts. Default is 0. If set
#'   to >= 1, the CMA-ES is restarted with a higher population size if one of the
#'   \code{restart.triggers} is activated.}
#'   \item{restart.multiplier [\code{numeric(1)}]}{Factor which is used to increase the population size after restart.}
#'   \item{stop.ons [\code{list}]}{List of stopping conditions. The default is to stop after 10 iterations or after a
#'   kind of a stagnation (see \code{\link{getDefaultStoppingConditions}})}.
#' }
#'
#' @references
#' [1] Auger and Hansen (2005). A Restart CMA Evolution Strategy With Increasing
#' Population Size. In IEEE Congress on Evolutionary Computation, CEC 2005, Proceedings,
#' pp. 1769-1776.
#' [2] N. Hansen (2006). The CMA Evolution Strategy: A Comparing Review. In J.A. Lozano,
#' P. Larranaga, I. Inza and E. Bengoetxea (Eds.). Towards a new evolutionary computation.
#' Advances in estimation of distribution algorithms. Springer, pp. 75-102.
#' [3] Hansen and Ostermeier (1996). Adapting arbitrary normal mutation distributions in evolution
#' strategies: The covariance matrix adaptation. In Proceedings of the 1996 IEEE
#' International Conference on Evolutionary Computation, pp. 312-317.
#'
#' @keywords optimize
#'
#' @param objective.fun [\code{smoof_function}]\cr
#'   Numerical objective function of type \code{smoof_function}. The function
#'   must expect a list of numerical values and return a scaler numerical value.
#' @param start.point [\code{numeric}]\cr
#'   Initial solution vector. If \code{NULL}, one is generated randomly within the
#'   box constraints offered by the paramter set of the objective function.
#'   Default is \code{NULL}.
#' @param monitor [\code{cma_monitor}]\cr
#'   Monitoring object.
#'   Default is \code{\link{makeSimpleMonitor}}.
#' @param control [\code{list}]\cr
#'   Futher paramters for the CMA-ES. See the details section for more in-depth
#'   information. Stopping conditions are also defined here.
#'   By default only some stopping conditions are passed. See \code{\link{getDefaultStoppingConditions}}.
#' @return [\code{CMAES_result}] Result object.
#'
#' @examples
#' # generate objective function from smoof package
#' fn = makeRosenbrockFunction(dimensions = 2L)
#' res = cmaes(
#'   fn,
#'   monitor = NULL,
#'   control = list(
#'     sigma = 1.5, lambda = 40,
#'     stop.ons = c(list(stopOnMaxIters(100L)), getDefaultStoppingConditions())
#'   )
#' )
#' print(res)
#'
#' @export
cmaes_custom = function(
  objective.fun,
  start.point = NULL,
  monitor = makeSimpleMonitor(),
  debug.logging = FALSE,
  control = list(
    stop.ons = c(
      getDefaultStoppingConditions()
    )
  )) {
  assertClass(objective.fun, "smoof_function")
  
  # extract relevant data
  par.set = getParamSet(objective.fun)
  lb = getLower(par.set); ub = getUpper(par.set)
  n = getNumberOfParameters(objective.fun)
  
  # sanity checks
  if (isNoisy(objective.fun)) {
    stopf("Noisy optimization is not supported at the moment.")
  }
  
  if (!isNumeric(par.set, include.int = FALSE)) {
    stopf("CMA-ES only works for objective functions with numeric parameters.")
  }
  
  if (isMultiobjective(objective.fun)) {
    stopf("CMA-ES can only handle single-objective functions.")
  }

  
  if (!is.null(monitor)) {
    assertClass(monitor, "cma_monitor")
  }
  
  # get stopping conditions
  stop.ons = getCMAESParameter(control, "stop.ons", NULL)
  if (is.null(stop.ons)) {
    stopf("There must be at least one stopping condition!")
  }
  assertList(stop.ons, min.len = 1L, types = "cma_stopping_condition")
  if (!is.null(start.point)) {
    assertNumeric(start.point, len = n, any.missing = FALSE)
  } else {
    if (!hasFiniteBoxConstraints(par.set)) {
      stopf("No start point provided. Cannot generate one, because parameter set cannot sample with Inf bounds!")
    }
    start.point = unlist(sampleValue(par.set))
  }
  # set initial distribution mean
  m = start.point
  
  # restart mechanism (IPOP-CMA-ES)
  restart.triggers = getCMAESParameter(control, "restart.triggers", character(0L))
  stop.ons.names = sapply(stop.ons, function(stop.on) stop.on$code)
  if(debug.logging == TRUE) print(collapse(stop.ons.names))
  if(debug.logging == TRUE) print(collapse(restart.triggers))
  if (!isSubset(restart.triggers, stop.ons.names)) {
    stopf("Only codes / short names of active stopping conditions allowed as restart trigger, but '%s' are no stopping conditions.", collapse(setdiff(restart.triggers, stop.ons.names), sep = ", "))
  }
  restart.multiplier = getCMAESParameter(control, "restart.multiplier", 2)
  assertNumber(restart.multiplier, lower = 1, na.ok = FALSE, finite = TRUE)
  max.restarts = getCMAESParameter(control, "max.restarts", 0L)
  assertInt(max.restarts)
  
  #FIXME: default value should be derived from bounds
  sigma = getCMAESParameter(control, "sigma", 0.5)
  assertNumber(sigma, lower = 0L, finite = TRUE)
  
  # Precompute E||N(0,I)||
  chi.n = sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n^2))
  
  # bookkeep best individual
  best.param = rep(NA, n)
  best.fitness = Inf
  
  # ======================================== added ====================================
  worst.fitness = Inf
  # ======================================== added ====================================

  # logs
  population.trace = list()
  
  # ======================================== added ====================================
  generation.bestfitness = list()
  evolutionPath = list()
  dispersion = list()
  # ======================================== added ====================================
  
  # init some termination criteria stuff
  iter = 0L
  n.evals = 0L
  start.time = Sys.time()
  
  # initialize stopped.on.t and stopped.on.chi that indicate the type of test which caused the termination of cma-es.
  # the stopping condition "stopOnOCD" sets the corresponding variable to "1" if that specific test has been significant.
  stopped.on.t = 0
  stopped.on.chi = 0
  
  result = callMonitor(monitor, "before")
  
  callMonitor(monitor, "before")
  
  # somehow dirty trick to "really quit" if stopping condition is met and
  # now more restart should be triggered.
  do.terminate = FALSE
  restarts = -1
  for (run in 0:max.restarts) {
    restarts = restarts + 1
    # population and offspring size
    if (run == 0) {
      lambda = getCMAESParameter(control, "lambda", 4L + floor(3 * log(n)))
      assertInt(lambda, lower = 4)
      mu = getCMAESParameter(control, "mu", floor(lambda / 2))
      assertInt(mu)
    } else {
      lambda = getCMAESParameter(control, "lambda", 4L + floor(3 * log(n)))
      # increase population size (IPOP-CMA-ES)
      lambda = ceiling(restart.multiplier^run * lambda)
      mu = floor(lambda / 2)
    }
    
    # path for covariance matrix C and stepsize sigma
    pc = rep(0, n)
    ps = rep(0, n)
    
    # initialize recombination weights
    weights = getCMAESParameter(control, "weights", log(mu + 0.5) - log(1:mu))
    if (any(weights < 0)) {
      stopf("All weights need to be positive, but there are %i negative ones.", sum(which(weights < 0)))
    }
    weights = weights / sum(weights)
    if (!(sum(weights) - 1.0) < .Machine$double.eps) {
      stopf("All 'weights' need to sum up to 1, but actually the sum is %f", sum(weights))
    }
    
    # variance-effectiveness / variance effective selection mass of sum w_i x_i
    mu.eff = sum(weights)^2 / sum(weights^2) # chosen such that mu.eff ~ lambda/4
    
    # step-size control
    cs = (mu.eff + 2) / (n + mu.eff + 5)
    ds = 1 + 2 * max(0, sqrt((mu.eff - 1) / (n + 1)) - 1) + cs # damping factor
    
    # covariance matrix adaption parameters
    cc = (4 + mu.eff / n) / (n + 4 + 2 * mu.eff / n)
    c1 = 2 / ((n + 1.3)^2 + mu.eff)
    alpha.mu = 2L
    cmu = min(1 - c1, alpha.mu * (mu.eff - 2 + 1/mu.eff) / ((n + 2)^2 + mu.eff))
    
    # covariance matrix
    sigma = getCMAESParameter(control, "sigma", 0.5)
    B = diag(n)
    D = diag(n)
    BD = B %*% D
    C = BD %*% t(BD) # C = B D^2 B^T = B B^T, since D equals I_n
    Cinvsqrt = B %*% diag(1 / sqrt(diag(D))) %*% t(B)
    
    
    # no restart trigger fired until now
    restarting = FALSE
    
    #restart iter logs the amount of iterations within the current restart loop (for OCD)
    restartIter = 0
    
    # break inner loop if terminating stopping condition active or
    # restart triggered
    while (!restarting) {
      iter = iter + 1L
      restartIter = restartIter + 1
      
      # create new population of search points
      arz = matrix(rnorm(n * lambda), ncol = lambda) # ~ N(0, I)
      ary = BD %*% arz # ~ N(0, C)
      arx = m + sigma * ary # ~ N(m, sigma^2 C)
      
      # Here we apply a penalization of violated bounds
      arx.repaired = ifelse(arx < lb, lb, ifelse(arx > ub, ub, arx))
      
      # Prepare penalization based on distance to repaired points (see Eq. 51)
      penalty.alpha = 1L
      penalties = penalty.alpha * colSums((arx - arx.repaired)^2)
      penalties[is.infinite(penalties)] = .Machine$double.max / 2
      
      # compute fitness values of repaired points
      fitn.repaired = if (isVectorized(objective.fun)) {
        objective.fun(arx.repaired)
      } else {
        apply(arx.repaired, 2L, function(x) objective.fun(x))
      }
      
      # apply penalization (see Eq. 51)
      fitn = fitn.repaired + penalties
      
      # update evaluation
      n.evals = n.evals + lambda
      
      # order fitness values
      fitn.ordered.idx = order(fitn, decreasing = FALSE)
      fitn.ordered = fitn[fitn.ordered.idx]
      
      # update best solution so far
      valid = (penalties == 0)
      if (any(valid)) {
        min.valid.idx = which.min(fitn.repaired[valid])
        if (fitn.repaired[valid][min.valid.idx] < best.fitness) {
          best.fitness = fitn.repaired[valid][min.valid.idx]
          best.param = arx[, valid, drop = FALSE][, min.valid.idx]
        }
      }
      
      # ======================================== added ====================================
      # update worst solution so far
      if (fitn.ordered[length(fitn.ordered)] > worst.fitness | is.infinite(worst.fitness)) {
        worst.fitness = fitn.ordered[length(fitn.ordered)]
      }
      # ======================================== added ====================================
      
      # update mean value / center of mass
      new.pop.idx = fitn.ordered.idx[1:mu]
      x.best = arx[, new.pop.idx, drop = FALSE]
      m.old = m
      m = drop(x.best %*% weights)
      
      y.best = ary[, new.pop.idx, drop = FALSE]
      y.w = drop(y.best %*% weights)
      z.best = arz[, new.pop.idx, drop = FALSE]
      z.w = drop(z.best %*% weights)
      
      # log population
      population.trace[[iter]] = x.best
      
      # Update evolution path with cumulative step-size adaption (CSA) / path length control
      # For an explanation of the last factor see appendix A in https://www.lri.fr/~hansen/cmatutorial.pdf
      ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mu.eff) * (Cinvsqrt %*% y.w)
      h.sigma = as.integer(norm2(ps) / sqrt(1 - (1 - cs)^(2 * (iter + 1))) < chi.n * (1.4 + 2 / (n + 1)))
      
      # Update covariance matrix
      pc = (1 - cc) * pc + h.sigma * sqrt(cc * (2 - cc) * mu.eff) * y.w
      y = BD %*% z.best
      delta.h.sigma = as.numeric((1 - h.sigma) * cc * (2 - cc) <= 1)
      C = (1 - c1 - cmu) * C + c1 * (pc %*% t(pc) + delta.h.sigma * C) + cmu * y %*% diag(weights) %*% t(y)
      
      # Update step-size sigma
      sigma = sigma * exp(cs / ds * ((norm2(ps) / chi.n) - 1))
      if (debug.logging) write(paste("Best fitness:", best.fitness), file = "debug.txt", append = TRUE)
      if (debug.logging) write(paste("sigma:", sigma), file = "debug.txt", append = TRUE)
      if (debug.logging) write(paste("Covmat:", sum(abs(C))), file = "debug.txt", append = TRUE)
      if (debug.logging) write(paste("Dispersion:", sum(abs(m-arx.repaired))), file = "debug.txt", append = TRUE)
      
      # Finally do decomposition C = B D^2 B^T
      e = eigen(C, symmetric = TRUE)
      B = e$vectors
      D = diag(sqrt(e$values))
      BD = B %*% D
      Cinvsqrt = B %*% diag(1 / diag(D)) %*% t(B) # update C^-1/2
      
      result = c(result, callMonitor(monitor, "step"))
      
      # escape flat fitness values
      if (fitn.ordered[1L] == fitn.ordered[ceiling(0.7 * lambda)]) {
        sigma = sigma * exp(0.2 + cs / ds)
        if (!is.null(monitor)) {
          warningf("Flat fitness values; increasing mutation step-size. Consider reformulating the objective!")
        }
      }
      
      ######### normalization and logging functionality for OCD ########
      if ("OCD" %in% stop.ons.names) {
        # get the call parameters from OCD needed for normalization
        param.set = stop.ons[[grep("OCD",stop.ons)]]$param.set
        # log best fitness value per generation
        if (!is.infinite(best.fitness)) generation.bestfitness[[iter]] = best.fitness
        else {
          print("infinite value caught")
          generation.bestfitness[[iter]] = .Machine$integer.max
        }
        # define upper and lower bound for normalization after nPreGen generations.
        # bounds are fixed once nPreGen generations are reached.
        if(iter == param.set[[2]]){
          upper.bound = worst.fitness
          lower.bound = min(unlist(generation.bestfitness))
        }
        # initialize list "evolutionPath" to be used as a performance indicator
        evolutionPath[[iter]] = sigma
        # initialize list "dispersion" to be used as a performance indicator
        dispersion[[iter]] = sum(abs(m-arx.repaired))
        # populate list with performance indicators for OCD.
        # if necessary, normalize the performance indicator of interest. 
        # For example, fitnessValue is normalized in the range upper.bound - lower.bound, 
        # i.e. the range of the objective values after nPreGen generations as defined above. This value is fixed for all upcomming generations
        performance.indicator = list("fitnessValue.iter" = if(iter < param.set[[2]]) best.fitness else (best.fitness)/(upper.bound-lower.bound),
                                    "fitnessValue.all" = if(iter < param.set[[2]]) generation.bestfitness[-length(generation.bestfitness)]
                                    else unlist(generation.bestfitness[-length(generation.bestfitness)])/(upper.bound-lower.bound),
                                    "dispersion.iter" = unlist(dispersion[length(dispersion)])/(sqrt(n)), 
                                    "dispersion.all" = unlist(dispersion[-length(dispersion)])/(sqrt(n)),
                                    "evolutionPath.iter" = unlist(evolutionPath[length(evolutionPath)])/(sqrt(n)), 
                                    "evolutionPath.all" = unlist(evolutionPath[-length(evolutionPath)])/(sqrt(n)))
        
      }
      
      
      
      # CHECK STOPPING CONDITIONS
      # =========================
      stop.obj = checkStoppingConditions(stop.ons)
      
      n.stop.codes = length(stop.obj$codes)
      if (max.restarts > 0L && any(stop.obj$codes %in% restart.triggers)) {
        if (!is.null(monitor)) {
          #messagef("Restart trigger fired! Restarting!!!")
        }
        n.stop.codes = sum(!(stop.obj$codes %in% restart.triggers))
        restarting = TRUE
      }
      
      # check if CMA-ES should really quit, i.e., is there a stopping condition,
      # that is active and does not trigger a restart?
      if (!restarting && (n.stop.codes > 0L)) {
        do.terminate = TRUE
        break
      }
    }
    
    # really quit without more restarts
    if (do.terminate) {
      break
    }
  }
  
  result = c(result, callMonitor(monitor, "after"))
  result = c(result, paste("-1", restarts))
  
  # log the type of test that caused the termination of cma-es in the output data.
  # "-2" indicated the termination based on the chi-squared test, "-3" indicates the termination based on the t-test
  if ("OCD" %in% stop.ons.names) {
    result = c(result, paste("-2", stopped.on.t))
    result = c(result, paste("-3", stopped.on.chi))
  }
  
  if (length(result) > 0) return(result)
  else {
  makeS3Obj(
      par.set = par.set,
      best.param = best.param,
      best.fitness = best.fitness,
      n.evals = n.evals,
      past.time = as.integer(difftime(Sys.time(), start.time, units = "secs")),
      n.iters = iter - 1L,
      n.restarts = run,
      population.trace = population.trace,
      message = stop.obj$stop.msgs,
      classes = "cma_result"
    )
  }
}

#' @export
print.cma_result = function(x, ...) {
  best.param = list(x$best.param)
  names(best.param) = getParamIds(x$par.set)
  catf("Best parameter      : %s", paramValueToString(x$par.set, best.param))
  catf("Best fitness value  : %.6g", x$best.fitness)
  catf("Termination         : %s", x$message)
  catf("  #Iterations       : %i", x$n.iters)
  catf("  #Evaluations      : %i", x$n.evals)
  catf("  Time              : %i (secs)", x$past.time)
}