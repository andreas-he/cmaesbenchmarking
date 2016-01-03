#' @title Stopping condition: maximal iterations.
#'
#' @description Stop on maximal number of iterations.
#'
#' @param max.iter [integer(1)]\cr
#'   Maximal number of iterations.
#'   Default is \code{100}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnMaxIters = function(max.iter = 100L) {
  assertInt(max.iter, na.ok = FALSE)
  force(max.iter)
  return(makeStoppingCondition(
    name = "maxIter",
    message = sprintf("MaxIter: reached maximal number of iterations/generations %i.", max.iter),
    stop.fun = function(envir = parent.frame()) {
      return(envir$iter > max.iter)
    }
  ))
}

#' @title Stopping condition: indefinite covariance matrix.
#'
#' @description Stop if covariance matrix is not positive definite anymore.
#'
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnIndefCovMat = function() {
  return(makeStoppingCondition(
    name = "indefCovMat",
    message = "Covariance matrix is not numerically positive definite.",
    stop.fun = function(envir = parent.frame()) {
      e.values = envir$e$values
      return(any(e.values <= sqrt(.Machine$double.eps) * abs(max(e.values))))
    }
  ))
}

#' @title Stopping condition: optimal params.
#'
#' @description Stop if euclidean distance of parameter is below
#' some tolerance value.
#'
#' @param opt.param [\code{numeric}]\cr
#'   Known optimal parameter settings.
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \eqn{1e^{-8}}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnOptParam = function(opt.param, tol = 1e-8) {
  assertNumeric(opt.param, any.missing = FALSE, all.missing = FALSE)
  assertNumber(tol, lower = 0, na.ok = FALSE, finite = TRUE)
  force(opt.param)
  force(tol)
  return(makeStoppingCondition(
    name = "optParamTol",
    message = sprintf("Optimal parameters approximated nicely (gap < %.2f).", tol),
    stop.fun = function(envir = parent.frame()) {
      return(sqrt(sum(envir$best.param - opt.param)^2) < tol)
    }
  ))
}

#' @title Stopping condition: optimal objective value.
#'
#' @description Stop if best solution is close to optimal objective value.
#'
#' @param opt.value [\code{numeric(1)}]\cr
#'   Known optimal objective function value.
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \eqn{1e^{-8}}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnOptValue = function(opt.value, tol = 1e-8) {
  assertNumber(opt.value, na.ok = FALSE)
  assertNumber(tol, lower = 0, na.ok = FALSE, finite = TRUE)
  force(opt.value)
  force(tol)
  return(makeStoppingCondition(
    name = "optValTol",
    message = sprintf("Optimal function value approximated nicely (gap < %.10f).", tol),
    stop.fun = function(envir = parent.frame()) {
      return(abs(envir$best.fitness - opt.value) < tol)
    }
  ))
}

#' @title Stopping condition: maximal time.
#'
#' @description Stop if maximal running time budget is reached.
#'
#' @param budget [\code{integer(1)}]\cr
#'   Time budget in seconds.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnTimeBudget = function(budget) {
  assertInt(budget, na.ok = FALSE, lower = 1L)
  force(budget)
  return(makeStoppingCondition(
    name = "timeBudget",
    message = sprintf("Time budget of %i [secs] reached.", budget),
    stop.fun = function(envir = parent.frame()) {
      return(difftime(Sys.time(), envir$start.time, units = "secs") > budget)
    }
  ))
}

#' @title Stopping condition: maximal funtion evaluations.
#'
#' @description Stop if maximal number of function evaluations is reached.
#'
#' @param max.evals [\code{integer(1)}]\cr
#'   Maximal number of allowed function evaluations.
#' @return [\code{cma_stopping_condition}]
#' @export
stopOnMaxEvals = function(max.evals) {
  assertInt(max.evals, na.ok = FALSE, lower = 1L)
  force(max.evals)
  return(makeStoppingCondition(
    name = "maxEvals",
    message = sprintf("Maximal number of %i function evaluations reached.", max.evals),
    stop.fun = function(envir = parent.frame()) {
      return(envir$n.evals >= max.evals)
    }
  ))
}

#' @title Stopping condition: low standard deviation.
#'
#' @description Stop if the standard deviation falls below a tolerance value
#' in all coordinates?
#'
#' @param tol [\code{integer(1)}]\cr
#'   Tolerance value.
#' @return [\code{cma_stopping_condition}]
#' @export
#FIXME: default value is 10^(-12) * sigma. Here we have no access to the sigma value.
stopOnTolX = function(tol = 10^(-12)) {
  assertInt(tol, na.ok = FALSE)
  force(tol)
  return(makeStoppingCondition(
    name = "tolX",
    message = sprintf("Standard deviation below tolerance in all coordinates."),
    stop.fun = function(envir = parent.frame()) {
      return(all(envir$D < tol) && all((envir$sigma * envir$p.c) < tol))
    }
  ))
}

#' @title Stopping condition: principal axis.
#'
#' @description Stop if addition of 0.1 * sigma in a principal axis
#' direction does not change mean value.
#'
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnNoEffectAxis = function() {
  return(makeStoppingCondition(
    name = "noEffectAxis",
    message = "Addition of 0.1 times sigma does not change mean value.",
    stop.fun = function(envir = parent.frame()) {
      ii = (envir$iter %% envir$n) + 1L
      ui = envir$e$vectors[, ii]
      if (any(envir$e$values[ii] < 0)) return(TRUE)
      lambdai = sqrt(envir$e$values[ii])
      m.old = envir$m.old
      return(sum((m.old - (m.old + 0.1 * envir$sigma * lambdai * ui))^2) < .Machine$double.eps)
    }
  ))
}

#' @title Stopping condition: standard deviation in coordinates.
#'
#' @description Stop if addition of 0.2 * standard deviations in any
#' coordinate does not change mean value.
#'
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnNoEffectCoord = function() {
  return(makeStoppingCondition(
    name = "noEffectCoord",
    message = "Addition of 0.2 times sigma in any coordinate does not change mean value.",
    stop.fun = function(envir = parent.frame()) {
      m.old = envir$m.old
      return(sum((m.old - (m.old + 0.2 * envir$sigma))^2) < .Machine$double.eps)
    }
  ))
}

#' @title Stopping condition: high condition number.
#'
#' @description Stop if condition number of covariance matrix exceeds
#' tolerance value.
#'
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \code{1e14}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnCondCov = function(tol = 1e14) {
  assertNumber(tol, na.ok = FALSE, lower = 0, finite = TRUE)
  force(tol)
  return(makeStoppingCondition(
    name = "conditionCov",
    message = sprintf("Condition number of covariance matrix exceeds %f", tol),
    stop.fun = function(envir = parent.frame()) {
    #C = covmat
    #catch invalid values for covmat
    if (any(is.na(envir$C) | is.nan(envir$C) | is.infinite(envir$C))) return(TRUE)
    else return(kappa(envir$C) > tol)

    }
  ))
}


#===============================================================================
#=============================Online Convergence Detection======================
#===============================================================================
#' @title Stopping Condition: Online Convergence Detection.
#'
#' @description Online Convergence Detection (OCD) is a technique for detecting convergence of an algorithm based on statistical testing [1].
#' A two sided t-test as well as a chi-squared variance test are performed and serve as stopping criterion if significant.
#' @details 
#' Basically, two different analyses are performed for detecting convergence.
#' A statistical chi-squared variance test is performed which checks whether the variance of a set of performance indicator values decreases 
#' below a predefined variance limit significantly. Additionally, a two sided t-test is performed in order to check if there is no significant 
#' linear trend of the performance indicator values. The significance level for both tests is fixed with alpha = 0.05. 
#' The algorithm execution is terminated if one of these conditions holds for the last i and second last (i - 1) generation. In this implementation,
#' the performance indicator of interest are: \code{fitnessValue, dispersion, evolutionPath}. A performance indicator value corresponds to the difference
#' between e.g. the best fitness value of the current generation and that of the last generation. Depending on the number \code{nPreGen}, a vector
#' \code{PI} of length \code{nPreGen} is computed internally, that stores those differences for each active performance indicator.
#' @references 
#' [1] Wagner and Trautmann (2009). OCD: Online Convergence Detection for Evolutionary Multi-Objective Algorithms Based on Statistical Testing.
#' In Lecture Notes in Computer Science, pp. 198-215.
#' @param varLimit 
#'   \code{cma_stopping_condition} specifies the variance limit passed to the chi-squared variance test. The null hypotheses of this test
#'   is: var(PI) >= VarLimit 
#' @param nPreGen 
#'   The number \code{nPreGen} specifies the number of preceding generations for which the performance indicator values should be computed.
#'   The statistical tests consider exactly \code{nPreGen} generation in their testing procedure.
#' @param maxGen
#'   The number of iteration that should be spent at the maximum. \code{maGen} is the upper iteration limit if neither the t-test nor the
#'   chi-squared variance test terminate the optimization. Default is \code{maGen = NULL}.
#' @param fitnessValue
#'   The logical parameter \code{fitnessValue} indicates if the best fitness value of a population should be used as a performance indicator or not.
#'   Default is \code{fitnessValue = TRUE}.
#' @param dispersion
#'   The logical parameter \code{dispersion} indicates if the dispersion of the population should be used as a performance indicator or not.
#'   Default is \code{dispersion = FALSE}.
#' @param evolutionPath
#'   The logical parameter \code{evolutionPath} indicates if the evolutionPath or the step size, i.e. a cmaes parameter that controls 
#'   the evolution of the population in the objective space, should be used as a performance indicator or not.
#'   Default is \code{evolutionPath = FALSE}.
#' @return \code{stopOnOCD} returns TRUE if the optimizer should terminate the execution or FALSE if not.
#' @family stopping conditions
#' @export
stopOnOCD = function(varLimit, nPreGen,maxGen = NULL, fitnessValue = TRUE, dispersion = FALSE, evolutionPath = FALSE)
{
  # Check if varLimit is a single numeric
  assertNumber(varLimit, na.ok = FALSE)
  # Check if nPreGen is a single integerish value
  assertInt(nPreGen, na.ok = FALSE)
  # initialize significane level alpha with default value 0.05
  alpha = 0.05
  # Check if maxGen is a single integerish value. If no value for maxGen is passed, set maxGen to "Inf"
  if(!is.null(maxGen)) {
    assertInt(maxGen, na.ok = FALSE)
  }else{
    maxGen = Inf
  }
  # initialize vector of reference measures for calculating the different performance indicators.
  PF_i = numeric()
  # initialize list of performance indicator values
  PI_all = list()
  # initialize list of performance indicator values that stores the values of the generations i-nPreGen to i-1.
  PI_current_gen = list()
  # initialize list of performance indicator values that stores the values of the generations i-(nPreGen+1) to i-2.
  PI_preceding_gen = list()
  # initialize list of indicators
  indicator = numeric()
  # fill list with enabled indicators. Example: fitnessValue.iter contains the best.fitness value of the current generation,
  # fitnessValue.all contains the fitness values of all generations except the current iteration. 
  if(fitnessValue == TRUE) indicator = c(indicator, "fitnessValue.iter", "fitnessValue.all")
  if(dispersion == TRUE) indicator = c(indicator, "dispersion.iter", "dispersion.all")
  if(evolutionPath == TRUE) indicator = c(indicator, "evolutionPath.iter", "evolutionPath.all")
  #print(indicator)
  # initialize p-values of Chi-squared variance test
  pvalue_current_gen_chi = list()
  pvalue_preceding_gen_chi = list()
  # initialize p-values of the t-test on the regression coefficient
  pvalue_current_gen_t = list()
  pvalue_preceding_gen_t = list()
  # return stopping condition being compatible with cma-es implementation by Jakob Bossek
  return(makeStoppingCondition(
    name = "OCD",
    message = sprintf("OCD successful: Variance limit %f", varLimit),
    param.set = list(varLimit, nPreGen),
    stop.fun = function(envir = parent.frame()) {
      # Check if the number of iterations exceeds the user-defined number of maxGen. If TRUE, stop cma-es
      if(envir$restartIter >= maxGen){
        return(envir$restartIter >= maxGen)
      }
      
      # Check if number of iterations is greater than user-defined nPreGen
      if(envir$restartIter > nPreGen){
        # PF_i is the indicator value of the current iteration (e.g. best fitness value, sigma, or the dispersion of the individuals)
        # PF_i is used as a reference value for calculating the indicator values (difference to this value) of the last nPreGen generations.
        for (i in seq(1,length(indicator),2)){
          PF_i = c(PF_i, get(indicator[i], envir$performance.indicator))
        }
        chi.test = FALSE
        t.test = FALSE
        for (i in seq(2, length(indicator), 2)){
          # PI_all is a vector with one entry for each generation, except the first generation.
          # PI_all stores the difference between the performance indicator values of the last nPreGen generations and the current generation i.
          PI_all[[(i/2)]] = sapply(get(indicator[i], envir$performance.indicator), function (x) abs(x-PF_i[[(i/2)]]), simplify = TRUE)
          # PI_current_gen is a subset of PI_all which stores the last nPreGen indicator values with respect to the current generation i.
          PI_current_gen[[i/2]] = PI_all[[i/2]][(envir$iter-nPreGen):(envir$iter -1)]
          # in the first test, PI_preceding_gen equals PI_current_gen as indexing in this iteration would be invalid
          if((envir$restartIter - nPreGen) <= 1){
            # PI_preceding_gen is a subset of PI_all which stores the last nPreGen indicator values with respect to the last generation i-1.
            PI_preceding_gen = PI_current_gen
          }else{
            PI_preceding_gen[[i/2]] =  PI_all[[i/2]][(envir$iter - (nPreGen+1)):(envir$iter - 2)]
          }
          # perform chi2 variance tests and return corresponding p-values for each active performance indicator
          pvalue_current_gen_chi[[i/2]] = pChi2(varLimit, PI_current_gen[[i/2]])
          pvalue_preceding_gen_chi[[i/2]] = pChi2(varLimit, PI_preceding_gen[[i/2]])
          # perform two-sided t-test and return corresponding p-values for each active performance indicator
          pvalue_current_gen_t[[i/2]] = pReg(PI_current_gen[[i/2]])
          pvalue_preceding_gen_t[[i/2]] = pReg(PI_preceding_gen[[i/2]])
        }
        # set chi.test to TRUE if the pvalue of the chi-squared variance test is below the significance level alpha, i.e. reject H0: var(PI) >= varLimit
        if (all(pvalue_current_gen_chi <= alpha && pvalue_preceding_gen_chi <= alpha)) chi.test = TRUE
        # set t.test to TRUE if the pvalue of the two sided t-test is above the significance level alpha, i.e. reject H0: beta=0
        if (all(pvalue_current_gen_t > alpha && pvalue_preceding_gen_t > alpha)) t.test = TRUE
        # return TRUE, i.e. stop cmaes exectuion, if p-value is below specified significance level alpha
        # log termination condition in cma_es
        if (chi.test == TRUE) envir$stopped.on.chi = envir$stopped.on.chi + 1
        if (t.test == TRUE) envir$stopped.on.t = envir$stopped.on.t + 1
        # return TRUE, i.e. stop cmaes exectuion, if p-value is below specified significance level alpha
        return (chi.test || t.test)
      }
      else{
        return(FALSE)
      }
    }
  ))
}


#' @export
pChi2 <- function (varLimit, PI) {
  # Determine degrees of freedom
  N = length(PI)-1
  # calculate test statistic
  Chi = (var(PI)*N)/varLimit
  # get p-value of corresponding chi2 test
  p = pchisq(Chi, N, lower.tail = TRUE)
  return (p)
}

#' @export
pReg <- function (PI) {
  # Determin degrees of freedom
  N = length(PI)-1
  # standardize PI
  if(sd(PI)==0) {
    return (1)
  }else{
    PI = (PI-mean(PI))/sd(PI)
  }
  # initialize X, i.e. a vector of the generations of PI
  X = seq(1,length(PI),1)
  # standardize X
  X = (X-mean(X))/sd(X)
  # linear regression without intercept
  beta = (solve(X%*%X))%*%(X%*%PI)
  # residuals
  residuals = PI - X*beta
  # mse
  mse = (residuals%*%residuals) / N
  # compute test statistic
  t = beta/sqrt(mse*solve(X%*%X))
  # look up t distribution for N degrees of freedom
  p_value = 2*pt(-abs(t), df=N)
  return (p_value)
}

