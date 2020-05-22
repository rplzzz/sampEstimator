#' Estimate the sample size required to achieve a specified precision
#' 
#' Calculate the number of samples needed to estimate the population prevalence
#' with a specified precision.  Precision is defined in terms of the quantiles
#' of the distribution of estimated prevalence values.  For example, we might say,
#' "If the true prevalence is 1%, the 5th percentile of estimated values must be
#' greater than 0.5%", which captures the idea that we want to be certain that 
#' we don't grossly underestimate the prevalence.
#' 
#' The characteristics of the system that are specified are:
#' \itemize{
#' \item{Total population size}
#' \item{Test characteristics}
#' \item{True population prevalence}
#' \item{Quantile to place a requirement on}
#' \item{Target value of that quantile}
#' }
#' 
#' The solver actually solves for the sample size that is as close as possible
#' to being equal to the required condition, so it is indifferent to whether 
#' we are setting an upper bound or a lower bound on the quantile.  So, in place
#' of the example above we could just as well say that we want the 95th percentile
#' of estimated values to be less than 2%.
#' 
#' Solution is carried out using a bounded secant solver.  No initial guess is
#' required because we can start with the basic assumption that the number of 
#' samples must be greater than zero and less than the total population.  Since
#' the number of samples must be a whole number, we shrink the interval until it's
#' less than a quarter-person wide, and then we return the ceiling of the upper
#' bound.  
#' 
#' @param Npop Population size
#' @param popprev True population prevalence
#' @param prob The quantile to place a restriction on (e.g., 0.05 for the 5th percentile)
#' @param targ Target value for \code{prob}
#' @param sensitivity Sensitivity of the test
#' @param specificity Specificity of the test
#' @export
sampEstimate <- function(Npop, popprev, prob, targ, sensitivity, specificity)
{
  ## Create the target function for the solver
  targfun <- function(nsamp) {
    nsamp <- round(nsamp)
    q <- qnpos(prob, nsamp, popprev, Npop, sensitivity, specificity) / nsamp
    q - targ
  }
  
  ## Start with some initial guesses.  These should easily bracket the true value
  x1 <- 5
  x2 <- sqrt(Npop)
  
  bssolv(x1, x2, targfun)  
}