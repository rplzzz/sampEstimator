
#' Calculate the probability of a specified net false negatives
#' 
#' The net false negatives are the number of false negatives minus the number
#' of false positives.  It is, in other words, the number of counts by which
#' the observed positives are shifted downward, relative to the actual number
#' of positive counts in the sample.  This function calculates the probability
#' that the net false negatives equal a specified value.
#' 
#' This probability can be calculated in terms of the probability mass function
#' for the binomial distribution, \eqn{\rho_{bn}} (\code{\link[stats]{dbinom}} in
#' R-lingo).
#' \deqn{P(NFN = \delta) = \sum_{j=0}^{k} \rho_{bn}(j, k, 1-\phi) \rho_{bn}(j-\delta, k, 1-\eta),}
#' where \eqn{\phi} is the sensitivity of the test and \eqn{\eta} is the specificity.
#' 
#' Note that this function is not vectorized.
#' 
#' @param k Number of trials
#' @param kpos Number of positive samples in the sample pool (before test error
#' is applied)
#' @param delta The number of net false negatives
#' @param phi Sensitivity (true positive rate) of the test.
#' @param eta Specificity (true negative rate) of the test.
#' @importFrom stats dbinom
#' @export
pnetfn <- function(k, kpos, delta, phi, eta)
{
  ## if either j or j+delta is outside of the range 0:k, then the term will be
  ## zero, so we can trim the range of the sum thus:
  j1 <- pmax(0, delta)
  
  ## A false positive can only trun a negative into a positive, and vice versa,
  ## so the number of trials is asymmetric.  This also limits the range of j
  kneg <- k-kpos
  j2 <- pmin(kneg+delta, kpos)
  
  if(j2 < j1) {
    ## no viable combination of false negatives and false positives
    return(0)
  }
  j <- seq(j1,j2)
  sum(
    ## j false negatives out of kpos potential positives &
    ## j-delta false positives out of kneg potential negatives
    dbinom(j, kpos, 1-phi) * dbinom(j-delta, kneg, 1-eta)
  )
}



#' Probability of \eqn{n} positive observations, corrected for imperfect testing
#' 
#' Compute the probability that we observe \eqn{n} positive test results, given a
#' population prevalence, total population size, and sensitivity and specificity
#' for the test.
#' 
#' The probability of observing X positive results is given by the probability
#' mass function for the hypergeometric distribution, \eqn{\rho_H}, which is
#' calculated with \code{\link[stats]{dhyper}}.  The correction for test errors
#' is calculated with \code{\link{pnetfn}}.  You get \eqn{n} positives if there
#' were \eqn{m} positives in your sample, and the net false negatives were equal
#' to \eqn{m-n}.  Therefore,
#' \deqn{P(X = n) = \sum_{i=0}^{k} \rho_H(i, M, N-M, k) * P(NFN=i-n).}
#' 
#' We allow the user to specify an arbitrary population prevalence; it will be
#' rounded to integer population counts.  
#' 
#' @param k Sample size
#' @param p Population prevalence
#' @param npos Target number of positive observations
#' @param Npop Total population
#' @param phi Test sensitivity
#' @param eta Test specificity
#' @importFrom stats dhyper
#' @export
pnpos <- function(k, p, npos, Npop, phi, eta)
{
  Mpos <- round(p*Npop)
  Mneg <- Npop - Mpos
  
  i <- seq(0,k)
  probihyper <- dhyper(i, Mpos, Mneg, k)
  probdelta <- sapply(i,
                      function(i) {
                        ## k samples, i potential positives, error shift of i-npos
                        pnetfn(k, i, i-npos, phi, eta)
                      })
  sum(probihyper * probdelta)
}
