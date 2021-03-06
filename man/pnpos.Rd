% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pnetfn.R
\name{pnpos}
\alias{pnpos}
\title{Probability of \eqn{n} positive observations, corrected for imperfect testing}
\usage{
pnpos(k, p, npos, Npop, phi, eta)
}
\arguments{
\item{k}{Sample size}

\item{p}{Population prevalence}

\item{npos}{Target number of positive observations}

\item{Npop}{Total population}

\item{phi}{Test sensitivity}

\item{eta}{Test specificity}
}
\description{
Compute the probability that we observe \eqn{n} positive test results, given a
population prevalence, total population size, and sensitivity and specificity
for the test.
}
\details{
The probability of observing X positive results is given by the probability
mass function for the hypergeometric distribution, \eqn{\rho_H}, which is
calculated with \code{\link[stats:Hypergeometric]{dhyper}}.  The correction for test errors
is calculated with \code{\link{pnetfn}}.  You get \eqn{n} positives if there
were \eqn{m} positives in your sample, and the net false negatives were equal
to \eqn{m-n}.  Therefore,
\deqn{P(X = n) = \sum_{i=0}^{k} \rho_H(i, M, N-M, k) * P(NFN=i-n).}

We allow the user to specify an arbitrary population prevalence; it will be
rounded to integer population counts.
}
\section{TODO}{

This function should really be called \code{dnpos}
}

