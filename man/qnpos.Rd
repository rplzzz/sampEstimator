% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pnetfn.R
\name{qnpos}
\alias{qnpos}
\title{Find a quantile of the distribution of positive measurements}
\usage{
qnpos(p, k, popprev, Npop, phi, eta)
}
\arguments{
\item{p}{Probability of the quantile being sought}

\item{k}{Sample size}

\item{popprev}{Population prevalence}

\item{Npop}{Total population}

\item{phi}{Test sensitivity}

\item{eta}{Test specificity}
}
\description{
For the probability density defined in \code{\link{pnpos}}, find a quantile.
The quantile is defined as the smallest \eqn{k_{pos}} such that 
\eqn{\sum_{i=0}^{k_{pos}} dnpos(k_{pos}) \ge p}, for some specified p.
}
\details{
Since we don't have a good way to calculate the sum in the definition, other
than by iterating over the terms, we don't bother with a secant search or 
anything like that; we just keep summing until we get there.
}
