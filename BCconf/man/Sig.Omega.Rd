% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConfounderPvalue.R
\name{Sig.Omega}
\alias{Sig.Omega}
\title{Determine the significance of Omega.ols}
\usage{
Sig.Omega(Corr.Object, Omega.ols = NULL)
}
\arguments{
\item{Corr.Object}{A list outputed by \code{Correction}}

\item{Omega.ols}{Tests the null hypothesis for this specified value of Omega.ols. The default is to use the estimate stored in \code{Corr.Object}}
}
\value{
A list \item{p.values}{A \code{length(ind.cov)} vector of P values} \item{proxy.n}{A proxy for the sample size}
}
\description{
This tests is the null hypothesis Omega.ols = 0. If you accept the null hypothesis DO NOT replace Omega.ols with 0 in your estimates for the main effect B. Omega.ols is still needed to account for the correlation amongst test statistics.
}
