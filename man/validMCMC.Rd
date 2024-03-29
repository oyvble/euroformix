% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validMCMC.R
\name{validMCMC}
\alias{validMCMC}
\title{validMCMC}
\usage{
validMCMC(mcmcfit, trace = TRUE, acf = TRUE)
}
\arguments{
\item{mcmcfit}{A object returned by contLikMCMC}

\item{trace}{Whether showing trace of samples.}

\item{acf}{Whether showing autocorrelation function of samples.}
}
\description{
Validates aposteriori samples from MCMC method
}
\details{
This function takes samples from the MCMC as given in a matrix and shows the aposterior functions.
}
\author{
Oyvind Bleka
}
