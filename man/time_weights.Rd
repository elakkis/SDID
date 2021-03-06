% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{time_weights}
\alias{time_weights}
\title{The function find optimal time weights}
\usage{
time_weights(Y00, Y01, dzeta, rho, tol, weight_type)
}
\arguments{
\item{Y00}{- T0 x N0 matrix of control unit outcomes before the treatment}

\item{Y01}{- T-T0 x N0 matrix of control unit outcomes after the treatment}

\item{dzeta}{- Ridge penalty parameter}

\item{rho}{- parameter used in ADMM}

\item{tol}{- Convergence tolerance for ADMM algorithm}

\item{weight_type}{- type of weight to be computed. Supported options are "SDID" and "Kernel"}
}
\value{
T + 1 x 1 matrix of time weights for the estimation
}
\description{
The function find optimal time weights
}
