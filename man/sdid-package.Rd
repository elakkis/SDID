\name{sdid-package}
\alias{sdid-package}
\alias{sdid}
\docType{package}
\title{
Synthetic Difference in Difference}
\description{
The package implements Synthetic Difference in Difference algorithm Proposed by Arkhangelsky et al (2019)}
\details{
Package includes functions (all further information is available with help option)
simulate_panel(...)
prepare_data(...)
estimate(...)
att(...)
inference(...)
admm_solver(...)
admm_time_solver(...)
}
\author{
Emil Lakkis
Maintainer: Emil Lakkis <elakkis@umich.edu>
}
\references{
Arkhangelsky et al (2019)
}

\examples{

#simulate the data
N <- 100
T <- 40
treat_effect <- 1
treated_units <- 5
treat_date <- 20
covariates <- NULL
df <- simulate_panel(N, T, treat_effect, treated_units, treat_date)

#estimate
result <- estimate(df, covariates, treatment_time, treated_units)

}

\keyword{ package }
