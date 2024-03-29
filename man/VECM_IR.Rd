% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VECM_IR.R
\name{VECM_IR}
\alias{VECM_IR}
\title{Impulse responses analysis}
\usage{
VECM_IR(
  object,
  horizon = 10,
  transform = c("none", "orthogonal"),
  method = c("none", "analytic")
)
}
\arguments{
\item{object}{An object of class \code{VECM}.}

\item{horizon}{The forecasting horizon.}

\item{transform}{...}

\item{method}{...}
}
\value{
An object of class \code{VECM_IR} is a list with the following components: \cr

\item{theta}{...}
\item{sigma_theta}{...}
\item{theta_cum}{...}
\item{sigma_theta_cum}{...}
\item{transform}{...}
\item{method}{...}
}
\description{
This function calculates impulse response functions for VEC models.
}
\examples{
#ciao

}
\references{
Lütkepohl, H., & Krätzig, M. (2004). \emph{Applied time series econometrics}. Cambridge: Cambridge University Press.

Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
}
\author{
Gianluca Carpigo \email{gianluca.carpigo@uniroma1.it}

Michele Cianfriglia \email{michele.cianfriglia@libero.it}
}
