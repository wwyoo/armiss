\name{armiss}
\alias{armiss-package}
\alias{armiss}
\docType{package}
\title{
Estimating parameters for an AR(p) process using missing data via
Gaussian imputation.
}
\description{
Given a time series with missing observations, this package will reorder
the series by dividing the data into the observed and missing parts
respectively using a permutation matrix. The missing values are then
imputed using conditional Gaussian distribution based on observed data. The reconstructed series is then used for parameter estimation.
}
\details{
\tabular{ll}{
Package: \tab armiss\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-27\cr
License: \tab GPL(>=2)\cr
}
The package has three functions: \code{ar.miss}, \code{elem} and \code{covmat}. The function
\code{ar.miss} implements the imputation and estimation procedures, and calls
upon \code{elem} and \code{covmat} to construct the permutation matrix and AR(p)
covariance matrix respectively. The function \code{ar.miss} accepts data in the
form of vector or time series object with numeric values and outputs the
parameter estimates. For futher detail, please the documentation page
for \code{ar.miss}.
}
\author{
William Weimin Yoo

Maintainer: William Weimin Yoo <wyoo@ncsu.edu>

}
\references{
Yoo, W. W. and Ghosh, S. K. (2013): "Gaussian Imputation of an ARMA
Process with Missing Observations", working paper.
}
\keyword{Gaussian imputation }
\keyword{Autoregressive process }
\seealso{
\code{\link{ar.mle}}, \code{\link{arima}} 
}
\examples{
#generate AR(1) with missing values
#true mean = 0, AR coefficient = 0.5, innovation variance = 1
set.seed(2345)
N <- 365
ar1 <- arima.sim(n = N, list(ar = 0.5, sd = 1))
index <- sample(1:N, (N - 200), replace = FALSE)  #165 missing values
ar1[index] <- NA

est <- ar.miss(ar1)
est
}
