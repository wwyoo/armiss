\name{covmat}
\alias{covmat}
\title{
Constructs the covariance matrix for an AR(p) process
}
\description{
Given the AR coefficients and the innovation variance, this function
constructs the corresponding covariance matrix for the process by
calling upon the \code{ARMAacf} function in the base \code{stats} package.
}
\usage{
covmat(phi, sigma2, N)
}
\arguments{
  \item{phi}{
A vector of AR coefficients, where they have to satisfy the necessary stationarity conditions.
}
  \item{sigma2}{
The innovation variance, must be greater then zero.
}
   \item{N}{
Positive integer for data length.
}
}
\details{
The method used can be found in Section 3.3 of Brockwell and Davis 1991.
}
\value{
Returns the covariance matrix of an AR(p) process
}
\references{
Brockwell, P.J. and Davis, R.A. (1991) Time Series: Theory and Methods, Second Edition, Springer.
}
\author{
William Weimin Yoo
}

\seealso{
\code{\link{ARMAacf}}
}
\examples{
ar <- c(0.5, 0.2)
sigma2 <- 1
Sigma <- covmat(phi = ar, sigma2 = sigma2, N = 5)
Sigma
}
\keyword{Autoregressive }
\keyword{Covariance }
