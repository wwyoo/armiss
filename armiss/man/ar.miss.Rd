\name{ar.miss}
\alias{ar.miss}
\title{
Gaussian Imputation of an AR(p) Process with Missing Observations
}
\description{
The function estimates parameters of an AR(p) process using data with missing observations. It first reorders the data into the observed and missing parts respectively using a permutation matrix. Imputation is done by generating observations from a conditional Gaussian distribution based on the observed data. The reconstructed series is then estimated using maximum likelihood with the order chosen by AIC. 
}
\usage{
ar.miss(data, epsilon = 0.001, order = NULL, max.iter = 100, sym = NA,
  control.optim = list(maxit = 200))
}
\arguments{
  \item{data}{
    A vector or time series object with numeric entries, can be a list if missing values are encoded with characters
}
  \item{epsilon}{
Controls the rate of convergence in the imputation and estimation algorithm. Defaults to 0.001. Please see the \code{R} tutorial vignette for more details.
}
  \item{order}{
Positive integer order of the AR process. If none is specified (default), then order is selected using AIC.
}
  \item{max.iter}{
The upper limit iteration for the imputation and estimation algorithm. Defaults to 100. Please see the \code{R} tutorial vignette for more details.
}
  \item{sym}{
The symbol used to represent missing values. Defaults to NA. Can be a character string or numeric 
}
  \item{control.optim}{
Control variables for \code{optim} function used to do likelihood maximization. Here we set the BFGS iteration to have an upper limit of 200. See the help page on optim for more details.
}
}
\details{
Please refer to the tutorial/manual \code{R} vignette for more details regarding the theorectical underpinings of the algorithm.}
\value{
\item{mu }{MLE for AR(p) process mean}
\item{phi }{A vector of MLE for AR(p) coefficients, where the order p
  was chosen by AIC}
\item{sigma2}{MLE for innovation variance}
}
\references{
Yoo, W. W. and Ghosh, S. K. (2013): "Gaussian Imputation of an ARMA
Process with Missing Observations", working paper.
}
\author{
William Weimin Yoo
}
\note{
This is still work in progress. More functionality will be added in the future. 
}

\seealso{
\code{\link{ar.mle}}, \code{\link{arima}}, \code{\link{optim}}
}
\examples{
#generate AR(2) with missing values
#true mean = 0, AR coefficient = c(0.5, 0.2), innovation variance = 1
set.seed(1234)
N <- 365  #data length
ar2 <- arima.sim(n = N, list(ar = c(0.5, 0.2), sd = 1))
index <- sample(1:N, (N - 200), replace = FALSE)  #165 missing data
ar2[index] <- NA

est <- ar.miss(ar2)
est 
}
\keyword{Gaussian imputation }
\keyword{Autoregressive process }
