covmat <-
function(phi, sigma2, N)
#compute AR(p) covariance matrix

#argumnet: phi
#input: AR coefficients, must satisfy stationarity constriants

#argument: sigma2
#input: innovation variance, must be greater than 0

#argument: N
#input: positive integer for data length

#output: gives the covariance matrix of AR(p)
{
  #if no AR order selected, then series is independent
  if(length(phi) == 0){
   Sigma <- sigma2 * diag(rep(1, N))
  }
  else{
   #ARMAacf gives correlation matrix
   Corr <- toeplitz(ARMAacf(ar = phi, lag.max = N - 1))
   gamma0 <- sigma2 / (1 - crossprod(phi, #unconditional variance
               Corr[1, 2:(1 + length(phi))]))
   Sigma <- as.numeric(gamma0) * Corr
  }
  return(Sigma)
}
