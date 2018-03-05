elem <-
function(data, sym = NA)
#compute permutation matrix

#argument: data
#input: vector or time series object of numeric values

#argument: sym
#input: symbol used to represent missing values
#       can be numeric or character
#defaults to NA

#output: permutation matrix that will seperate the
#        observed and missing data
{
  if(is.na(sym)){
    oidx <- which(!is.na(data))  #observed data index
  }
  else{
    oidx <- which(data!=sym)
  }
  #get order of observed and missing values
  N <- length(data)
  permu <- c(oidx, (1:N)[-oidx])
  pmat <- cbind(1:N, permu)
  P <- matrix(0, N, N)
  P[pmat] <- 1
  return(P)
}
