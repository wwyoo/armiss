\name{ar1sim}
\alias{ar1sim}
\docType{data}
\title{
Simulated AR(1) time series with random missing observations
}
\description{
This is a simulated AR(1) time series. The true process mean is zero, the AR coefficient was set to 0.5 and the innovation variance is one. The length of this series is 365. A total of 165 observations were randomly sampled without replacement to denote as missing values, by replacing their values with NA. Hence we have a total of 200 observed data. 
}
\usage{data(ar1sim)}
\format{
 Time-Series of length 365 with numeric entries for observed values but
 NA for missing values.\\
 First 5 entries: NA NA 0.9289 NA 0.0119 ...
}

\examples{
data(ar1sim)
ar1sim[1:20]  #take a look at the data
}
\keyword{datasets}
