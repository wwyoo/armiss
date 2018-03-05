\name{elem}
\alias{elem}
\title{
Permutation matrix
}
\description{
Given a time series with missing observations, this function will reorder its elements into the observed and missing parts respectively. 
}
\usage{
elem(data, sym = NA)
}
\arguments{
  \item{data}{
A vector or time series object with numeric entries. A list is also accepted if missing values are encoded with characters.
}
  \item{sym}{
Symbol used to represent missing values. Defaults to NA. Numeric and characters are allowed.
}
}

\value{
Returns the permutation matrix that will seperate observed and missing values in a time series.
}
\references{
Yoo, W. W. and Ghosh, S. K. (2013): "Gaussian Imputation of an ARMA
Process with Missing Observations", working paper.
}
\author{
William Weimin Yoo
}

\examples{
#example data with 8888 representing missing values
data <- c(1, 2, 8888, 4, 8888)
P <- elem(data = data, sym = 8888)
P %*% data
}

\keyword{Permutation matrix }
