\name{scores}
\alias{scores}
\title{Compare real to calculated data}
\description{
	Compares calculated paths of \code{\link{viterbi}} with the original path.
}
\usage{
scores(org, calc)
}
\arguments{
  \item{org}{Original path of jumps.}
  \item{calc}{Calculated path of jumps via \code{\link{viterbi}}.}
}

\value{
	No data is returned. However, this function will tell you all necessary quantities.
}
\examples{
	## Data simulation
	theta <- c(0.2, 0.1, 0.05, -0.5, 0.7) # vector for data simulation
	X <- simulate(theta, 100) # simulates 100 log-returns

	test <- c(0.3, -0.1, 0.1, 0.3, 1) # vector to test with viterbi
	
	## Viterbi-algorithm calculation
	calcpath <- viterbi(test, X[,1])
	plot(calcpath, main="jumps", xlab="time", ylab="value")

	## Comparison of original path with calculated path
	scores(X[,2],calcpath)
}
\keyword{arith}