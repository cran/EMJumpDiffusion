\name{viterbi}
\alias{viterbi}
\title{Viterbi-algorithm}
\description{
	A Viterbi-algorithm to calculate the most probable arisal of jumps given a vector of log-returns and a parameter vector.
}
\usage{
	viterbi(theta, X)
}
\arguments{
  \item{theta}{Parameter vector.}
  \item{X}{Vector of log-returns.}
}
\details{
	The vector for the starting parameters must have 5 entries. The first one denotes the Bernoulli jump probability.
	The second and third entry have to be drift and volatility (must be >0) of the diffusion process.
	The last two are the mean and variance (must be >0) of the jump sizes.
	Initial guesses should not be too far off, or else calculated probabilities might become numerically=0.
  
	X should be a vector of log-returns, because the theory for this program was delevoped using their unique characteristics.
}
\value{
	The most probable path of jumps for log-returns X under parameter vector theta is returned.
	A 0 denotes no jump, a 1 denotes a jump.
}
\references{
	Cappe, O., Moulines, E. and Ryden, T. (2005). Interference in Hidden Markov Models. Springer Series in Statistics.
}
\examples{
	## Data simulation
	theta <- c(0.2, 0.1, 0.05, -0.5, 0.7) # vector for data simulation
	X <- simulate(theta, 100) # simulates 100 log-returns

	test <- c(0.3, -0.1, 0.1, 0.3, 1) # vector to test with viterbi
	
	## Viterbi-algorithm calculation
	calcpath <- viterbi(test, X[,1])
	plot(calcpath, main="jumps", xlab="time", ylab="value")
}
\keyword{optimize}