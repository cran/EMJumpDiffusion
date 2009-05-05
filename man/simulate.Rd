\name{simulate}
\alias{simulate}
\title{Data simulation for use with EM-algorithm.}
\description{
	Simulates data of log-returns for use with \code{\link{EMjump}}.
}
\usage{
simulate(theta, T)
}
\arguments{
  \item{theta}{Parameter vector for simulated data.}
  \item{T}{Length of data to be simulated.}
}
\details{
  	The vector for the starting parameters must have 5 entries. The first one denotes the Bernoulli jump probability.
  	The second and third entry have to be drift and volatility (must be >0) of the diffusion process.
	The last two are the mean and variance (must be >0) of the jump sizes.
	
	T can be any positive integer number.
}
\value{
	Returns a matrix of dimension T x 2.
	The first row corresponds to the actual log-returns.
	The second row corresponds to jumps made, 0 for no jump and 1 for a jump.
}

\examples{
	## Data simulation
	theta <- c(0.2, 0.1, 0.05, -0.5, 0.7)
	X <- simulate(theta, 30) # simulates 30 log-returns
	plot(X[,1],type="h",main="log returns",xlab="time",ylab="value")
}
\keyword{datagen}