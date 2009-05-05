\name{EM-JumpDiffusion-package}
\alias{EM-JumpDiffusion-package}
\alias{EM-JumpDiffusion}
\docType{package}
\title{
	Expectation-Maximization-Algorithm for Jump Diffusion Processes
}
\description{
	Calculates parameters for jump diffusion processes via EM-algorithm.
	The jumps-times are considered to be Bernoulli distributed (jump probability \eqn{\lambda}{lambda})
	with normally distributed jump-sizes (mean \eqn{\mu_d}{mu_d} and variance \eqn{\sigma_d}{sigma_d}).
	The diffusion process uses drift \eqn{\mu}{mu} and volatility \eqn{\sigma}{sigmad}.
	Is \eqn{S_t} a given set of stock price data, then it's log-returns are considered to follow the distribution:
	\deqn{log{\frac{S_t}{S_{t-1}}}\sim \left(1-\lambda\right)\mathcal{N}\left(\mu-\frac{\sigma^2}{2},\sigma^2\right)+\lambda\,\mathcal{N}\left(\mu-\frac{\sigma^2}{2}+\mu_d,\sigma^2+\sigma_d^2\right)}{
	log(S_t/S_(t-1))~(1-lambda) N(mu-sigma/2,sigma^2)+lambda N(mu-sigma/2+mu_d,sigma^2+sigma_d^2)}.
}
\details{
\tabular{ll}{
Package: \tab EM-JumpDiffusion\cr
Type: \tab Package\cr
Version: \tab 1.4.1\cr
Date: \tab 2009-05-05\cr
License: \tab GPL (version 2 or later)
}
	Use \code{\link{simulate}} to create log-returns or use real data.
	
	Use \code{\link{EMjump}} to calculate to optimize parameters of the jump diffusion process.
	
	The Viterbi-algorithm \code{\link{viterbi}} will calculate the most probable jump-times, given a set of log-returns.
	
	When using simulated data, \code{\link{scores}} can compare the original path of jumps with the calculated path of \code{\link{viterbi}}.
}
\references{
	Matthias Graser (2009). Expectation-Maximization-Algorithmus fuer Sprung-Diffusionsmodelle. Diploma thesis, Stuttgart.
	
	Rama Cont, Peter Tankov (2003). Financial modelling with Jump Processes. Chapman & Hall / CRC Press.
	
	Olivier Cappe, Eric Moulines, Tobias Ryden (2005). Interference in Hidden Markov Models. Springer Series in Statistics.
}
\author{
	Matthias Graser, <matthias.graser@gmx.de>
}
\keyword{ package }
\examples{
	## Will use every function in this package


	## Data simulation
	theta <- c(0.2, 0.1, 0.05, -0.5, 0.7)
	X <- simulate(theta, 30) # simulates 30 log-returns and uses function rbern
	plot(X[,1],type="h",main="log returns",xlab="time",ylab="value")

	start<-c(0.3, -0.1, 0.1, 0.3, 1) # starting vector

	## Parameter calculation
	calc <- EMjump(start, X[,1])

	## Viterbi-algorithm calculation
	calcpath <- viterbi(calc, X[,1])
	
	## Comparison of original path with calculated path
	scores(X[,2],calcpath)
}
\keyword{documentation}