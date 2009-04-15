\name{EMjump}
\alias{EMjump}
\title{EM-algorithm for Jump Diffusion processes}
\description{
	Calculates parameters for jump diffusion processes via EM-algorithm.
	The jumps-times are considered to be Bernoulli distributed (jump probability \eqn{\lambda}{lambda})
	with normally distributed jump-sizes (mean \eqn{\mu_d}{mu_d} and variance \eqn{\sigma_d}{sigma_d}).
	The diffusion process uses drift \eqn{\mu}{mu} and volatility \eqn{\sigma}{sigmad}.
	Is \eqn{S_t} a given set of stock price data, then it's log-returns are considered to follow the distribution:
	\deqn{log{\frac{S_t}{S_{t-1}}}\sim \left(1-\lambda\right)\mathcal{N}\left(\mu-\frac{\sigma^2}{2},\sigma^2\right)+\lambda\,\mathcal{N}\left(\mu-\frac{\sigma^2}{2}+\mu_d,\sigma^2+\sigma_d^2\right)}{
	log(S_t/S_(t-1))~(1-lambda) N(mu-sigma/2,sigma^2)+lambda N(mu-sigma/2+mu_d,sigma^2+sigma_d^2)}.
}
\usage{
EMjump(theta, X, stop)
}
\arguments{
  \item{theta}{Vector for starting parameters.}
  \item{X}{Vector of log-returns.}
  \item{stop}{Stopping criteria}
}
\details{
	The vector for the starting parameters must have 5 entries.
	The first one denotes the Bernoulli jump probability.
	The second and third entry have to be drift and volatility (must be >0) of the diffusion process.
	The last two are the mean and variance (must be >0) of the jump sizes.
	Initial guesses should not be too far off, or else calculated probabilities might become numerically=0.
  
	X should be a vector of log-returns, because the theory for this program was delevoped using their unique characteristics.
  
	If stop is missing, a default value of 1e-16 will be used.
}
\value{
	A vector for the updated parameters will be returned in the same order as theta above.
}
\references{
	Matthias Graser (2009). Expectation-Maximization-Algorithmus fuer Sprung-Diffusionsmodelle. Diploma thesis, Stuttgart.
	
	Rama Cont, Peter Tankov (2003). Financial modelling with Jump Processes. Chapman & Hall / CRC Press.
	
	Olivier Cappe, Eric Moulines, Tobias Ryden (2005). Interference in Hidden Markov Models. Springer Series in Statistics.
}

\examples{
## Data simulation
theta <- c(0.2, 0.1, 0.05, -0.5, 0.7)
X <- simulate(theta, 30) # simulates 30 log-returns
plot(X[,1],type="h",main="log returns",xlab="time",ylab="value")

start <- c(0.3, -0.1, 0.1, 0.3, 1) # starting vector

## parameter calculation
EMjump(start, X[,1])
}
\keyword{optimize}