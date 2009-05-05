\name{rbern}
\alias{rbern}
\title{Bernoulli distribution}
\description{
	Creates Bernoulli distributed random variables via uniformly distrubuted random variables.
}
\usage{
	rbern(n, lambda)
}
\arguments{
  \item{n}{Number of observations.}
  \item{lambda}{Probability for observations.}
}
\details{
	The Bernoulli distribution generates a 0 with probability 1-lambda and a 1 with probability lambda.
	Since lambda is a probability, it must be between 0 and 1. 
}
\value{
	Vector of length n with values 0 and 1.
}
\examples{
	##Create a vector of 10 Bernoulli dristributed random variables with probability 0.2 for a 1
	u <- rbern(10,0.2)
}
\keyword{datagen}