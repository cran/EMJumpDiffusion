`viterbi` <-
function(theta,X)
{

lambda<-theta[1]
mu<-theta[2]
sigma<-theta[3]
mud<-theta[4]
sigmad<-theta[5]
T<-length(X)


# initializing matrices and comparison parameters
comp <- -Inf
max <- -Inf
s <- 0
path <- c(1:T) # matrix to store the most most propable path
rem <- matrix(data=0,2,T) # matrix to remember the probabilities of the path
viterbi <- matrix(data=0,2,T) # matrix with Viterbi-probabilities
P <- matrix(data=0, 2) # transition matrix of jumps
F <- matrix(data=0, 2, T) # transition matrix of observed information


P[1] <- log(1-lambda)
P[2] <- log(lambda)

# calculate densities of observed log-returns
for (t in 1:T)
{
F[1,t] <- -0.5*log(2*pi*sigma^2) - (X[t]-mu+0.5*sigma^2)^2 /(2*sigma^2)
F[2,t] <- -0.5*log(2*pi*(sigma^2+sigmad^2)) - (X[t]-mud-mu+0.5*sigma^2)^2 /(2*(sigma^2+sigmad^2))
}



# initial distribution
viterbi[1,1] <- P[1]+F[1,1]
viterbi[2,1] <- P[2]+F[2,1]

# calculates viterbi-probabilities and the matrix rem
for( t in 2:T)
{
for( i in 1:2)
{
# iterative calculation of the maximum 
for(j in 1:2)
{
# comparison with preseding values
comp <- viterbi[j,t-1]+P[i]
if(comp >= max)
{
max <- comp
rem[i,t-1] <- j
}
}
# the calculated maximum equals the next viterbi-probability
viterbi[i,t] <- max+F[i,t]
max <- -Inf# reset for next iteration
}
}

# comparison of the last viterbi-probabilities to calculate the most probable path
comp <- -Inf# reset of coparison parameter 
for(i in 1:2)
{
if(viterbi[i,T]>=comp)
{
comp <- viterbi[i,T]
path[T] <- i-1
}
}


s<-path[T]+1
# backtracing
for(t in T:2)
{
path[t-1] <- rem[s,t-1]-1# following the most probable path
s <- rem[s,t-1]# reference to the next most probable path
}

return(path)# returns most probable path
}

