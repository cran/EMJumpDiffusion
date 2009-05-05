`EMjump` <-
function(theta,X,stop=10^-16)
{
lambda<-theta[1]
mu<-theta[2]
sigma<-theta[3]
mud<-theta[4]
sigmad<-theta[5]
T<-length(X)
if(stop==0)
{
stop<- 10^-16
}
it<-0

# parameter sigma needs to be calculated via optim
# this function represents the iterative calculation of the expected value of the log-likelihood within the EM-algorithm
# all constant addends have been removed

expval <- function (sigma)
{
expval <-0;
for(t in 1:T)
{
for(i in 0:1)
{
expval <- expval+ (phisum[i+1,t]*(-0.5*log((sigma^2)+i*(sigmad^2))-((X[t]-i*mud-mu+0.5*(sigma^2))^2 )/(2*((sigma^2)+i*(sigmad^2)))))
}
}
expval
}

# memory allocation for matrices and vectors
Q <- matrix(data=0, 2) # transition matrix of jumps
G <- matrix(data=0, 2, T) # transition matrix of observed information
c<-matrix(data=0,T) # scaling vector for forward-/backward-probabilities
A <-matrix(data=0, 2, T) # matrix for forward probabilities
Anorm <-matrix(data=0, 2, T) # matrix for scaled forward probabilities
B <- matrix(data=0, 2, T) # matrix for backward probabilities
Bnorm <- matrix(data=0, 2, T) # matrix for scalded backward probabilities
phi <- array(0,c(2,2,T)) # 3-dim matrix for probabilities phi(i,j,t)
phisum<-matrix(data=0,2,T) # matrix for probabilities phi(i,t)
done<-FALSE


# iteration of EM-algorithm
while(done==FALSE)
{

# vector to store the last parameter vector
old<-c(lambda,mu,sigma,mud,sigmad)

# initializing transition matrix of jumps
Q[1] <- 1-lambda
Q[2] <- lambda

# calculation of densities of observed log-returns
for (t in 1:T)
{
G[1,t] <- 0.001*1/(sqrt(2*pi*sigma^2))*exp(-(X[t]-mu+0.5*sigma^2)^2 /(2*sigma^2))
G[2,t] <- 0.001*1/(sqrt(2*pi*(sigma^2+sigmad^2)))*exp(-(X[t]-mud-mu+0.5*sigma^2)^2 /(2*(sigma^2+sigmad^2)))
}

# initializing and calculation of forward-probabilities

A[1,1]<-Q[1]*G[1,1]# Intialisierung der alpha_Tilde
A[2,1]<-Q[2]*G[2,1]# Intialisierung der alpha_Tilde
c[1]<-1/(A[1,1]+A[2,1]) # Initialisierung des Skalierungsvektors
Anorm[1,1]<-A[1,1]*c[1] # Initialisierung der alpha_Dach
Anorm[2,1]<-A[2,1]*c[1] # Initialisierung der alpha_Dach

# iteration to calculate forward-probabilities, their scaled versions and the scaling vector

for(t in 1:(T-1))
{
# calculation of forward-probabilities
A[1,t+1]<-Anorm[1,t]*Q[1]*G[1,t+1]+Anorm[2,t]*Q[1]*G[1,t+1]
A[2,t+1]<-Anorm[1,t]*Q[2]*G[2,t+1]+Anorm[2,t]*Q[2]*G[2,t+1]

# calculation of scaling vektor
c[t+1] <- 1/(sum(A[1:2,t+1]))

# calculation of scaled forward-probabilities
Anorm[1,t+1]<-c[t+1]*A[1,t+1]
Anorm[2,t+1]<-c[t+1]*A[2,t+1]
}

# initializing and calculation of backward-probabilities

B[1,T]<-1# Intialisierung der beta_Tilde
B[2,T]<-1# Intialisierung der beta_Tilde
Bnorm[1,T]<-c[T]*B[1,T]# Intialisierung der beta_Dach
Bnorm[2,T]<-c[T]*B[2,T] # Intialisierung der beta_Dach


# iteration to calculate backward-probabilities and their scaled versions

for(t in T:2)
{
# calculation of backward-probabilities
B[1,t-1]<-Q[1]*G[1,t]*Bnorm[1,t]+Q[2]*G[2,t]*Bnorm[2,t];
B[2,t-1]<-Q[1]*G[1,t]*Bnorm[1,t]+Q[2]*G[2,t]*Bnorm[2,t];

# calculation of scaled backward-probabilities
Bnorm[1,t-1]<-c[t-1]*B[1,t-1];
Bnorm[2,t-1]<-c[t-1]*B[2,t-1];
};

# calculating phi(i,j,t) and phi(i,t)

for(i in 1:2)
{
for(j in 1:2)
{
for(t in 1:(T-1))
{
phi[i,j,t]<-Anorm[i,t]*Q[j]*G[j,(t+1)]*Bnorm[j,(t+1)];
};
};
};


for(t in 1:(T-1))
{
phisum[1,t]<-phi[1,1,t]+phi[1,2,t];
phisum[2,t]<-phi[2,1,t]+phi[2,2,t];
};


# declaraton of additional variables to calculate the update equations
a1<-0
b1<-0
c1<-0
d1<-0

# update for lambda

for (t in 1:T)
{
a1<-a1+phi[1,1,t]
b1<-b1+phi[1,2,t]
c1<-c1+phi[2,1,t]
d1<-d1+phi[2,2,t]
}
lambda<-(b1+d1)/(a1+b1+c1+d1)

# reset
a1<-0
b1<-0
c1<-0
d1<-0
e1<-0

# calculation of updated mu-0.5*sigma^2 =: help1, mu_d and sigma^2+sigmad^2=:help2

for (t in 1:T)
{
a1<-a1+phisum[1,t]*X[t]# used for help1
b1<-b1+phisum[1,t]# used for help1
c1<-c1+phisum[2,t]*X[t]# used for mud
d1<-d1+phisum[2,t]# used for mud and help2
}

help1 <- a1/b1

for (t in 1:T)
{
e1<-e1+phisum[2,t]*(X[t]-help1-mud)^2# used for help2
}

mud<-c1/d1 - help1
help2<-e1/d1

# optimizing sigma from old parameters and new mud
# optim with method L-BFGS-B calculates its own gradient and works with
# upper and lower bounds. The lower bound prevents sigma from being equal to 0,
# the upper bound prevents sigma^2 from becoming greater than help2=sigma^2+sigmad^2
rtg<-optim(sigma,expval,NULL,method="L-BFGS-B",lower=0.0001,upper=sqrt(help2-0.0001),control=list( fnscale=-1))
sigma<-rtg$par[1]

# calculates sigmad and mu from help2 and help1
sigmad<-sqrt(help2-sigma^2)
mu<-help1+0.5*sigma^2

if( abs(lambda-old[1]) <= stop)
{
if( abs(mu-old[2]) <= stop)
{
if( abs(sigma-old[3]) <= stop)
{
if( abs(mud-old[4]) <= stop)
{
if( abs(sigmad-old[5]) <= stop)
{
if( it>100)
{
done<-TRUE
}
}
}
}
}
}
it<-it+1
}

cat("Estimated parameters of the jump diffusion process", "\n")
cat("Jump probability:\t", lambda,"\n")
cat("Drift:\t\t\t", mu,"\n")
cat("Volatility:\t\t", sigma,"\n")
cat("Mean of jumps:\t\t", mud,"\n")
cat("Variance of jumps:\t",sigmad,"\n")

theta<-c(lambda,mu,sigma,mud,sigmad)

return(theta)#returns optimized parameters
}

