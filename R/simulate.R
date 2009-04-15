`simulate` <-
function(theta,T)
{
lambda<-theta[1]
mu<-theta[2]
sigma<-theta[3]
mud<-theta[4]
sigmad<-theta[5]

jump<-matrix(data=0,T);
jumpdiff<-matrix(data=0,T);
back<-matrix(data=0,T,2);
Abgleich<-0
for (t in 1:T)
{
bern<-rbern(1,lambda)
if(bern==0)
{
jumpdiff[t]<-(mu-0.5*sigma^2)+rnorm(1,0,sigma)
jump[t]<-0
}
else
{
jumpdiff[t]<-(mu-0.5*sigma^2)+rnorm(1,0,sigma)+rnorm(1,mud,sigmad)
jump[t]<-1
}
}

back[,1]<-jumpdiff
back[,2]<-jump
return(back)
}

