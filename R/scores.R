`scores` <-
function(org, calc)
{
hit<-0
miss<-0
total<-0
toomuch<-0
T<-length(org)
for (t in 2:(T-1))
{
if(org[t]==1)
{
total<-total+1
if(calc[t]==1)
{
hit<-hit+1
}
else
{
miss<-miss+1
}
}
else
{
if(calc[t]==1)
{
toomuch<-toomuch+1
}
}
} 
cat("Comparison of original jumps with estimated jumps \n")
cat("Total jumps:\t\t", total, "\n")
cat("Discovered jumps:\t", hit, "\n")
cat("Missed jumps:\t\t", miss, "\n")
cat("Additional jumps:\t", toomuch, "\n")
}

