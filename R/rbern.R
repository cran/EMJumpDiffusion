rbern <- function(n,lambda)
{
	try <- runif(n);
	print(try)
	for (i in 1:n)
	{
		if(lambda<try[i])
		{
			try[i]<-0
		}
		else
		{
			try[i]<-1
		}
	}
	return(try)
}