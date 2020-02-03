#Log-logistic distribution density function.
dllogis <- function(xdata, shapeParameter, scaleParameter, log = FALSE)
{
	if(shapeParameter <= 0)
	{
		stop("shapeParameter must be > 0.")
	}
	
	if(scaleParameter <= 0)
	{
		stop("scaleParameter must be > 0.")
	}
	
	a <- scaleParameter
	b <- shapeParameter
	x <- xdata
	
	if(log == FALSE)
	{
		f <- (b/a)*(x/a)^(b-1)/(1+(x/a)^b)^2
		return(f)
	}
	else
	{
		f <- (log(b/a) + (b-1)*log(x/a)) - 2*log(1 + (x/a)^b)
		return(f)
	}
}