#A gamma distribution that can be shifted along the abscissa (x-axis).
#Otherwise is the same as an ordinary gamma distribution.
dshiftedgamma <- function(xdata, shapeParameter, scaleParameter, shiftParameter, log = FALSE)
{
	x <- xdata
	k <- shapeParameter
	theta <- scaleParameter
	x0 <- shiftParameter
	if(x - x0 <= 0)
	{
		print("negative logarithm")
		print(x[which(x - x0 <= 0)])
		print(x)
		stop()
	}
	if(log == FALSE)
	{
		f <- 1/(gamma(k)*theta^k)*(x - x0)^(k-1)*exp(-(x - x0)/theta)
	}
	else
	{
		f <- (k - 1)*log(x - x0) - log(gamma(k)) - k*log(theta) - (x - x0)/theta
	}
	return(f)
}