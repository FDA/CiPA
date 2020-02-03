#Just a shorter way of testing if LB <= x <= UB.
is_between <- function(LB,x,UB)
{
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(LB))
  stopifnot(is.numeric(UB))
  
  stopifnot(length(x) == 1)
  stopifnot(length(LB) == 1)
  stopifnot(length(UB) == 1)

  stopifnot(LB <= UB)
  
  if(x >= LB & x <= UB)
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}