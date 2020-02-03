get_start_time <- function()
{
  starttime <- unclass(as.POSIXct(strptime(date(),"%c")))[1]
  return(starttime)
}