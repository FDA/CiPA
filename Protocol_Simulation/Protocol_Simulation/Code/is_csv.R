is_csv <- function(myPath)
{
  stopifnot(is.character(myPath) == TRUE)
  iscsv <- grepl("\\.csv$",myPath)
  return(iscsv)
}