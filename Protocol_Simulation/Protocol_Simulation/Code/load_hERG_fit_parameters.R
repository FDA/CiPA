load_hERG_fit_parameters <- function(filePath)
{
  stopifnot(is.character(filePath) == TRUE)
  stopifnot(is_csv(filePath))
  #print(getwd())
  #print(filePath)
  myTable <- read.csv(filePath,stringsAsFactors = FALSE)
  
  requiredNames <- c("Kmax",
                     "Ku",
                     "n",
                     "halfmax",
                     "Vhalf",
                     "slope")
  
  actualNames <- names(myTable)
  stopifnot(all_field_names_match(requiredNames,actualNames) == TRUE)
  
  stopifnot(is.numeric(myTable$Kmax) == TRUE)
  stopifnot(is.numeric(myTable$Ku) == TRUE)
  stopifnot(is.numeric(myTable$n) == TRUE)
  stopifnot(is.numeric(myTable$halfmax) == TRUE)
  stopifnot(is.numeric(myTable$Vhalf) == TRUE)
  stopifnot(is.numeric(myTable$slope) == TRUE)
  
  return(myTable)
}