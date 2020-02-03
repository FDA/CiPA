#Ensures all of the field names from a dataframe or a vector with named elements have a 1-to-1 correspondence with a vector of requiredNames.

all_field_names_match <- function(requiredNames,actualNames)
{
  stopifnot(is.character(requiredNames))
  stopifnot(is.character(actualNames))
  mustMatch1 <- all(requiredNames %in% actualNames)
  mustMatch2 <- all(actualNames %in% requiredNames)
  didMatch <- mustMatch1 & mustMatch2
  return(didMatch)
}