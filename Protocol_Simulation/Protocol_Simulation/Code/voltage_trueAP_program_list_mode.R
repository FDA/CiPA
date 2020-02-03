voltage_trueAP_program_list_mode <- function(inputList)
{
  stopifnot(is.list(inputList))
  
  expected_list_fields <- c("Time",
                            "V")
  
  fields_match <- all_field_names_match(names(inputList),expected_list_fields)
  if(fields_match == FALSE)
  {
    #print("EXPECTED LIST FIELDS:")
    #print(expected_list_fields)
    #print("INPUT LIST FIELDS:")
    #print(names(inputList))
    stop("At least one input field is missing or unexpected for the input voltage program and the input data file for that voltage program.")
  }
  
  #Error if the time vector is not strictly monotonically increasing.
  stopifnot(all(diff(inputList$Time) > 0))
  
  myList <- list()
  myList$timepoints <- inputList$Time * 1000 #convert to milliseconds.
  myList$voltage_program <- inputList$V
  
  #warning if the time vector does not appear to be re-scaled to milliseconds.
  if(max(myList$timepoints) < 50)
  {
    warning("Total length of input time is being read as under 50 milliseconds. Have the time data been converted to milliseconds?")
  }
  
  #warning if the time vector goes beyond 60,000 milliseconds (1 minute).
  if(max(myList$timepoints) > 60000)
  {
    warning("Total length of input time is beyond 60,000 milliseconds, or beyond 1 minute. This is probably too long. Has a mistake been made?")
  }
  
  if(max(abs(myList$voltage_program)) > 1000)
  {
   warning("Maximum absolute value of observed voltage in the action potential raw data exceeds 1000 millivolts, which does not appear biologically possible. Has a mistake been made?") 
  }
  
  if(max(abs(myList$voltage_program)) < 1.0)
  {
    warning("Maximum absolute value of observed voltage in the action potential raw data is under 1.0 millivolts, which is far below that expected of an action potential. Has a mistake been made?") 
  }
  
  stopifnot(all(is.numeric(myList$timepoints)))
  stopifnot(all(is.numeric(myList$voltage_program)))
  stopifnot(all(myList$timepoints >= 0))
  return(myList)
}