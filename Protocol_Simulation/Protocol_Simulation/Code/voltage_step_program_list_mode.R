voltage_step_program_list_mode <- function(inputList)
{
  stopifnot(is.list(inputList))
  expected_list_fields <- c("tfinal",
                            "tstart",
                            "tholding",
                            "trelax",
                            "deltaV",
                            "voltageOffset")

  fields_match <- all_field_names_match(names(inputList),expected_list_fields)
  if(fields_match == FALSE)
  {
    #print("EXPECTED LIST FIELDS:")
    #print(expected_list_fields)
    #print("INPUT LIST FIELDS:")
    #print(names(inputList))
    stop("At least one input field is missing or unexpected for the input voltage program and the input data file for that voltage program.")
  }

  #function(tinput,tstart,tholding,trelax,deltaV)
  tfinal <- inputList$tfinal
  tstart <- inputList$tstart
  tholding <- inputList$tholding
  trelax <- inputList$trelax
  deltaV <- inputList$deltaV
  voltageOffset <- inputList$voltageOffset
  
  if(tfinal < tstart + tholding + trelax)
  {
    stop("tfinal must be longer than the sum tstart + tholding + trelax.")
  }
  
  tinput <- 0:tfinal
  
  # varNameList <- c("tfinal","tstart","tholding","trelax","deltaV")
  # varList <- c(tinput,tstart,tholding,trelax,deltaV)
  # names(varList) <- varNameList

  Voutput <- NULL
  
  stopifnot(tinput >= 0)
  stopifnot(tstart >= 0)
  stopifnot(tholding >= 0)
  stopifnot(trelax >= 0)
  stopifnot(deltaV >= 0)
  
  for(t.time in tinput)
  {
    if(is_between(0,t.time,tstart))
    {
      V <- 0
    }
    
    if(is_between(tstart,t.time,tstart + tholding))
    {
     V <- deltaV
    }
    
    if(is_between(tstart + tholding,t.time,tstart + tholding + trelax))
    {
      V <- -(deltaV/trelax)*(t.time - (tstart + tholding)) + deltaV
    }
    
    if(t.time > tstart + tholding + trelax)
    {
      V <- 0
    }
    Voutput <- c(Voutput,V)
  }
  myList <- list()
  myList$timepoints <- 0:tfinal
  myList$voltage_program <- Voutput + voltageOffset
  stopifnot(all(is.numeric(myList$timepoints)))
  stopifnot(all(is.numeric(myList$voltage_program)))
  return(myList)
}