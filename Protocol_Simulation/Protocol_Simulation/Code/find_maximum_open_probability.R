find_maximum_open_probability <- function(inputDataFrame,
                                          voltage_program)
# tstart,
# tholding,
# trelax,
# deltaV,
# voltageOffset
{
  stopifnot(is.data.frame(inputDataFrame))
  stopifnot(is.numeric(voltage_program))
  #stopifnot(is.list(voltage_program_inputList))

  correctStateNames <- c("time",
                         "IC1",
                         "IC2",
                         "C1",
                         "C2",
                         "O",
                         "IO",
                         "IObound",
                         "Obound",
                         "Cbound",
                         "D",
                         "V")
  stopifnot(all_field_names_match(correctStateNames,names(inputDataFrame)))
  
  # tRampStart <- voltage_program_inputList$tstart + voltage_program_inputList$tholding
  # tRampEnd <- voltage_program_inputList$tstart + 
  #             voltage_program_inputList$tholding +
  #             voltage_program_inputList$trelax
  #print(tRampStart)
  #print(tRampEnd)
  #timesOfInterest <- inputDataFrame$time >= tRampStart & inputDataFrame$time <= tRampEnd
  
  #inputDataFrame <- inputDataFrame[timesOfInterest,]
  #print("P_OPEN:")
  #print(inputDataFrame$O)
  maximumOpenProbability <- max(inputDataFrame$O)
  #timeOfMaximum <- inputDataFrame$time[which.max(inputDataFrame$O)]
  indexOfMaximum_pOpen <- which.max(inputDataFrame$O)
  
  #voltageAtMaximum <- voltage_program[indexOfMaximum_pOpen]
  
  #Vmax_data <- voltage_step_program_list_mode(voltage_program_inputList)
  voltageAtMaximum <- voltage_program[indexOfMaximum_pOpen]
  
  stopifnot(is.numeric(maximumOpenProbability))
  stopifnot(is.numeric(voltageAtMaximum))

  outputValues <- list()
  outputValues$maximumOpenProbability <- maximumOpenProbability
  outputValues$voltageAtMaximum <- voltageAtMaximum
  return(outputValues)
}