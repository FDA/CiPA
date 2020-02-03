#From a given voltage protocol, this function is used to convert the voltageProgram into an events data frame that will work with lsoda().

generate_events_voltage_dataframe <- function(t.time,voltageProgram)
{
  stopifnot(is.numeric(t.time))
  stopifnot(is.numeric(voltageProgram))
  stopifnot(all(t.time >= 0))
  N1 <- length(t.time)
  N2 <- length(voltageProgram)
  stopifnot(N1 == N2)
  
  eventData <- data.frame(var = rep("V",N1),
                          time = t.time,
                          value = voltageProgram,
                          method = rep("replace",N1))
  currentVoltage <- Inf
  
  myOutput <- data.frame(var = character(),
                         time = numeric(),
                         value = numeric(),
                         method = character())
  
  for(i in 1:nrow(eventData))
  {
    if(currentVoltage == eventData$value[i])
    {
      next
    }
    else
    {
      mySubframe <- data.frame(var = "V",
                               time = t.time[i],
                               value = eventData$value[i],
                               method = "replace")
      myOutput <- rbind(myOutput,mySubframe)
      currentVoltage <- eventData$value[i]
    }
  }
  eventData <- myOutput
  
  return(eventData)
}