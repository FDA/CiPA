drug_blockade_experiment <- function(modelname,
                                     parameterVector,
                                     initialStates,
                                     drugConc,
                                     voltage_program,
									                   timepoints,
                                     breakTolerance = 0.001,
                                     maxPulses = 100L,
                                     minimumPulses = 10L,
                                     mydrug,
                                     myTemperature,
									                   diagnosticPlots = FALSE)
{
  myVoltage <- voltage_program
  
  myEvents <- generate_events_voltage_dataframe(t.time = timepoints,
                                                voltageProgram = myVoltage)
  
  initialStates["V"] <- myVoltage[1]
  myEvents <- list(data = myEvents)
  
  
  #CONTROL PULSES
  initialStates["D"] <- 0 #Control pulses use no drug
  
  controlData <- blockade_calculation(maxPulses = maxPulses,
                                      minimumPulses = minimumPulses,
                                      breakTolerance = breakTolerance,
                                      modelname = modelname,
                                      initstates = initialStates,
                                      pars = parameterVector,
                                      timepoints = timepoints,
                                      events = myEvents,
                                      voltage_program = voltage_program)
  
  controlStateTrajectory <- controlData$stateTrajectory
  
  
  initialStates <- controlData$finalStateVector
  initialStates["D"] <- drugConc
  
  drugBlockadeData <- blockade_calculation(maxPulses = maxPulses,
                                           minimumPulses = minimumPulses,
                                           breakTolerance = breakTolerance,
                                           modelname = modelname,
                                           initstates = initialStates,
                                           pars = parameterVector,
                                           timepoints = timepoints,
                                           events = myEvents,
                                           voltage_program = voltage_program)
  
  drugBlockadeStateTrajectory <- drugBlockadeData$stateTrajectory

  const.Faraday <- 96485 #Couloumbs per mol.
  TKelvin <- myTemperature + 273.15 #Kelvin
  const.UGC <- 8.314 #J/mol/Kelvin
  Ki <- 150 #Units unknown.
  Ko <-5 #Units unknown.
  voltage.reversal <-((const.UGC*TKelvin)/const.Faraday)*log(Ko/Ki) * 1000 #convert to millivolts

  current.control <- controlStateTrajectory$O * (controlStateTrajectory$V - voltage.reversal)
  I.control.max <- max(current.control)

  current.drug <- drugBlockadeStateTrajectory$O * (drugBlockadeStateTrajectory$V - voltage.reversal)
  I.drug.max <- max(current.drug)

  outputList <- list()
  outputList$fractionBlockade <- 1 - I.drug.max/I.control.max
  outputList$timeOfMaxDrugCurrent <- drugBlockadeStateTrajectory$time[which.max(current.drug)]
  outputList$timeOfRisingEdge <- drugBlockadeStateTrajectory$time[which.max(diff(current.drug))]
  outputList$adjustedTimeOfMaxDrugCurrent <- outputList$timeOfMaxDrugCurrent - outputList$timeOfRisingEdge
  
  if(diagnosticPlots == TRUE)
  {
    plot(controlStateTrajectory$time,
         current.control / max(current.control),
         col = "blue",
         type = "l",
         xlab = "Time, t [ms]",
         ylab = "Normalized Current (I/Imaxcontrol) or Popen",
         main = sprintf("Drug = %s, ctrl converged in %d iterations, conc = %0.5g nM.",
                        mydrug,
                        controlData$convergedIn,
                        drugConc),
         sub = sprintf("Drug block converged in %d iterations.",
                       drugBlockadeData$convergedIn)
    )
    grid()
    
    lines(drugBlockadeStateTrajectory$time,
          current.drug / max(current.control), #We divide by max(current.control) to see how much reduction we obtained.
          col = "red")
    
    lines(controlStateTrajectory$time,
          controlStateTrajectory$O,
          col = "green")
    
    lines(controlStateTrajectory$time,
          drugBlockadeStateTrajectory$O,
          col = "green",
          lty = 2)
    
    legend("topright", legend=c("I, Final Control Pulse", "I, Final Drug Pulse", "Popen, Final Control Pulse", "Popen, Final Drug Pulse"),
           col=c("blue", "red","green"), lty=c(1,1,1,2), cex=0.75)
    
    string1 <- sprintf("fractionBlockade = %0.5g",outputList$fractionBlockade)
    string2 <- sprintf("timeOfMaxDrugCurrent = %0.5g",outputList$timeOfMaxDrugCurrent)
    string3 <- sprintf("timeOfRisingEdge = %0.5g",outputList$timeOfRisingEdge)
    string4 <- sprintf("adjustedTimeOfMaxDrugCurrent = %0.5g",outputList$adjustedTimeOfMaxDrugCurrent)
    
    myLabels <- c(string1,string2,string3,string4)
    text(x = 0.75 * max(controlStateTrajectory$time),
         y = c(0.40,0.30,0.20,0.10),
         labels = myLabels)
  }
  
  return(outputList)
}