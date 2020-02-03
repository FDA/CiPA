#Computes the block for a given drug concentration.
#Pulses are input until the maximum open probability changes by 0.1% or less from one pulse to the next, or at maxPulses, whichever comes first.
#At least minimumPulses must occur.
#modelname is from dyn.load() of the .dll or .so compiled hergmod.c code.
#initstates is the initial values of the Markov model state vector given in ../Input/hERG_parms/hergmod_states.txt.
#pars is from ../Input/hERG_parms/hergmod_pars.txt
#timepoints is from t = 0 to t = tfinal for some desired tfinal.
#voltage_program is used to input an arbitrary voltage protocol. Used together with events.

blockade_calculation <- function(maxPulses = 50L,
                                 minimumPulses = 10L,
                                 breakTolerance = 0.001,
                                 modelname,
                                 initstates,
                                 pars,
                                 timepoints,
                                 voltage_program,
                                 events)
{
  stopifnot(breakTolerance > 0)
  stopifnot(is.numeric(breakTolerance))
  stopifnot(is.integer(maxPulses))
  stopifnot(maxPulses > 0)
  stopifnot(is.integer(minimumPulses))
  stopifnot(minimumPulses >= 0)
  
  priorMaxProbability <- 0
  for(i in 1:maxPulses)
  {
    myOutput <- run_simulation(modelname = modelname,
                               initstates = initstates,
                               pars = pars,
                               timepoints,
                               events = events)
    
    myOutput <- as.data.frame(myOutput)
    
    pData <- find_maximum_open_probability(inputDataFrame = myOutput,
                                           voltage_program = voltage_program)
    
    saveNames <- names(myOutput)
    initstates <- as.numeric(myOutput[nrow(myOutput),])
    names(initstates) <- saveNames
    initstates <- initstates[!(names(initstates) == "time")]
    
    relativeChange <- abs(priorMaxProbability - pData$maximumOpenProbability)/priorMaxProbability
    if(relativeChange <= breakTolerance && i >= minimumPulses)
    {
      break
    }
    else
    {
      priorMaxProbability <- pData$maximumOpenProbability
    }
  }
  
  if(markov_states_greater_than_zero(myOutput) == FALSE)
  {
    return(NULL)
  }
  
  if(markov_states_sum_to_one(myOutput) == FALSE)
  {
    return(NULL)
  }
  
  myOutputData <- list()
  myOutputData$stateTrajectory <- myOutput
  myOutputData$convergedIn <- i
  myOutputData$finalStateVector <- initstates
  myOutputData$maximumValues <- pData
  return(myOutputData)
}