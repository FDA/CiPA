#Runs the hERG model simulation of the Markov model by integrating the differential equations.
#modelname is from dyn.load() command from loading the .dll or .so compiled from hergmod.c.
#initstates is the initial values at t = 0 of the Markov model state vector.
#pars are from ../Input/hERG_parms/hergmod_pars.txt
#timepoints are the time vector, from t = 0 to t = tfinal for some desired tfinal.
#events is used with other functions to input an arbitrary voltage protocol.

run_simulation <- function(modelname,
                           initstates,
                           pars,
                           timepoints,
                           events = NULL)
  {
  stopifnot(all(is.numeric(pars)))
  
  stopifnot(is.character(modelname))
  expected_state_names <-  c("IC1",
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
  input_state_names <- names(initstates)
  stopifnot(all_field_names_match(requiredNames = expected_state_names,
                                  actualNames = input_state_names))

  stopifnot(markov_states_sum_to_one(data.frame(as.list(initstates))))
  stopifnot(markov_states_greater_than_zero(data.frame(as.list(initstates))))
  
  expected_par_names <- c("A1","B1","q1","A2","B2",
                          "q2","A3","B3","q3","A4",
                          "B4","q4","A11","B11","q11",
                          "A21","B21","q21","A31","B31",
                          "q31","A41","B41","q41","A51",
                          "B51","q51","A52","B52","q52",
                          "A53","B53","q53","A61","B61",
                          "q61","A62","B62","q62","A63",
                          "B63","q63","Kmax","Ku","n",
                          "halfmax","Kt","Vhalf","T","timeout",
                          "starttime")
  
  input_par_names <- names(pars)
  stopifnot(all_field_names_match(requiredNames = expected_par_names,
                                  actualNames = input_par_names))
  
  
  pars["starttime"] <- get_start_time()
  try(
    {
      outputStates <- ode(y = initstates,
                          times = timepoints,
                          func = "derivs",
                          parms = pars,
                          dllname = modelname,
                          initfunc = "initmod",
                          nout = 0,
                          rtol = 1e-3,
                          atol = 1e-6,
                          method = "lsoda",
                          events = events)
    }
  )
  
  doesNotExist <- !exists("outputStates")
  solverError <- inherits(outputStates,"try-error")
  wrongSolutionLength <- length(outputStates[,1]) != length(timepoints)
  incorrectTimePoints <- !all(outputStates[,1] == timepoints)
  nanValues <- any(is.nan(outputStates))
  
  generalFailure <- doesNotExist ||
    solverError ||
    wrongSolutionLength ||
    incorrectTimePoints ||
    nanValues
  
  if(generalFailure == TRUE)
  {
    return(NULL)
  }
  return(outputStates)
}