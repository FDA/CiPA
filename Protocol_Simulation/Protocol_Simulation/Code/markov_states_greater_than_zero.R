#Tests if the sum of the probability states in the Markov model are all positive to within a tolerance of negativeTol.

markov_states_greater_than_zero <- function(stateDataFrame,
                                            negativeTol = 1e-3)
{
  stopifnot(is.data.frame(stateDataFrame))
  stopifnot(is.numeric(negativeTol))
  stopifnot(length(negativeTol) == 1)
  stopifnot(negativeTol >= 0)
  
  markovStateFields <- c("IC1",
                         "IC2",
                         "C1",
                         "C2",
                         "O",
                         "IO",
                         "IObound",
                         "Obound",
                         "Cbound")
  allGreaterThanZero <- NULL
  for(markovState in markovStateFields)
  {
    allGreaterThanZero[markovState] <- all(stateDataFrame[markovState] > (0 - negativeTol))
  }
  return(all(allGreaterThanZero) == TRUE)
}