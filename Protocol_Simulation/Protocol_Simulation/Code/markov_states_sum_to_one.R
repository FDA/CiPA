#Verifies that the states of the Markov model sum to exactly 1 to within a tolerance of +/- sumTol.

markov_states_sum_to_one <- function(stateDataFrame,
                                     sumTol = 1e-3)
{
  stopifnot(is.data.frame(stateDataFrame))
  stopifnot(is.numeric(sumTol))
  stopifnot(length(sumTol) == 1)
  markovStateFields <- c("IC1",
                         "IC2",
                         "C1",
                         "C2",
                         "O",
                         "IO",
                         "IObound",
                         "Obound",
                         "Cbound")
  
  sumVector <- rowSums(stateDataFrame[,markovStateFields])
  
  summedToOne <- NULL
  for(mySum in sumVector)
  {
    testValue <- is_between(1 - sumTol,mySum,1 + sumTol)
    summedToOne <- c(summedToOne,testValue)
  }
  return(all(summedToOne == TRUE))
}