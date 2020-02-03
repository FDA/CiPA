context("Testing find_maximum_open_probability.R")

test_that("Test that correct open probability is calculated from example data frame.",
          {
            state.time    <- c(0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.)
            
            state.IC1     <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
            state.IC2     <- c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
            state.C1      <- c(0.2,0.3,0.5,0.2,0.1,0.0,0.4,0.7,0.9,0.2,0.4)
            state.C2      <- c(0.2,0.1,0.2,0.3,0.2,0.1,0.2,0.2,0.0,0.2,0.1)
            state.O       <- c(0.1,0.1,0.1,0.3,0.5,0.6,0.3,0.0,0.0,0.2,0.2)
            state.IO      <- c(0.1,0.1,0.1,0.1,0.0,0.1,0.0,0.0,0.0,0.1,0.0)
            state.IObound <- c(0.1,0.1,0.0,0.0,0.1,0.1,0.0,0.0,0.0,0.1,0.0)
            state.Obound  <- c(0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1,0.1)
            state.Cbound  <- c(0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1)
            state.D       <- rep(100,11)
            state.V       <- rep(-80,11) #V state is irrelevant because the correct V is calculated internally.
            
            inputDataFrame <- data.frame(time = state.time,
                                         IC1 = state.IC1,
                                         IC2 = state.IC2,
                                         C1 = state.C1,
                                         C2 = state.C2,
                                         O = state.O,
                                         IO = state.IO,
                                         IObound = state.IObound,
                                         Obound = state.Obound,
                                         Cbound = state.Cbound,
                                         D = state.D,
                                         V = state.V)
            inputList <- list()
            inputList$tfinal <- 10
            inputList$tstart <- 2
            inputList$tholding <- 2
            inputList$trelax <- 3
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            
            voltageData <- voltage_step_program_list_mode(inputList = inputList)
            
            testData <- find_maximum_open_probability(inputDataFrame = inputDataFrame,
                                                      voltage_program = voltageData$voltage_program)
            pOMax <- testData$maximumOpenProbability
            VatMax <- testData$voltageAtMaximum
            
            expect_equal(0.6,pOMax)
            expect_equal(0.0,VatMax)
          })

test_that("Assertions trigger.",
          {
            state.time    <- c(0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.)
            
            state.IC1     <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
            state.IC2     <- c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
            state.C1      <- c(0.2,0.3,0.5,0.2,0.1,0.0,0.4,0.7,0.9,0.2,0.4)
            state.C2      <- c(0.2,0.1,0.2,0.3,0.2,0.1,0.2,0.2,0.0,0.2,0.1)
            state.O       <- c(0.1,0.1,0.1,0.3,0.5,0.6,0.3,0.0,0.0,0.2,0.2)
            state.IO      <- c(0.1,0.1,0.1,0.1,0.0,0.1,0.0,0.0,0.0,0.1,0.0)
            state.IObound <- c(0.1,0.1,0.0,0.0,0.1,0.1,0.0,0.0,0.0,0.1,0.0)
            state.Obound  <- c(0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1,0.1)
            state.Cbound  <- c(0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1)
            state.D       <- rep(100,11)
            state.V       <- rep(-80,11) #V state is irrelevant because the correct V is calculated internally.
            

            inputList <- list()
            inputList$tfinal <- 75
            inputList$tstart <- 2
            inputList$tholding <- 2
            inputList$trelax <- 3
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80

            voltageData <- voltage_step_program_list_mode(inputList = inputList)
            
            expect_error(
              find_maximum_open_probability(inputDataFrame = "a",
                                            voltage_program = voltageData$voltage_program)
              )
            
            inputDataFrame <- data.frame(time = state.time,
                                         IC1 = state.IC1,
                                         IC2 = state.IC2,
                                         C1 = state.C1,
                                         C2 = state.C2,
                                         O = state.O,
                                         IO = state.IO,
                                         IObound = state.IObound,
                                         Obound = state.Obound,
                                         Cbound = state.Cbound,
                                         D = state.D,
                                         V = state.V)
            
            expect_error(
              find_maximum_open_probability(inputDataFrame = inputDataFrame,
                                            voltage_program = "a")
            )
            
          })