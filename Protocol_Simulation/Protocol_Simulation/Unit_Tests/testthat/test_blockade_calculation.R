context("Testing blockade_calculation.R.")
test_that("Test that code runs with reasonable inputs for the control case.",
          {
            modelname <- "hergmod"
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            
            inputList <- list()
            inputList$tfinal <- 1500
            inputList$tstart <- 300
            inputList$tholding <- 500
            inputList$trelax <- 100
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            
            myVoltage <- voltage_step_program_list_mode(inputList)
            
            myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                          voltageProgram = myVoltage$voltage_program)
            
            initstates["V"] <- myVoltage$voltage_program[1]
            myEvents <- list(data = myEvents)
            
            
            myTemperature <- 37
            pars["T"] <- myTemperature
            
            mydrug <- "astemizole"
            drugConc <- 0 #nM of drug
            
            myhERGparams <- load_hERG_fit_parameters(
              sprintf("../../../hERG_fitting/results/%s_boot_pars.csv",mydrug)
            )
            new_hERG_pars <- myhERGparams[7,]
            pars["Kmax"] <- new_hERG_pars$Kmax
            pars["Ku"] <- new_hERG_pars$Ku
            pars["n"] <- new_hERG_pars$n
            pars["halfmax"] <- new_hERG_pars$halfmax
            pars["Vhalf"] <- new_hERG_pars$Vhalf
            
            
            # myOutputData <- blockade_calculation(maxPulses = 100L,
            #                                  minimumPulses = 10L,
            #                                  breakTolerance = 0.01,
            #                                  modelname = modelname,
            #                                  initstates = initstates,
            #                                  pars = pars,
            #                                  timepoints = timepoints,
            #                                  events = myEvents,
            #                                  tstart = tstart,
            #                                  tholding = tholding,
            #                                  trelax = trelax,
            #                                  deltaV = deltaV,
            #                                  voltageOffset = voltageOffset)
            
            myOutputData <- blockade_calculation(maxPulses = 100L,
                                                 minimumPulses = 10L,
                                                 breakTolerance = 0.01,
                                                 modelname = modelname,
                                                 initstates = initstates,
                                                 pars = pars,
                                                 timepoints = myVoltage$timepoints,
                                                 voltage_program = myVoltage$voltage_program,
                                                 events = myEvents)
            
            expect_true(is.list(myOutputData))
            expect_true(is.data.frame(myOutputData$stateTrajectory))
            expect_true(myOutputData$convergedIn >= 10L)
          })


test_that("Test that code runs with reasonable inputs when drug is added.",
          {
            modelname <- "hergmod"
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            
            inputList <- list()
            inputList$tfinal <- 1500
            inputList$tstart <- 300
            inputList$tholding <- 500
            inputList$trelax <- 100
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            
            myVoltage <- voltage_step_program_list_mode(inputList)
            
            myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                          voltageProgram = myVoltage$voltage_program)
            
            initstates["V"] <- myVoltage$voltage_program[1]
            myEvents <- list(data = myEvents)
            
            
            myTemperature <- 37
            pars["T"] <- myTemperature
            
            mydrug <- "astemizole"

            
            myhERGparams <- load_hERG_fit_parameters(
              sprintf("../../../hERG_fitting/results/%s_boot_pars.csv",mydrug)
            )
            new_hERG_pars <- myhERGparams[7,]
            pars["Kmax"] <- new_hERG_pars$Kmax
            pars["Ku"] <- new_hERG_pars$Ku
            pars["n"] <- new_hERG_pars$n
            pars["halfmax"] <- new_hERG_pars$halfmax
            pars["Vhalf"] <- new_hERG_pars$Vhalf
            
            drugConc <- 30 #nM of drug
            initstates["D"] <- drugConc
            myOutputData <- blockade_calculation(maxPulses = 100L,
                                                 minimumPulses = 10L,
                                                 breakTolerance = 0.01,
                                                 modelname = modelname,
                                                 initstates = initstates,
                                                 pars = pars,
                                                 timepoints = myVoltage$timepoints,
                                                 voltage_program = myVoltage$voltage_program,
                                                 events = myEvents)
            
            expect_true(is.list(myOutputData))
            expect_true(is.data.frame(myOutputData$stateTrajectory))
            expect_true(myOutputData$convergedIn >= 10L)
          })


test_that("Test that assertions trigger",
          {
            modelname <- "hergmod"
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            
            inputList <- list()
            inputList$tfinal <- 1500
            inputList$tstart <- 300
            inputList$tholding <- 500
            inputList$trelax <- 100
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            
            myVoltage <- voltage_step_program_list_mode(inputList)
            
            myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                          voltageProgram = myVoltage$voltage_program)
            
            initstates["V"] <- myVoltage$voltage_program[1]
            myEvents <- list(data = myEvents)
            
            
            myTemperature <- 37
            pars["T"] <- myTemperature
            
            mydrug <- "astemizole"
            
            
            myhERGparams <- load_hERG_fit_parameters(
              sprintf("../../../hERG_fitting/results/%s_boot_pars.csv",mydrug)
            )
            new_hERG_pars <- myhERGparams[7,]
            pars["Kmax"] <- new_hERG_pars$Kmax
            pars["Ku"] <- new_hERG_pars$Ku
            pars["n"] <- new_hERG_pars$n
            pars["halfmax"] <- new_hERG_pars$halfmax
            pars["Vhalf"] <- new_hERG_pars$Vhalf
            
            drugConc <- 30 #nM of drug
            initstates["D"] <- drugConc
           #max pulses is not an integer here.
           expect_error(blockade_calculation(maxPulses = 100.00000,
                                             minimumPulses = 10L,
                                             breakTolerance = 0.01,
                                             modelname = modelname,
                                             initstates = initstates,
                                             pars = pars,
                                             timepoints = myVoltage$timepoints,
                                             voltage_program = myVoltage$voltage_program,
                                             events = myEvents))
            
            #maxPulses is a negative integer here.
            expect_error(blockade_calculation(maxPulses = -100L,
                                              minimumPulses = 10L,
                                              breakTolerance = 0.01,
                                              modelname = modelname,
                                              initstates = initstates,
                                              pars = pars,
                                              timepoints = myVoltage$timepoints,
                                              voltage_program = myVoltage$voltage_program,
                                              events = myEvents))
            
            #min pulses is not an integer here
            expect_error(blockade_calculation(maxPulses = 100L,
                                              minimumPulses = 10.0000,
                                              breakTolerance = 0.01,
                                              modelname = modelname,
                                              initstates = initstates,
                                              pars = pars,
                                              timepoints = myVoltage$timepoints,
                                              voltage_program = myVoltage$voltage_program,
                                              events = myEvents))
            
            #minimum pulses is a negative integer here.
            expect_error(blockade_calculation(maxPulses = 100L,
                                              minimumPulses = -10L,
                                              breakTolerance = 0.01,
                                              modelname = modelname,
                                              initstates = initstates,
                                              pars = pars,
                                              timepoints = myVoltage$timepoints,
                                              voltage_program = myVoltage$voltage_program,
                                              events = myEvents))
            
            #breakTolerance is a string here.
            expect_error(blockade_calculation(maxPulses = 100L,
                                              minimumPulses = 10L,
                                              breakTolerance = "lol wrong input",
                                              modelname = modelname,
                                              initstates = initstates,
                                              pars = pars,
                                              timepoints = myVoltage$timepoints,
                                              voltage_program = myVoltage$voltage_program,
                                              events = myEvents)) 
          })