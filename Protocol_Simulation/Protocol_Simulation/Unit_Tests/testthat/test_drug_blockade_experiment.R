context("Testing drug_blockade_experiment.R")
dir.create("../../Output/test_drug_blockade_experiment_plots",
           showWarnings = FALSE)
test_that("Test that the code runs silently for the 5-second ramp protocol.",
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
            inputList <- read.csv("../../Input/voltage_program_data_files/5_second_ramp_protocol.csv",
                                  stringsAsFactors = FALSE)
            
            expect_silent(myVoltage <- voltage_step_program_list_mode(inputList))
            
            expect_silent(myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                          voltageProgram = myVoltage$voltage_program))
            
            initstates["V"] <- myVoltage$voltage_program[1]
            myEvents <- list(data = myEvents)
            
            
            myTemperature <- 37
            pars["T"] <- myTemperature
            

            mydrug <- "astemizole"
            expect_silent(myhERGparams <- load_hERG_fit_parameters(
              sprintf("../../../hERG_fitting/results/%s_boot_pars.csv",mydrug)
            ))
            new_hERG_pars <- myhERGparams[7,]
            pars["Kmax"] <- new_hERG_pars$Kmax
            pars["Ku"] <- new_hERG_pars$Ku
            pars["n"] <- new_hERG_pars$n
            pars["halfmax"] <- new_hERG_pars$halfmax
            pars["Vhalf"] <- new_hERG_pars$Vhalf
            
            pdf("../../Output/test_drug_blockade_experiment_plotstest_drug_blockade_experiment_5_second_ramp_protocol.pdf")
            expect_silent(
              blockadeData <- drug_blockade_experiment(modelname = "hergmod",
                                                   parameterVector = pars,
                                                   initialStates = initstates,
                                                   drugConc = 30.0,
                                                   voltage_program = myVoltage$voltage_program,
                                                   timepoints = myVoltage$timepoints,
                                                   breakTolerance = 0.001,
                                                   maxPulses = 100L,
                                                   minimumPulses = 10L,
                                                   mydrug = "astemizole",
                                                   myTemperature = 37.0,
                                                   diagnosticPlots = TRUE)
              
            )
            dev.off()
            expect_true(is.numeric(blockadeData$fractionBlockade))
            expect_lt(blockadeData$fractionBlockade,1.0 + 1e-3)
            expect_gt(blockadeData$fractionBlockade,0.0 - 1e-3)
            expect_true(is.numeric(blockadeData$timeOfMaxDrugCurrent))
            expect_true(is.numeric(blockadeData$timeOfRisingEdge))
            expect_true(is.numeric(blockadeData$adjustedTimeOfMaxDrugCurrent))
            expect_true(blockadeData$adjustedTimeOfMaxDrugCurrent <= blockadeData$timeOfMaxDrugCurrent)
          })

test_that("Test that the code runs silently for the 30-second ramp protocol.",
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
            inputList <- read.csv("../../Input/voltage_program_data_files/30_second_ramp_protocol.csv",
                                  stringsAsFactors = FALSE)
            
            expect_silent(myVoltage <- voltage_step_program_list_mode(inputList))
            
            expect_silent(myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                                        voltageProgram = myVoltage$voltage_program))
            
            initstates["V"] <- myVoltage$voltage_program[1]
            myEvents <- list(data = myEvents)

            myTemperature <- 37
            pars["T"] <- myTemperature
            mydrug <- "astemizole"
            expect_silent(myhERGparams <- load_hERG_fit_parameters(
              sprintf("../../../hERG_fitting/results/%s_boot_pars.csv",mydrug)
            ))
            new_hERG_pars <- myhERGparams[7,]
            pars["Kmax"] <- new_hERG_pars$Kmax
            pars["Ku"] <- new_hERG_pars$Ku
            pars["n"] <- new_hERG_pars$n
            pars["halfmax"] <- new_hERG_pars$halfmax
            pars["Vhalf"] <- new_hERG_pars$Vhalf
            
            pdf("../../Output/test_drug_blockade_experiment_plotstest_drug_blockade_experiment_30_second_ramp_protocol.pdf")
            expect_silent(
              blockadeData <- drug_blockade_experiment(modelname = "hergmod",
                                                   parameterVector = pars,
                                                   initialStates = initstates,
                                                   drugConc = 30.0,
                                                   voltage_program = myVoltage$voltage_program,
                                                   timepoints = myVoltage$timepoints,
                                                   breakTolerance = 0.001,
                                                   maxPulses = 100L,
                                                   minimumPulses = 10L,
                                                   mydrug = "astemizole",
                                                   myTemperature = 37.0,
                                                   diagnosticPlots = TRUE)
              
            )
            dev.off()
            expect_true(is.numeric(blockadeData$fractionBlockade))
            expect_lt(blockadeData$fractionBlockade,1.0 + 1e-3)
            expect_gt(blockadeData$fractionBlockade,0.0 - 1e-3)
            expect_true(is.numeric(blockadeData$timeOfMaxDrugCurrent))
            expect_true(is.numeric(blockadeData$timeOfRisingEdge))
            expect_true(is.numeric(blockadeData$adjustedTimeOfMaxDrugCurrent))
            expect_true(blockadeData$adjustedTimeOfMaxDrugCurrent <= blockadeData$timeOfMaxDrugCurrent)
          })

test_that("Test that the code runs silently for the experimental action potential.",
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
            inputList <- read.table("../../Input/voltage_program_data_files/AP2Hz",
                                  header = TRUE)

            expect_silent(myVoltage <- voltage_trueAP_program_list_mode(inputList))
            expect_silent(myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                                        voltageProgram = myVoltage$voltage_program))
            
            initstates["V"] <- myVoltage$voltage_program[1]
            myEvents <- list(data = myEvents)
            
            myTemperature <- 37.0
            pars["T"] <- myTemperature
            mydrug <- "astemizole"
            expect_silent(myhERGparams <- load_hERG_fit_parameters(
              sprintf("../../../hERG_fitting/results/%s_boot_pars.csv",mydrug)
            ))
            new_hERG_pars <- myhERGparams[7,]
            pars["Kmax"] <- new_hERG_pars$Kmax
            pars["Ku"] <- new_hERG_pars$Ku
            pars["n"] <- new_hERG_pars$n
            pars["halfmax"] <- new_hERG_pars$halfmax
            pars["Vhalf"] <- new_hERG_pars$Vhalf
            
            pdf("../../Output/test_drug_blockade_experiment_plotstest_drug_blockade_experiment_trueAP_protocol.pdf")
            expect_silent(
              blockadeData <- drug_blockade_experiment(modelname = "hergmod",
                                                   parameterVector = pars,
                                                   initialStates = initstates,
                                                   drugConc = 30.0,
                                                   voltage_program = myVoltage$voltage_program,
                                                   timepoints = myVoltage$timepoints,
                                                   breakTolerance = 0.001,
                                                   maxPulses = 100L,
                                                   minimumPulses = 10L,
                                                   mydrug = "astemizole",
                                                   myTemperature = 37.0,
                                                   diagnosticPlots = TRUE)
              
            )
            dev.off()
            expect_true(is.numeric(blockadeData$fractionBlockade))
            expect_lt(blockadeData$fractionBlockade,1.0 + 1e-3)
            expect_gt(blockadeData$fractionBlockade,0.0 - 1e-3)
            expect_true(is.numeric(blockadeData$timeOfMaxDrugCurrent))
            expect_true(is.numeric(blockadeData$timeOfRisingEdge))
            expect_true(is.numeric(blockadeData$adjustedTimeOfMaxDrugCurrent))
            expect_true(blockadeData$adjustedTimeOfMaxDrugCurrent <= blockadeData$timeOfMaxDrugCurrent)
          })
