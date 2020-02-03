context("Testing run_simulation.R.")
test_that("Simulation will run with default parameters.",
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
            timepoints <- 0:1000
            
            
            expect_silent(myOutput <- run_simulation(modelname,
                                                     initstates,
                                                     pars,
                                                     timepoints,
                                                     events = NULL))
            expect_false(is.null(myOutput))
            myOutputDataFrame <- as.data.frame(myOutput)
            expect_true(markov_states_sum_to_one(myOutputDataFrame,1e-6))
            expect_true(markov_states_greater_than_zero(myOutputDataFrame,1e-6))
          })

test_that("Simulation will run with a random row of hERG parameters for astemizole at 1x Cmax.",
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
            timepoints <- 0:1000
            astemizolehERG <- load_hERG_fit_parameters("../../../hERG_fitting/results/astemizole_boot_pars.csv")
            new_hERG_pars <- astemizolehERG[7,]
            
            parsNew <- pars
            parsNew["Kmax"] <- new_hERG_pars$Kmax
            parsNew["Ku"] <- new_hERG_pars$Ku
            parsNew["n"] <- new_hERG_pars$n
            parsNew["halfmax"] <- new_hERG_pars$halfmax
            parsNew["Vhalf"] <- new_hERG_pars$Vhalf
            
            drugCmaxTable <- read.csv("../../Input/drug_cmaxes/cipa_drug_information.csv",
                                      stringsAsFactors = FALSE)
            cmaxes <- drugCmaxTable$therapeutic
            names(cmaxes) <- drugCmaxTable$drug
            Cmax.astemizole <- cmaxes["astemizole"]
            initstatesNew <- initstates
            initstatesNew["D"] <- Cmax.astemizole
            
            myOutput <- run_simulation(modelname,
                                       initstatesNew,
                                       parsNew,
                                       timepoints,
                                       events = NULL)
            expect_false(is.null(myOutput))
            expect_true(markov_states_sum_to_one(as.data.frame(myOutput),1e-6))
            expect_true(markov_states_greater_than_zero(as.data.frame(myOutput),1e-6))
            myOutputFrame <- as.data.frame(myOutput)
          })

test_that("Simulation with the ramp protocol.",
          {
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            timepoints <- 0:1000
            astemizolehERG <- load_hERG_fit_parameters("../../../hERG_fitting/results/astemizole_boot_pars.csv")
            new_hERG_pars <- astemizolehERG[7,]
            
            parsNew <- pars
            parsNew["Kmax"] <- new_hERG_pars$Kmax
            parsNew["Ku"] <- new_hERG_pars$Ku
            parsNew["n"] <- new_hERG_pars$n
            parsNew["halfmax"] <- new_hERG_pars$halfmax
            parsNew["Vhalf"] <- new_hERG_pars$Vhalf
            
            drugCmaxTable <- read.csv("../../Input/drug_cmaxes/cipa_drug_information.csv",
                                      stringsAsFactors = FALSE)
            cmaxes <- drugCmaxTable$therapeutic
            names(cmaxes) <- drugCmaxTable$drug
            Cmax.astemizole <- cmaxes["astemizole"]
            initstatesNew <- initstates
            initstatesNew["D"] <- Cmax.astemizole
            
            inputList <- list()
            inputList$tfinal <- 3000
            inputList$tstart <- 300
            inputList$tholding <- 500
            inputList$trelax <- 500
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            myVoltage <- voltage_step_program_list_mode(inputList)
            
            myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                          voltageProgram = myVoltage$voltage_program)
            
            initstatesVoltageRamp <- initstatesNew
            initstatesVoltageRamp["V"] <- myVoltage$voltage_program[1]
            myEvents <- list(data = myEvents)
            
            myOutput <- run_simulation(modelname = "hergmod",
                                       initstates = initstatesVoltageRamp,
                                       pars = parsNew,
                                       timepoints = myVoltage$timepoints,
                                       events = myEvents)
            
            myOutputDataFrame <- as.data.frame(myOutput)
            expect_false(is.null(myOutput))
            expect_true(markov_states_sum_to_one(myOutputDataFrame,1e-6))
            expect_true(markov_states_greater_than_zero(myOutputDataFrame,1e-6))
          })

test_that("Expect error if state names do not exactly match what is expected.",
          {
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            timepoints <- 0:1000
            astemizolehERG <- load_hERG_fit_parameters("../../../hERG_fitting/results/astemizole_boot_pars.csv")
            new_hERG_pars <- astemizolehERG[7,]
            
            parsNew <- pars
            parsNew["Kmax"] <- new_hERG_pars$Kmax
            parsNew["Ku"] <- new_hERG_pars$Ku
            parsNew["n"] <- new_hERG_pars$n
            parsNew["halfmax"] <- new_hERG_pars$halfmax
            parsNew["Vhalf"] <- new_hERG_pars$Vhalf
            
            drugCmaxTable <- read.csv("../../Input/drug_cmaxes/cipa_drug_information.csv",
                                      stringsAsFactors = FALSE)
            cmaxes <- drugCmaxTable$therapeutic
            names(cmaxes) <- drugCmaxTable$drug
            Cmax.astemizole <- cmaxes["astemizole"]
            initstates["D"] <- Cmax.astemizole
            inputList <- list()
            inputList$tfinal <- 3000
            inputList$tstart <- 300
            inputList$tholding <- 500
            inputList$trelax <- 500
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            myVoltage <- voltage_step_program_list_mode(inputList)
            
            myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                          voltageProgram = myVoltage$voltage_program)

            myEvents <- list(data = myEvents)
            
            #One bad state name, but all others correct.
            bad_initstates_1 <- data.frame(time_bad = 0.0,
                                         IC1 = 0.0,
                                         IC2 = 0.0,
                                         C1 = 1.0,
                                         C2 = 0.0,
                                         O = 0.0,
                                         IO = 0.0,
                                         IObound = 0.0,
                                         Obound = 0.0,
                                         Cbound = 0.0,
                                         D = 0.5,
                                         V = -40.0)
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = bad_initstates_1,
                             parsNew,
                             timepoints,
                             events = myEvents))
            
            #One extra state name, but all others correct.
            bad_initstates_2 <- data.frame(time = 0.0,
                                         IC1 = 0.0,
                                         IC2 = 0.0,
                                         C1 = 1.0,
                                         C2 = 0.0,
                                         O = 0.0,
                                         IO = 0.0,
                                         IObound = 0.0,
                                         Obound = 0.0,
                                         Cbound = 0.0,
                                         D = 0.5,
                                         V = -40.0,
                                         extra_nonsense = "stuff")
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = bad_initstates_2,
                             parsNew,
                             timepoints,
                             events = myEvents))
            
            #Completely bogus state names
            bad_initstates_3 <- data.frame(xyz1 = 1,
                                           xyz2 = 3,
                                           xyz3 = "a",
                                           xyz4 = 2 + 3i)
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = bad_initstates_3,
                             parsNew,
                             timepoints,
                             events = myEvents))
            
            #One missing state name, but all others are correct.
            bad_initstates_4 <- data.frame(time = 0.0,
                                           IC1 = 0.0,
                                           IC2 = 0.0,
                                           #C1 = 1, C1 will be missing.
                                           C2 = 0.0,
                                           O = 0.0,
                                           IO = 0.0,
                                           IObound = 0.0,
                                           Obound = 0.0,
                                           Cbound = 0.0,
                                           D = 0.5,
                                           V = -40.0)
            
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = bad_initstates_4,
                             parsNew,
                             timepoints,
                             events = myEvents))
            
          })

test_that("Expect error if pars is not numeric.",
          {
            modelname = "hergmod"
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            timepoints <- 0:1000
            astemizolehERG <- load_hERG_fit_parameters("../../../hERG_fitting/results/astemizole_boot_pars.csv")
            new_hERG_pars <- astemizolehERG[7,]
            
            parsNew <- pars
            parsNew["Kmax"] <- new_hERG_pars$Kmax
            parsNew["Ku"] <- new_hERG_pars$Ku
            parsNew["n"] <- new_hERG_pars$n
            parsNew["halfmax"] <- new_hERG_pars$halfmax
            parsNew["Vhalf"] <- new_hERG_pars$Vhalf
            
            drugCmaxTable <- read.csv("../../Input/drug_cmaxes/cipa_drug_information.csv",
                                      stringsAsFactors = FALSE)
            cmaxes <- drugCmaxTable$therapeutic
            names(cmaxes) <- drugCmaxTable$drug
            Cmax.astemizole <- cmaxes["astemizole"]
            initstates["D"] <- Cmax.astemizole
            inputList <- list()
            inputList$tfinal <- 3000
            inputList$tstart <- 300
            inputList$tholding <- 500
            inputList$trelax <- 500
            inputList$deltaV <- 120.0
            inputList$voltageOffset <- -80.0
            myVoltage <- voltage_step_program_list_mode(inputList)
            
            myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                          voltageProgram = myVoltage$voltage_program)
            
            myEvents <- list(data = myEvents)
            
            oldPars <- pars
            
            #One par field is not numeric.
            pars_bad_1 <- oldPars
            pars_bad_1["A1"] <- "stuff"
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = pars_bad_1,
                             timepoints,
                             events = myEvents))
            
            #All fields are not numeric.
            pars_bad_2 <- oldPars
            pars_bad_2[1:51] <- as.character(1:51)
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = pars_bad_2,
                             timepoints,
                             events = myEvents))
          })

test_that("Expect error if pars fields do not match the expected field names.",
          {
            modelname = "hergmod"
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            timepoints <- 0:1000
            astemizolehERG <- load_hERG_fit_parameters("../../../hERG_fitting/results/astemizole_boot_pars.csv")
            new_hERG_pars <- astemizolehERG[7,]
            
            parsNew <- pars
            parsNew["Kmax"] <- new_hERG_pars$Kmax
            parsNew["Ku"] <- new_hERG_pars$Ku
            parsNew["n"] <- new_hERG_pars$n
            parsNew["halfmax"] <- new_hERG_pars$halfmax
            parsNew["Vhalf"] <- new_hERG_pars$Vhalf
            
            drugCmaxTable <- read.csv("../../Input/drug_cmaxes/cipa_drug_information.csv",
                                      stringsAsFactors = FALSE)
            cmaxes <- drugCmaxTable$therapeutic
            names(cmaxes) <- drugCmaxTable$drug
            Cmax.astemizole <- cmaxes["astemizole"]
            initstates["D"] <- Cmax.astemizole
            inputList <- list()
            inputList$tfinal <- 3000
            inputList$tstart <- 300
            inputList$tholding <- 500
            inputList$trelax <- 500
            inputList$deltaV <- 120.0
            inputList$voltageOffset <- -80.0
            myVoltage <- voltage_step_program_list_mode(inputList)
            
            myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                          voltageProgram = myVoltage$voltage_program)
            
            myEvents <- list(data = myEvents)
            
            oldPars <- pars
            
            #One fieldname is incorrect, but all others are correct.
            parnames <- names(oldPars)
            bad_par_names_1 <- parnames
            bad_par_names_1[1] <- "bad_field_name"
            bad_pars_1 <- oldPars
            names(bad_pars_1) <- bad_par_names_1
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = bad_pars_1,
                             timepoints,
                             events = myEvents))
            
            #All fieldnames are incorrect except one.
            parnames <- names(oldPars)
            bad_par_names_2 <- parnames
            bad_par_names_2[2:51] <- paste0(as.character(2:51),rep("xyz",51-2+1))
            bad_pars_2 <- oldPars
            names(bad_pars_2) <- bad_par_names_2
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = bad_pars_2,
                             timepoints,
                             events = myEvents))
            
            #All fieldnames are correct, but there is one extra field that are unexpected.
            bad_pars_3 <- oldPars
            bad_pars_3["tomsawyer"] <- 4
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = bad_pars_3,
                             timepoints,
                             events = myEvents))
            
            #All fieldnames are correct, but there is MULTIPLE extra fields that are unexpected of various data types.
            bad_pars_4 <- oldPars
            bad_pars_4["tomsawyer"] <- 4.31
            bad_pars_4["huckleberryfinn"] <- 1.0
            bad_pars_4["papfinn"] <- -3.14159
            bad_pars_4["jim"] <- -10.0
            bad_pars_4["beckythatcher"] <- 1111.5
            bad_pars_4["marktwain"] <- 0.0
            bad_pars_4["auntpolly"] <- -1.0
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = bad_pars_4,
                             timepoints,
                             events = myEvents))
            
            #None of the fieldnames are correct, and there are fewer than the correct number.
            bad_pars_5 <- c(romeo = 1.0,
                            juliet = 2.0,
                            tybalt = 1e-9,
                            paris = 300.0,
                            benvolio = 3.14159,
                            mercutio = -77.1,
                            shakespeare = 0.0,
                            ladymontague = -1e-9,
                            ladycapulet = -1.1,
                            lordmontague = 1+1e-9,
                            lordcapulet = 1-1e-9)
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = bad_pars_5,
                             timepoints,
                             events = myEvents))
            
            #None of the fieldnames are correct, but their quantity is exactly equal to the correct number.
            bad_pars_6 <- c(hamlet = 0.0,
                            claudius = -10.0,
                            gertrude = 15.0,
                            polonius = -1000.0,
                            horatio = 1000.0,
                            ophelia = -1e-9,
                            laertes = 123.5,
                            fortinbras = 1.0,
                            rosencrantz = -100.5,
                            guildenstern = -33.0,
                            osric = 1e-9,
                            voltimand = 0.5 - 1e-9,
                            cornelius = 1 - 1e-9,
                            marcellus = -1 - 1e-9,
                            bernardo = -1 + 1e-9,
                            francisco = 0.5,
                            reynaldo = -0.5)
            
            bad_pars_6_append <- 1:(51 - 17)
            names(bad_pars_6_append) <- paste0(rep("filler",51-17),1:(51-17))
            bad_pars_6 <- c(bad_pars_6,bad_pars_6_append)
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = bad_pars_6,
                             timepoints,
                             events = myEvents))
            
            #None of the fieldnames are correct, and there is extra incorrect fieldnames.
            bad_pars_7 <- c(countdracula = 0.0,
                            minaharker = 1.0,
                            vanhelsing = -1.0,
                            jonathanharker = -1e-9,
                            lucywestenra = 1e-9,
                            renfield = 1 - 1e-9,
                            franksteinsmonster = 100.0,
                            johnseward = -123.5,
                            quinceymorris = 1 + 1e-9,
                            victorfrankenstein = 0.5,
                            bramstoker = 1000.0)
            bad_pars_7_append <- 1:(100 - 11)
            names(bad_pars_7_append) <- paste0(rep("filler",100 - 11),1:(100 - 11))
            bad_pars_7 <- c(bad_pars_7,bad_pars_7_append)
            
            expect_error(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = bad_pars_7,
                             timepoints = myVoltage$timepoints,
                             events = myEvents))
          })

test_that("Error if the modelname is not a character string.",
          {
            modelname <- 123.5
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            timepoints <- 0:1000
            
            expect_error(run_simulation(modelname,
                                        initstates,
                                        pars,
                                        timepoints,
                                        events = NULL))
          })

test_that("Expect that we can run the simulation without error with the 5-second ramp protocol.",
          {
            modelname = "hergmod"
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            astemizolehERG <- load_hERG_fit_parameters("../../../hERG_fitting/results/astemizole_boot_pars.csv")
            new_hERG_pars <- astemizolehERG[7,]
            
            parsNew <- pars
            parsNew["Kmax"] <- new_hERG_pars$Kmax
            parsNew["Ku"] <- new_hERG_pars$Ku
            parsNew["n"] <- new_hERG_pars$n
            parsNew["halfmax"] <- new_hERG_pars$halfmax
            parsNew["Vhalf"] <- new_hERG_pars$Vhalf
            pars <- parsNew
            drugCmaxTable <- read.csv("../../Input/drug_cmaxes/cipa_drug_information.csv",
                                      stringsAsFactors = FALSE)
            cmaxes <- drugCmaxTable$therapeutic
            names(cmaxes) <- drugCmaxTable$drug
            Cmax.astemizole <- cmaxes["astemizole"]
            initstates["D"] <- Cmax.astemizole
            inputList <- list()
            inputList <- read.csv("../../Input/voltage_program_data_files/5_second_ramp_protocol.csv",stringsAsFactors = FALSE)
            expect_silent(myVoltage <- voltage_step_program_list_mode(inputList))
            
            expect_silent(myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                          voltageProgram = myVoltage$voltage_program))

            expect_silent(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = pars,
                             timepoints = myVoltage$timepoints,
                             events = myEvents))
            })

test_that("Expect that we can run the simulation without error with the 30-second ramp protocol.",
          {
            modelname = "hergmod"
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            astemizolehERG <- load_hERG_fit_parameters("../../../hERG_fitting/results/astemizole_boot_pars.csv")
            new_hERG_pars <- astemizolehERG[7,]
            
            parsNew <- pars
            parsNew["Kmax"] <- new_hERG_pars$Kmax
            parsNew["Ku"] <- new_hERG_pars$Ku
            parsNew["n"] <- new_hERG_pars$n
            parsNew["halfmax"] <- new_hERG_pars$halfmax
            parsNew["Vhalf"] <- new_hERG_pars$Vhalf
            pars <- parsNew
            drugCmaxTable <- read.csv("../../Input/drug_cmaxes/cipa_drug_information.csv",
                                      stringsAsFactors = FALSE)
            cmaxes <- drugCmaxTable$therapeutic
            names(cmaxes) <- drugCmaxTable$drug
            Cmax.astemizole <- cmaxes["astemizole"]
            initstates["D"] <- Cmax.astemizole
            inputList <- list()
            inputList <- read.csv("../../Input/voltage_program_data_files/30_second_ramp_protocol.csv",stringsAsFactors = FALSE)
            expect_silent(myVoltage <- voltage_step_program_list_mode(inputList))
            
            expect_silent(myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                                        voltageProgram = myVoltage$voltage_program))
            
            expect_silent(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = pars,
                             timepoints = myVoltage$timepoints,
                             events = myEvents))
          })

test_that("Expect that we can run the simulation without error with the experimental action potential.",
          {
            modelname = "hergmod"
            initstates <- read.table("../../Input/hERG_parms/hergmod_states.txt")
            stateNames <- initstates$V1
            initstates <- initstates$V2
            names(initstates) <- stateNames
            pars <- read.table("../../Input/hERG_parms/hergmod_pars.txt")
            parNames <- pars$V1
            pars <- pars$V2
            names(pars) <- parNames
            astemizolehERG <- load_hERG_fit_parameters("../../../hERG_fitting/results/astemizole_boot_pars.csv")
            new_hERG_pars <- astemizolehERG[7,]
            
            parsNew <- pars
            parsNew["Kmax"] <- new_hERG_pars$Kmax
            parsNew["Ku"] <- new_hERG_pars$Ku
            parsNew["n"] <- new_hERG_pars$n
            parsNew["halfmax"] <- new_hERG_pars$halfmax
            parsNew["Vhalf"] <- new_hERG_pars$Vhalf
            pars <- parsNew
            drugCmaxTable <- read.csv("../../Input/drug_cmaxes/cipa_drug_information.csv",
                                      stringsAsFactors = FALSE)
            cmaxes <- drugCmaxTable$therapeutic
            names(cmaxes) <- drugCmaxTable$drug
            Cmax.astemizole <- cmaxes["astemizole"]
            initstates["D"] <- Cmax.astemizole
            
            inputList <- list()
            inputList <- read.table("../../Input/voltage_program_data_files/AP2Hz",header = TRUE)

            expect_silent(myVoltage <- voltage_trueAP_program_list_mode(inputList))
            expect_silent(myEvents <- generate_events_voltage_dataframe(t.time = myVoltage$timepoints,
                                                                        voltageProgram = myVoltage$voltage_program))
            
            expect_silent(
              run_simulation(modelname = "hergmod",
                             initstates = initstates,
                             pars = pars,
                             timepoints = myVoltage$timepoints,
                             events = myEvents))
          })