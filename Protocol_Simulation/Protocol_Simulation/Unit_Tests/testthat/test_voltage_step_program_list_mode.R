context("testing voltage_step_program_list_mode.R.")
test_that("output has same length equal to the length of 0:tfinal.",
          {
            inputList <- list()
            inputList$tfinal <- 1700
            inputList$tstart <- 100
            inputList$tholding <- 100
            inputList$trelax <- 500
            inputList$deltaV <- 120.0
            inputList$voltageOffset <- -80.0
            expect_silent(myData <- voltage_step_program_list_mode(inputList))
            expect_equal(length(0:inputList$tfinal),length(myData$voltage_program))
          })

test_that("expect errors if there are missing or extra inputList fields.",
          {
            #Missing one required field.
            inputList <- list()
            #inputList$tfinal <- 1700
            inputList$tstart <- 100
            inputList$tholding <- 100
            inputList$trelax <- 500
            inputList$deltaV <- 120.0
            inputList$voltageOffset <- -80.0
            expect_error(voltage_step_program_list_mode(inputList))
            
            #One extra field.
            inputList <- list()
            inputList$tfinal <- 1700
            inputList$tstart <- 100
            inputList$tholding <- 100
            inputList$trelax <- 500
            inputList$deltaV <- 120.0
            inputList$voltageOffset <- -80.0
            inputList$DonaldTrump <- "The code just got ten lines longer!"
            expect_error(voltage_step_program_list_mode(inputList))
          })

test_that("expect an error if the inputList is not a list.",
          {
            expect_error(voltage_step_program_list_mode("stuff"))
            expect_error(voltage_step_program_list_mode(1.0))
            expect_error(voltage_step_program_list_mode(2L))
            expect_error(voltage_step_program_list_mode(Inf))
            expect_error(voltage_step_program_list_mode(-Inf))
            expect_error(voltage_step_program_list_mode(373.19))
            expect_error(voltage_step_program_list_mode(-1000.5))
            expect_error(voltage_step_program_list_mode(1e-9))
            expect_error(voltage_step_program_list_mode(-1e-9))
            expect_error(voltage_step_program_list_mode(1 - 1e-9))
            expect_error(voltage_step_program_list_mode(1 + 1e-9))
            expect_error(voltage_step_program_list_mode(-1 - 1e-9))
            expect_error(voltage_step_program_list_mode(-1 + 1e-9))
            expect_error(voltage_step_program_list_mode(0.0))
            expect_error(voltage_step_program_list_mode(0L))
            expect_error(voltage_step_program_list_mode(1L))
            expect_error(voltage_step_program_list_mode(-1L))
          })

test_that("expect errors if we try to put in negative numbers for times.",
          {
            inputList <- list()
            inputList$tfinal <- -1700
            inputList$tstart <- 100
            inputList$tholding <- 100
            inputList$trelax <- 500
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            expect_error(voltage_step_program_list_mode(inputList))
            
            inputList <- list()
            inputList$tfinal <- 1700
            inputList$tstart <- -100
            inputList$tholding <- 100
            inputList$trelax <- 500
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            expect_error(voltage_step_program_list_mode(inputList))
            
            inputList <- list()
            inputList$tfinal <- 1700
            inputList$tstart <- 100
            inputList$tholding <- -100
            inputList$trelax <- 500
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            expect_error(voltage_step_program_list_mode(inputList))
            
            inputList <- list()
            inputList$tfinal <- 1700
            inputList$tstart <- 100
            inputList$tholding <- 100
            inputList$trelax <- -500
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            expect_error(voltage_step_program_list_mode(inputList))

          })

test_that("expect error if tfinal < tstart + tholding + trelax.",
          {
            inputList <- list()
            inputList$tfinal <- 70
            inputList$tstart <- 100
            inputList$tholding <- 100
            inputList$trelax <- 500
            inputList$deltaV <- 120
            inputList$voltageOffset <- -80
            expect_error(voltage_step_program_list_mode(inputList))
          })

test_that("expect that we can read the 5-second csv and load it into the voltage program.",
          {
            inputList <- list()
            inputList <- read.csv("../../Input/voltage_program_data_files/5_second_ramp_protocol.csv")
            expect_silent(mystuff <- voltage_step_program_list_mode(inputList))
            expect_equal(length(mystuff$voltage_program),5001)
            expect_equal(max(mystuff$voltage_program),40)
            expect_equal(min(mystuff$voltage_program),-80)
          })

test_that("expect that we can read the 30-second csv and load it into the voltage program.",
          {
            inputList <- list()
            inputList <- read.csv("../../Input/voltage_program_data_files/30_second_ramp_protocol.csv")
            expect_silent(mystuff <- voltage_step_program_list_mode(inputList))
            expect_equal(length(mystuff$voltage_program),30001)
            expect_equal(max(mystuff$voltage_program),40)
            expect_equal(min(mystuff$voltage_program),-80)
          })

test_that("expect that we can read the ramp protocol for unit-testing in and get correct outputs.",
          {
            inputList <- list()
            inputList <- read.csv("../../Input/voltage_program_data_files/ramp_protocol_for_testthat.csv")
            expect_silent(mystuff <- voltage_step_program_list_mode(inputList))
            expect_equal(length(mystuff$voltage_program),17129 + 1)
            expect_equal(max(mystuff$voltage_program),196 - 72)
            expect_equal(min(mystuff$voltage_program),-72)
          })