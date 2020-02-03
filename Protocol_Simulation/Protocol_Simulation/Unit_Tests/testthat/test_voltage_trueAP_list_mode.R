context("testing voltage_trueAP_list_mode.R")
test_that("Test that the output is of length 101 for a known input.",
          {
            inputList <- list()
            inputList$Time <- 0:10
            inputList$V <- 0.1 * inputList$Time^2
            expect_silent(trueAP_voltage_data <- voltage_trueAP_program_list_mode(inputList))
            expect_equal(11, length(trueAP_voltage_data$timepoints))
            expect_equal(11, length(trueAP_voltage_data$voltage_program))
          })

test_that("expect errors if there are missing or extra inputList fields.",
          {
            #Missing one required field.
            inputList <- list()
            #inputList$Time <- 0:2000
            inputList$V <- 1:100
            expect_error(voltage_step_program_list_mode(inputList))
            
            #One extra field.
            inputList <- list()
            inputList$Time <- 0:20
            inputList$V <- 0:20
            inputList$DonaldTrump <- "The code just got ten lines longer!"
            expect_error(voltage_step_program_list_mode(inputList))
          })

test_that("expect an error if the inputList is not a list.",
          {
            expect_error(voltage_trueAP_program_list_mode("stuff"))
            expect_error(voltage_trueAP_program_list_mode(1.0))
            expect_error(voltage_trueAP_program_list_mode(2L))
            expect_error(voltage_trueAP_program_list_mode(Inf))
            expect_error(voltage_trueAP_program_list_mode(-Inf))
            expect_error(voltage_trueAP_program_list_mode(373.19))
            expect_error(voltage_trueAP_program_list_mode(-1000.5))
            expect_error(voltage_trueAP_program_list_mode(1e-9))
            expect_error(voltage_trueAP_program_list_mode(-1e-9))
            expect_error(voltage_trueAP_program_list_mode(1 - 1e-9))
            expect_error(voltage_trueAP_program_list_mode(1 + 1e-9))
            expect_error(voltage_trueAP_program_list_mode(-1 - 1e-9))
            expect_error(voltage_trueAP_program_list_mode(-1 + 1e-9))
            expect_error(voltage_trueAP_program_list_mode(0.0))
            expect_error(voltage_trueAP_program_list_mode(0L))
            expect_error(voltage_trueAP_program_list_mode(1L))
            expect_error(voltage_trueAP_program_list_mode(-1L))
            expect_error(voltage_trueAP_program_list_mode(1:100))
          })

test_that("Expect error if the time vector is not strictly monotonically increasing.",
          {
            #Entire vector is monotonically decreasing.
            inputList <- list()
            inputList$Time <- 100:0
            inputList$V <- 0.1 * inputList$Time^2
            expect_error(voltage_trueAP_program_list_mode(inputList))
            
            #A single pair of adjacent time points in a monotonically increasing sequence are flipped.
            inputList <- list()
            inputList$Time <- 0:100
            inputList$Time[2] <- 0
            inputList$Time[1] <- 1
            inputList$V <- 0.1 * inputList$Time^2
            expect_error(voltage_trueAP_program_list_mode(inputList))
          })

test_that("Expect that we can get the voltage program for the experimental 2-second action potential data file.",
          {
            inputList <- list()
            inputList <- read.table("../../Input/voltage_program_data_files/AP2Hz",header = TRUE)
            
            expect_silent(mystuff <- voltage_trueAP_program_list_mode(inputList))
            expect_equal(length(mystuff$voltage_program),1002)
            expect_equal(max(mystuff$voltage_program),44.8)
            expect_equal(min(mystuff$voltage_program),-94.3)
          })

test_that("Expect errors if the input file is erroneous.",
          {
            #Error if the time vector is not monotonically increasing.
            inputList <- list()
            inputList <- read.table("../../Input/voltage_program_data_files/AP2Hz_bad_data_for_unit_testing_1",
                                    header = TRUE)
            
            expect_error(voltage_trueAP_program_list_mode(inputList))
            
            #Error if the time vector is monotonically decreasing.
            inputList <- list()
            inputList <- read.table("../../Input/voltage_program_data_files/AP2Hz_bad_data_for_unit_testing_2",
                                    header = TRUE)
            
            expect_error(voltage_trueAP_program_list_mode(inputList))
            
            #Error if any time value is negative.
            inputList <- list()
            inputList <- read.table("../../Input/voltage_program_data_files/AP2Hz_bad_data_for_unit_testing_3",
                                    header = TRUE)
            
            expect_error(voltage_trueAP_program_list_mode(inputList))
          })

test_that("Expect warnings if the time does not appear to be in milliseconds.",
          {
            #Expect warning if the time is in seconds and not milliseconds.
            inputList <- list()
            inputList <- read.table("../../Input/voltage_program_data_files/AP2Hz",header = TRUE)
            inputList$Time <- inputList$Time/1000
            expect_warning(voltage_trueAP_program_list_mode(inputList),
                           "Total length of input time is being read as under 50 milliseconds. Have the time data been converted to milliseconds?")
            
            #Expect warning if the time is extremely long (greater than 1 minute).
            inputList <- list()
            inputList <- read.table("../../Input/voltage_program_data_files/AP2Hz",header = TRUE)
            inputList$Time <- inputList$Time * 1000
            expect_warning(voltage_trueAP_program_list_mode(inputList),
                           "Total length of input time is beyond 60,000 milliseconds, or beyond 1 minute. This is probably too long. Has a mistake been made?")
            
            #Expect warning if the absolute value of the action potential voltage is very high and appears to not be expresed in millivolts.
            inputList <- list()
            inputList <- read.table("../../Input/voltage_program_data_files/AP2Hz",header = TRUE)
            inputList$Time <- inputList$Time
            inputList$V <- inputList$V * 1000
            expect_warning(voltage_trueAP_program_list_mode(inputList),
                           "Maximum absolute value of observed voltage in the action potential raw data exceeds 1000 millivolts, which does not appear biologically possible. Has a mistake been made?")
            
            #Expect warning if the absolute value of the action potential voltage is very low and appears to not be expresed in millivolts (probably in volts instead).
            inputList <- list()
            inputList <- read.table("../../Input/voltage_program_data_files/AP2Hz",header = TRUE)
            inputList$Time <- inputList$Time
            inputList$V <- inputList$V / 1000
            expect_warning(voltage_trueAP_program_list_mode(inputList),
                           "Maximum absolute value of observed voltage in the action potential raw data is under 1.0 millivolts, which is far below that expected of an action potential. Has a mistake been made?")
          })