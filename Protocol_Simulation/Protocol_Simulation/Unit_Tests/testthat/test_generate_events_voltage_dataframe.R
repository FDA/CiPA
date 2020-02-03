context("Testing generate_events_voltage_dataframe.R.")

t.time <- 0:10
voltageProgram <- t.time^2
test_that("Test that assertions trigger with bad inputs.",
          {
            expect_error(generate_events_voltage_dataframe("badstring",voltageProgram))
            expect_error(generate_events_voltage_dataframe(-t.time,voltageProgram))
            expect_error(generate_events_voltage_dataframe(2+3i,voltageProgram))
            
            expect_error(generate_events_voltage_dataframe(t.time,"badstring"))
            expect_error(generate_events_voltage_dataframe(t.time,2+3i))
          })

test_that("Expect error if t.time and voltageProgram are of different lengths.",
          {
            expect_error(generate_events_voltage_dataframe(c(1,2,3),c(30,30)))
            expect_error(generate_events_voltage_dataframe(c(1,2),c(30,30,30)))
          })

requiredFieldNames <- c("var","time","value","method")
myFrame <- generate_events_voltage_dataframe(t.time,voltageProgram)
test_that("Expect fieldnames to match the required field names.",
          {
            expect_true(all_field_names_match(requiredFieldNames,names(myFrame)))
          })

test_that("Expect number of rows to be equal to length of input vectors.",
          {
            expect_true(nrow(myFrame) == length(t.time))
            expect_true(nrow(myFrame) == length(voltageProgram))
          })

test_that("Expect an input voltage program to produce a known a priori events table.",
          {
            inputVoltageProgram <- c(1,1,1,2,2,3,4,5,6,6,6,-7)
            t.time <- 0:(length(inputVoltageProgram)-1)
            predictedFrame <- data.frame(var = rep("V",7),
                                         time = c(0,3,5,6,7,8,11),
                                         value = c(1,2,3,4,5,6,-7),
                                         method = rep("replace",7))
            
            testedFrame <- generate_events_voltage_dataframe(t.time,inputVoltageProgram)
            expect_true(all.equal(predictedFrame,testedFrame))
          })