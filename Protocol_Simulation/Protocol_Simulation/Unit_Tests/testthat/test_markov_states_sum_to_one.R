context("Testing markov_states_sum_to_one.R")

test_that("Expect errors if stateDataFrame is not a dataframe.",
          {
            expect_error(markov_states_sum_to_one("a",1e-6))
            expect_error(markov_states_sum_to_one(23L,1e-6))
            expect_error(markov_states_sum_to_one(1.5423,1e-6))
            expect_error(markov_states_sum_to_one(-5.293,1e-6))
            expect_error(markov_states_sum_to_one(0,1e-6))
          })

stateDataFrame <- data.frame(IC1 = 0.5,
                             IC2 = 0.25,
                             C1 = 0.25,
                             C2 = 0,
                             O = 0,
                             IO = 0,
                             IObound = 0,
                             Obound = 0,
                             Cbound = 0)

test_that("Expect errors if sumTol is not a numeric single number.",
          {
            expect_error(markov_states_sum_to_one(stateDataFrame,"a"))
            expect_error(markov_states_sum_to_one(stateDataFrame,23+56i))
            expect_error(markov_states_sum_to_one(stateDataFrame,c(1,2)))
          })

stateDataFrameLong <- data.frame(IC1 = c(0.50,1,0,0,1 - 1e-9),
                                 IC2 = c(0.25,0,0.01,0,0),
                                 C1 = c(0.25,0,0.99,0,0),
                                 C2 = c(0,0,0,0.1,0),
                                 O = c(0,0,0,0.1,0),
                                 IO = c(0,0,0,0.1,0),
                                 IObound = c(0,0,0,0.1,0),
                                 Obound = c(0,0,0,0.1,0),
                                 Cbound = c(0,0,0,0.5,0)
                                 )

test_that("Expect TRUE if all rows sum to 1 within sumTol.",
          {
            expect_true(markov_states_sum_to_one(stateDataFrameLong,1e-6))
          })

stateDataFramePoorData1 <- data.frame(IC1 = c(0.50,1,0,0,1 - 1e-9),
                                 IC2 = c(0.25,0,0.01,0,0),
                                 C1 = c(0.25,0,0.99,0,0),
                                 C2 = c(0.1,0,0,0.1,0),
                                 O = c(0,0,0,0.1,0),
                                 IO = c(0,0,0,0.1,0),
                                 IObound = c(0,0,0,0.1,0),
                                 Obound = c(0,0,0,0.1,0),
                                 Cbound = c(0,0,0,0.5,0)
                                 )

stateDataFramePoorData2 <- data.frame(IC1 = c(0.50,1,0,0,1 - 1e-3),
                                      IC2 = c(0.25,0,0.01,0,0),
                                      C1 = c(0.25,0,0.99,0,0),
                                      C2 = c(0.0,0,0,0.1,0),
                                      O = c(0,0,0,0.1,0),
                                      IO = c(0,0,0,0.1,0),
                                      IObound = c(0,0,0,0.1,0),
                                      Obound = c(0,0,0,0.1,0),
                                      Cbound = c(0,0,0,0.5,0)
                                      )

stateDataFramePoorData3 <- data.frame(IC1 = 1,
                                      IC2 = 2,
                                      C1 = 3,
                                      C2 = 4,
                                      O = 5,
                                      IO = 6,
                                      IObound = 7,
                                      Obound = 8,
                                      Cbound = 9
                                      )

test_that("Expect FALSE if any row does not sum to 1 within sumTol.",
          {
            expect_false(markov_states_sum_to_one(stateDataFramePoorData1,1e-6))
            expect_false(markov_states_sum_to_one(stateDataFramePoorData2,1e-6))
            expect_false(markov_states_sum_to_one(stateDataFramePoorData3,1e-6))
          })