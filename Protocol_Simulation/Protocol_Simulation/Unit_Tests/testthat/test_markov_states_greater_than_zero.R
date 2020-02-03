context("Testing markov_states_greater_than_zero.R")

blankFrame <- data.frame()
test_that("Expect errors if inputs are of incorrect type or incorrect domain.",
          {
            expect_error(markov_states_greater_than_zero("a",1e-6))
            expect_error(markov_states_greater_than_zero(0L,1e-6))
            expect_error(markov_states_greater_than_zero(0.0,1e-6))
            expect_error(markov_states_greater_than_zero(5.325,1e-6))
            expect_error(markov_states_greater_than_zero(-1L,1e-6))
            expect_error(markov_states_greater_than_zero(-1.234,1e-6))
            expect_error(markov_states_greater_than_zero(Inf,1e-6))
            expect_error(markov_states_greater_than_zero(-Inf,1e-6))
            expect_error(markov_states_greater_than_zero("",1e-6))
            
            expect_error(markov_states_greater_than_zero(blankFrame,"a"))
            expect_error(markov_states_greater_than_zero(blankFrame,"1e-6"))
            expect_error(markov_states_greater_than_zero(blankFrame,c(1,2,3)))
            
            expect_error(markov_states_greater_than_zero(blankFrame,-1.2))
          })

goodFrame <- data.frame(IC1 = 0.5,
                        IC2 = 0.25,
                        C1 = 0.25,
                        C2 = 0,
                        O = 0,
                        IO = 0,
                        IObound = 0,
                        Obound = 0,
                        Cbound = 0)

goodFrameLong <- data.frame(IC1 = c(0.5,-1e-9,-1e-9),
                            IC2 = c(0.25,1,-1e-9),
                            C1 = c(0.25,0,-1e-9),
                            C2 = c(0,0,-1e-9),
                            O = c(0,0,-1e-9),
                            IO = c(0,0,-1e-9),
                            IObound = c(0,0,-1e-9),
                            Obound = c(0,0,-1e-9),
                            Cbound = c(0,0,1))

test_that("Expect TRUE if all values are greater than 0 within negativeTol.",
          {
            expect_true(markov_states_greater_than_zero(goodFrameLong,1e-6))
          })

badFrameLong1 <- data.frame(IC1 = c(-0.5,-1e-9,-1e-9),
                            IC2 = c(0.25,1,-1e-9),
                            C1 = c(0.25,0,-1e-9),
                            C2 = c(0,0,-1e-9),
                            O = c(0,0,-1e-9),
                            IO = c(0,0,-1e-9),
                            IObound = c(0,0,-1e-9),
                            Obound = c(0,0,-1e-9),
                            Cbound = c(0,0,1))

badFrameLong2 <- data.frame(IC1 = c(0.5,-1e-9,-1e-6),
                            IC2 = c(0.25,1,-1e-6),
                            C1 = c(0.25,0,-1e-6),
                            C2 = c(0,0,-1e-6),
                            O = c(0,0,-1e-6),
                            IO = c(0,0,-1e-6),
                            IObound = c(0,0,-1e-6),
                            Obound = c(0,0,-1e-6),
                            Cbound = c(0,0,1))

test_that("Expect FALSE if any element is >= (0 - negativeTol).",
          {
            expect_false(markov_states_greater_than_zero(badFrameLong1,1e-6))
            expect_false(markov_states_greater_than_zero(badFrameLong2,1e-6))
          })