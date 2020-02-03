context("testing is_between.R")

test_that("Test that all errors trip as they should.",
          {
            expect_error(is_between("a",1,3))
            expect_error(is_between(2,"a",3))
            expect_error(is_between(2,1,"a"))
            q <- c(1,2,3)
            expect_error(is_between(q,1,3))
            expect_error(is_between(2,q,3))
            expect_error(is_between(2,1,q))
  
            expect_error(is_between(testValue = 2,
                                   leftBoundary = 5,
                                   rightBoundary = 1))
          })

test_that("Test that the function works on some corner cases.",
          {
            expect_true(is_between(1,1,3))
            expect_true(is_between(2,2,3))
            expect_true(is_between(1,1,1))
            
          })

test_that("Test that function gives a correct result for some example cases.",
          {
            expect_true(is_between(2,2.5,3))
            expect_true(is_between(0.1,1.0,5))
            expect_true(is_between(-2,-1,0))
            expect_true(is_between(-2000,-1000,500))
            
            expect_false(is_between(2,1,5))
            expect_false(is_between(-1,-2,1))
            expect_false(is_between(0.01,0.001,0.05))
          })