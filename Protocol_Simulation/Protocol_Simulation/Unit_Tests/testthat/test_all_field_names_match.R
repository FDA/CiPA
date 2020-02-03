context("Testing all_field_names_match.R")

requiredNames <- c("a","b","c")
actualNames <- c("q","r","s")
test_that("Expect FALSE when names do not match at all.",
          {
            expect_false(all_field_names_match(requiredNames,actualNames))
          })

t1 <- c("a","b","c")
t2 <- c("a","b")
test_that("Expect FALSE when a field is missing from either character vector.",
          {
            expect_false(all_field_names_match(t1,t2))
            expect_false(all_field_names_match(t2,t1))
          })

t1 <- ""
t2 <- ""
test_that("Expect TRUE if both fields are blank.",
          {
            expect_true(all_field_names_match(t1,t2))
          })

t1 <- c("a","b","c")
t2 <- ""
test_that("Expect FALSE if one field is blank",
          {
            expect_false(all_field_names_match(t1,t2))
            expect_false(all_field_names_match(t2,t1))
          })

t1 <- c("a","b","c")
t2 <- c("a","b","c")
test_that("Expect TRUE if the fields are all the same.",
          {
            expect_true(all_field_names_match(t1,t2))
          })

t1 <- c("d","c","a","b")
t2 <- c("a","b","c","d")
test_that("Expect TRUE is the fields are all the same, but disordered.",
          {
            expect_true(all_field_names_match(t1,t2))
          })

test_that("Expect error if neither input is a character vector.",
          {
            expect_error(all_field_names_match(0,c("a","b","c")))
            expect_error(all_field_names_match(c("a","b","c"),0))
            expect_error(all_field_names_match(23L,-1))
            expect_error(all_field_names_match(1.5+2.7i,c("")))
          })