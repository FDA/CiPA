context("Testing is_csv.R")

test_that("Expect FALSE when file does not have a .csv file extension.",
          {
            expect_false(is_csv("c.sv"))
            expect_false(is_csv("1.c"))
            expect_false(is_csv(".cs"))
            expect_false(is_csv(""))
            expect_false(is_csv("backwards.vsc"))
            expect_false(is_csv("csv.svc"))
            expect_false(is_csv("/home/directory/onemore/csv.vcs"))
          })

test_that("Expect TRUE when file does have .csv file extension.",
          {
            expect_true(is_csv(".csv"))
            expect_true(is_csv("1.csv"))
            expect_true(is_csv("csv.csv"))
            expect_true(is_csv("vcs.csv"))
            expect_true(is_csv("a.csv"))
            expect_true(is_csv("123.csv"))
            expect_true(is_csv("321abc.csv"))
            expect_true(is_csv("0.csv"))
            expect_true(is_csv("/home/somewhere/playerunknown/battlegrounds/0.csv"))
          })

test_that("Expect an error if not a string input.",
          {
            expect_error(is_csv(1))
            expect_error(is_csv(0))
            expect_error(is_csv(-1))
            expect_error(is_csv(23L))
            expect_error(is_csv(1+2i))
            expect_error(is_csv(NA))
            expect_error(is_csv(NULL))
          })