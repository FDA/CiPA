context("Testing hERG fitting parameters loading function.")

test_that("Expect a 2000-row table of numeric values from the good astemizole boot_pars.csv dataset.",
          {
            myTable <- load_hERG_fit_parameters("../../../hERG_fitting/results/astemizole_boot_pars.csv")
            expect_true(is.data.frame(myTable))
            expect_true(is.numeric(myTable$Kmax))
            expect_true(is.numeric(myTable$Ku))
            expect_true(is.numeric(myTable$n))
            expect_true(is.numeric(myTable$halfmax))
            expect_true(is.numeric(myTable$Vhalf))
            expect_true(is.numeric(myTable$slope))
            expect_true(nrow(myTable) == 2000L)
          })

test_that("Expect error from a 10-row table of garbage values from the crap astemizole boot_pars.csv dataset.",
          {
            myTable <- expect_error(load_hERG_fit_parameters("../../Input/bad_boot_pars_for_testing/astemizole_boot_pars_crap.csv"))
          })

test_that("Expect errors if not a string input.",
          {
            expect_error(load_hERG_fit_parameters(1))
            expect_error(load_hERG_fit_parameters(-1))
            expect_error(load_hERG_fit_parameters(0))
            expect_error(load_hERG_fit_parameters(23L))
            expect_error(load_hERG_fit_parameters(1+2.333i))
          })

test_that("Expect errors if not a .csv file.",
          {
            expect_error(load_hERG_fit_parameters("astemizole.cvs"))
            expect_error(load_hERG_fit_parameters("/home/stuff/0.c"))
            expect_error(load_hERG_fit_parameters("1.vcs"))
            expect_error(load_hERG_fit_parameters("csv.svc"))
            expect_error(load_hERG_fit_parameters("/home/mystuff/table-1.vsc"))
          })