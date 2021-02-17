context("test-gravity.R")

test_that("DDM returns a valid output", {
  # Example for CRAN checks:
  # Executable in < 5 sec
  library(dplyr)
  data("gravity_no_zeros")

  # Choose 5 countries for testing
  countries_chosen <- c("AUS", "CHN", "GBR", "BRA", "CAN")
  grav_small <- filter(gravity_no_zeros, iso_o %in% countries_chosen)

  fit <- ddm(
    dependent_variable = "flow",
    distance = "distw",
    additional_regressors = c("rta", "comcur", "contig"),
    code_origin = "iso_o",
    code_destination = "iso_d",
    robust = FALSE,
    data = grav_small
  )

  expect_is(fit, "lm")

  fit2 <- ddm(
    dependent_variable = "flow",
    distance = "distw",
    code_origin = "iso_o",
    code_destination = "iso_d",
    robust = FALSE,
    data = grav_small
  )

  expect_is(fit2, "lm")
})
