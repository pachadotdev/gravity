context("test-stata_summary.R")

test_that("Stata summary returns a valid output", {
  # Example for CRAN checks:
  # Executable in < 5 sec
  library(dplyr)
  data("gravity_no_zeros")

  # Choose 5 countries for testing
  countries_chosen <- c("AUS", "CHN", "GBR", "BRA", "CAN")
  grav_small <- filter(gravity_no_zeros, iso_o %in% countries_chosen)

  # Using OLS for testing
  fit <- ols(
    dependent_variable = "flow",
    distance = "distw",
    additional_regressors = c("rta", "contig", "comcur"),
    income_origin = "gdp_o",
    income_destination = "gdp_d",
    code_origin = "iso_o",
    code_destination = "iso_d",
    uie = FALSE,
    robust = FALSE,
    data = grav_small
  )

  fit2 <- hm_summary(fit, robust = FALSE)

  expect_is(fit2, "summary.lm")
})
