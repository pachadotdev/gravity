## ----setup, cache = FALSE, echo = FALSE, message = FALSE, warning = FALSE, tidy = FALSE----
knitr::opts_chunk$set(eval = TRUE, message = FALSE, warning = FALSE)

## ----ddm----------------------------------------------------------------------
library(gravity)

fit <- ddm(
    dependent_variable = "flow",
    distance = "distw",
    additional_regressors = c("rta", "comcur", "contig"),
    code_origin = "iso_o",
    code_destination = "iso_d",
    data = gravity_no_zeros
  )

summary(fit)

## ----bvu----------------------------------------------------------------------
fit2 <- bvu(
    dependent_variable = "flow",
    distance = "distw",
    additional_regressors = c("rta", "contig", "comcur"),
    income_origin = "gdp_o",
    income_destination = "gdp_d",
    code_origin = "iso_o",
    code_destination = "iso_d",
    data = gravity_no_zeros
  )

summary(fit2)

## ----bvw----------------------------------------------------------------------
fit3 <- bvw(
    dependent_variable = "flow",
    distance = "distw",
    additional_regressors = c("rta", "comcur", "contig"),
    income_origin = "gdp_o",
    income_destination = "gdp_d",
    code_origin = "iso_o",
    code_destination = "iso_d",
    data = gravity_no_zeros
  )

summary(fit3)

## ----ppml---------------------------------------------------------------------
fit4 <- ppml(
    dependent_variable = "flow",
    distance = "distw",
    additional_regressors = c("rta", "comcur", "contig"),
    data = gravity_no_zeros
  )

summary(fit4)

## ----ppmlr--------------------------------------------------------------------
fit4r <- ppml(
    dependent_variable = "flow",
    distance = "distw",
    additional_regressors = c("rta", "comcur", "contig"),
    robust = TRUE,
    data = gravity_no_zeros
  )

summary(fit4r)

## ----gpml---------------------------------------------------------------------
fit5 <- gpml(
  dependent_variable = "flow",
  distance = "distw",
  additional_regressors = c("rta", "comcur", "contig"),
  robust = TRUE,
  data = gravity_no_zeros
)

summary(fit5)

## ----nbpml--------------------------------------------------------------------
fit6 <- nbpml(
  dependent_variable = "flow",
  distance = "distw",
  additional_regressors = c("rta", "comcur", "contig"),
  robust = TRUE,
  data = gravity_no_zeros
)

summary(fit6)

## ----tetrads------------------------------------------------------------------
fit8 <- tetrads(
  dependent_variable = "flow",
  distance = "distw",
  additional_regressors = "rta",
  code_origin = "iso_o",
  code_destination = "iso_d",
  filter_origin = "CHN",
  filter_destination = "USA",
  data = gravity_no_zeros
)

summary(fit8)

## ----tetrads2-----------------------------------------------------------------
fit8 <- tetrads(
  dependent_variable = "flow",
  distance = "distw",
  additional_regressors = "rta",
  code_origin = "iso_o",
  code_destination = "iso_d",
  filter_origin = "CHN",
  filter_destination = "USA",
  multiway = TRUE,
  data = gravity_no_zeros
)

summary(fit8)

