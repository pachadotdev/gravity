#' @title \insertCite{Eaton1995;textual}{gravity} threshold Tobit model (ET Tobit)
#'
#' @description \code{et_tobit} estimates gravity models in their additive form
#' by conducting a left-censored regression.
#'
#' @details \code{et_tobit} represents the \insertCite{Eaton1995;textual}{gravity} Tobit model
#' which is often used when several gravity models are compared, instead of adding number \code{1} to the dependent 
#' variable as done in \code{\link[gravity]{tobit}}, the constant added to the data is estimated and interpreted as a 
#' threshold.
#'
#' When taking the log of the gravity equation flows equal to zero constitute a problem as their
#' log is not defined. Therefore, a constant is added to the flows.
#'
#' Compared to the usual ET-Tobit approaches, in this package, the estimation
#' of the threshold is done before the other parameters are estimated.
#'
#' We follow \insertCite{Carson2007;textual}{gravity}, who show that taking the minimum
#' positive flow value as an estimate of the threshold is super-consistent and that
#' using this threshold estimate ensures that the parameter MLE are asymptotically normal with
#' the asymptotic variance identical to the variance achieved when the threshold is known. Hence, first
#' the threshold is estimated as the minimum positive flow. This threshold is added to the flow variable,
#' it is logged afterwards and taken as the dependent variable.
#'
#' The Tobit estimation is then conducted using the
#' \code{\link[censReg]{censReg}} function and setting the lower bound
#' equal to the log of the minimum positive flow value which was added to all
#' observations.
#'
#' A Tobit regression represents a combination of a binary and a
#' linear regression. This procedure has to be taken into consideration when
#' interpreting the estimated coefficients.
#'
#' The marginal effects of an explanatory variable on the expected value of
#' the dependent variable equals the product of both the probability of the
#' latent variable exceeding the threshold and the marginal effect of the
#' explanatory variable of the expected value of the latent variable.
#'
#' For a more elaborate Tobit function, see \code{\link[gravity]{ek_tobit}}
#' for the Eaton and Kortum (2001) Tobit model where each zero trade volume
#' is assigned a country specific interval with the upper
#' bound equal to the minimum positive trade level of the respective
#' importing country.
#'
#' The function is designed for cross-sectional data, but can be extended to panel data using the
#' \code{\link[censReg]{censReg}} function.
#'
#' A robust estimations is not implemented to the present
#' as the \code{\link[censReg]{censReg}} function is not
#' compatible with the \code{\link[sandwich]{vcovHC}} function.
#'
#' @param dependent_variable (Type: character) name of the dependent variable. Following 
#' \insertCite{Carson2007;textual}{gravity}, the smallest positive flow value is used as an estimate of the threshold, this value is is added to the \code{dependent_variable},
#' the result is logged and taken as the dependent variable in the Tobit estimation with
#' lower bound equal to the log of the smallest possible flow value.
#'
#' @param distance (Type: character) name of the distance variable that should be taken as the key independent variable 
#' in the estimation. The distance is logged automatically when the function is executed.
#'
#' @param additional_regressors (Type: character) names of the additional regressors to include in the model (e.g. a dummy
#' variable to indicate contiguity). Unilateral metric variables such as GDP should be inserted via the arguments 
#' \code{income_origin} and \code{income_destination}.
#'
#' Write this argument as \code{c(contiguity, common currency, ...)}. By default this is set to \code{NULL}.
#'
#' @param data (Type: data.frame) the dataset to be used.
#'
#' @param ... Additional arguments to be passed to the function.
#' 
#' @references
#' For more information on gravity models, theoretical foundations and
#' estimation methods in general see
#'
#' \insertRef{Anderson1979}{gravity}
#'
#' \insertRef{Anderson2001}{gravity}
#'
#' \insertRef{Anderson2010}{gravity}
#'
#' \insertRef{Baier2009}{gravity}
#'
#' \insertRef{Baier2010}{gravity}
#' 
#' \insertRef{Feenstra2002}{gravity}
#'
#' \insertRef{Head2010}{gravity}
#'
#' \insertRef{Head2014}{gravity}
#'
#' \insertRef{Santos2006}{gravity}
#'
#' and the citations therein.
#'
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#'
#' For estimating gravity equations using panel data see
#'
#' \insertRef{Egger2003}{gravity}
#'
#' \insertRef{Gomez-Herrera2013}{gravity}
#'
#' and the references therein.
#'
#' @examples
#' # Example for CRAN checks:
#' # Executable in < 5 sec
#' library(dplyr)
#' data("gravity_no_zeros")
#' 
#' # Choose 5 countries for testing
#' countries_chosen <- c("AUS", "CHN", "GBR", "BRA", "CAN")
#' grav_small <- filter(gravity_no_zeros, iso_o %in% countries_chosen)
#' 
#' grav_small <- grav_small %>%
#'   mutate(
#'     flow = ifelse(flow < 5, 0, flow), # cutoff for testing purposes
#'     lgdp_o = log(gdp_o),
#'     lgdp_d = log(gdp_d)
#'   )
#' 
#' fit <- et_tobit(
#'   dependent_variable = "flow",
#'   distance = "distw",
#'   additional_regressors = c("rta", "lgdp_o", "lgdp_d"),
#'   data = grav_small
#' )
#'
#' @return
#' The function returns the summary of the estimated gravity model as a
#' \code{\link[censReg]{censReg}}-object.
#'
#' @seealso \code{\link[censReg]{censReg}}, \code{\link[gravity]{et_tobit}}
#'
#' @export

et_tobit <- function(dependent_variable,
                     distance,
                     additional_regressors = NULL,
                     data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))

  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)

  stopifnot(is.character(distance), distance %in% colnames(data), length(distance) == 1)

  if (!is.null(additional_regressors)) {
    stopifnot(is.character(additional_regressors), all(additional_regressors %in% colnames(data)))
  }

  # Discarding unusable observations ----------------------------------------
  d <- data %>%
    filter_at(vars(!!sym(distance)), any_vars(. > 0)) %>%
    filter_at(vars(!!sym(distance)), any_vars(is.finite(.)))

  # Transforming data, logging distances ---------------------------------------
  d <- d %>%
    mutate(
      dist_log = log(!!sym(distance))
    )

  # Transforming data, logging flows -------------------------------------------
  flow_min_log <- filter_at(d, vars(!!sym(dependent_variable)), any_vars(. > 0))

  d <- d %>%
    mutate(
      y_log_et = ifelse(!!sym(dependent_variable) > 0, log(!!sym(dependent_variable)), NA)
    )

  # Transforming data, logging flows, distances --------------------------------
  d <- d %>%
    mutate(
      y2 = ifelse(!!sym(dependent_variable) > 0, !!sym(dependent_variable), NA),
      y2_log = log(!!sym("y2"))
    )

  y2min <- min(select(d, !!sym("y2")), na.rm = TRUE)
  y2min_log <- log(y2min)

  d <- d %>%
    rowwise() %>%
    mutate(y_cens_log_et = log(sum(!!sym(dependent_variable), y2min, na.rm = TRUE))) %>%
    ungroup()

  # Model -------------------------------------------------------------------
  if (!is.null(additional_regressors)) {
    vars <- paste(c("dist_log", additional_regressors), collapse = " + ")
  } else {
    vars <- "dist_log"
  }
  
  form <- stats::as.formula(paste("y_cens_log_et", "~", vars))
  
  model_et_tobit <- censReg::censReg(
    formula = form,
    left = y2min_log,
    right = Inf,
    data = d,
    start = rep(0, 3 + length(additional_regressors)),
    method = "BHHH"
  )
  
  model_et_tobit$call <- form
  return(model_et_tobit)
}
