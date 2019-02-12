#' @title Left-censored Tobit model with known threshold
#'
#' @description \code{tobit} estimates gravity models in their additive form
#' by conducting a left-censored regression, which, after adding the
#' constant \code{1} to the dependent variable, utilizes \code{log(1) = 0}
#' as the censoring value.
#'
#' @details \code{tobit} represents the left-censored tobit \insertCite{Tobin1958;textual}{gravity}
#' approach utilizing a known censoring threshold
#' which is often used when several gravity models are compared.
#'
#' When taking the log of the gravity equation flows equal to zero
#' constitute a problem as their log is not defined.
#'
#' Therefore, in the execution of the function the number \code{1}
#' is added to all flows and the \code{log(flows+1)} is
#' taken as the dependent variable.
#'
#' The tobit estimation is conducted using the \code{\link[censReg]{censReg}}
#' function and setting the lower bound equal to \code{0} as
#' \code{log(1)=0} represents the smallest flows in the transformed
#' variable.
#'
#' A tobit regression represents a combination of a binary and a
#' linear regression.
#'
#' This procedure has to be taken into consideration when
#' interpreting the estimated coefficients.
#'
#' The marginal effects of an explanatory variable on the expected value of
#' the dependent variable equals the product of both the probability of the
#' latent variable exceeding the threshold and the marginal effect of the
#' explanatory variable of the expected value of the latent variable.
#'
#' The function is designed for cross-sectional data,
#' but can be easily extended to panel data using the
#' \code{\link[censReg]{censReg}} function.
#'
#' A robust estimations is not implemented to the present
#' as the \code{\link[censReg]{censReg}} function is not
#' compatible with the \code{\link[sandwich]{vcovHC}} function.
#'
#' For a more elaborate Tobit function, see \code{\link[gravity]{ek_tobit}}
#' for the Eaton and Kortum (2001) Tobit model where each zero trade volume
#' is assigned a country specific interval with the upper
#' bound equal to the minimum positive trade level of the respective
#' importing country.
#'
#' @param dependent_variable (Type: character) name of the dependent variable. The number \code{1} is added and the 
#' transformed variable is logged and taken as the dependent variable in the tobit estimation with lower bound
#' equal to \code{0} as \code{log(1) = 0} represents the smallest flows
#' in the transformed variable.
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
#' @param added_constant (Type: numeric) the constant to be added to the dependent variable. The default value
#' is \code{1}. The minimum of \code{log(y + added_constant)} is taken as the
#' left boundary in the Tobit model.
#'
#' In the often used case of \code{added_constant = 1}, the dependent variable is left-censored at value \code{0}
#' as \code{log(1) = 0}.
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
#'     lgdp_o = log(gdp_o),
#'     lgdp_d = log(gdp_d)
#'   )
#' 
#' fit <- tobit(
#'   dependent_variable = "flow",
#'   distance = "distw",
#'   additional_regressors = c("rta", "lgdp_o", "lgdp_d"),
#'   added_constant = 1,
#'   data = grav_small
#' )
#'
#' @return
#' The function returns the summary of the estimated gravity model as a
#' \code{\link[censReg]{censReg}}-object.
#'
#' @seealso \code{\link[censReg]{censReg}}
#'
#' @export

tobit <- function(dependent_variable,
                  distance,
                  additional_regressors = NULL,
                  added_constant = 1,
                  data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))

  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)

  stopifnot(is.character(distance), distance %in% colnames(data), length(distance) == 1)

  if (!is.null(additional_regressors)) {
    stopifnot(is.character(additional_regressors), all(additional_regressors %in% colnames(data)))
  }

  stopifnot(is.numeric(added_constant), length(added_constant) == 1)

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
  d <- d %>%
    rowwise() %>%
    mutate(y_cens_log_tobit = log(sum(!!sym(dependent_variable), added_constant, na.rm = TRUE))) %>%
    ungroup()

  ypc_log_min <- min(d %>% select(!!sym("y_cens_log_tobit")), na.rm = TRUE)

  # Model ----------------------------------------------------------------------
  if (!is.null(additional_regressors)) {
    vars <- paste(c("dist_log", additional_regressors), collapse = " + ")
  } else {
    vars <- "dist_log"
  }
  
  form <- stats::as.formula(paste("y_cens_log_tobit", "~", vars))
  
  model_tobit <- censReg::censReg(
    formula = form,
    left = ypc_log_min,
    right = Inf,
    data = d,
    start = rep(0, 3 + length(additional_regressors)),
    method = "BHHH"
  )

  model_tobit$call <- form
  return(model_tobit)
}
