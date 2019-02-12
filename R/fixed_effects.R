#' @title Fixed Effects
#'
#' @description \code{fixed_effects} estimates gravity models via
#' OLS and fixed effects for the countries of origin and destination.
#'
#' @details To account for MR terms, \insertCite{Feenstra2002;textual}{gravity} and 
#' \insertCite{Feenstra2004;textual}{gravity} propose to use
#' importer and exporter fixed effects. Due to the use of these effects, all
#' unilateral influences such as GDPs can no longer be estimated.
#'
#' A disadvantage of the use of \code{fixed_effects} is that, when applied to
#' panel data, the number of country-year or country-pair fixed effects can be
#' too high for estimation. In addition, no comparative statistics are
#' possible with \code{fixed_effects} as the Multilateral Resistance terms are not estimated
#' explicitly. Nevertheless, \insertCite{Head2014;textual}{gravity} highlight the importance of
#' the use of fixed effects.
#'
#' By including country specific fixed effects, all monadic effects
#' are captured, including Multilateral Resistance terms.
#' Therefore, no other unilateral variables such as GDP can be
#' included as independent variables in the estimation.
#'
#' \code{fixed_effects} estimation can be used for both, cross-sectional as well as
#' panel data.
#'
#' Nonetheless, the function is designed to be consistent with the
#' Stata code for cross-sectional data provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' The function \code{fixed_effects} was therefore tested for
#' cross-sectional data. Its up to the user to ensure that the functions can be applied
#' to panel data.
#'
#' Also, note that by including bilateral fixed effects such as country-pair
#' effects, the coefficients of time-invariant observables such as distance
#' can no longer be estimated.
#'
#' Depending on the specific model, the code of the
#' respective function might have to be changed in order to exclude the distance
#' variable from the estimation.
#'
#' At the very least, the user should take special
#' care with respect to the meaning of the estimated coefficients and variances
#' as well as the decision about which effects to include in the estimation.
#'
#' When using panel data, the parameter and variance estimation of the models
#' may have to be changed accordingly.
#'
#' For a comprehensive overview of gravity models for panel data
#' see \insertCite{Egger2003;textual}{gravity}, \insertCite{Gomez-Herrera2013;textual}{gravity} and
#' \insertCite{Head2010;textual}{gravity} as well as the references therein.
#'
#' @param dependent_variable (Type: character) name of the dependent variable. This variable is logged and then used as 
#' the dependent variable in the estimation.
#'
#' @param distance (Type: character) name of the distance variable that should be taken as the key independent variable 
#' in the estimation. The distance is logged automatically when the function is executed.
#'
#' @param additional_regressors (Type: character) names of the additional regressors to include in the model (e.g. a dummy
#' variable to indicate contiguity). Unilateral metric variables such as GDPs can be added but those variables have to be 
#' logged first.
#'
#' Write this argument as \code{c(contiguity, common currency, ...)}. By default this is set to \code{NULL}.
#'
#' @param code_origin (Type: character) country of origin variable (e.g. ISO-3 country codes). The variables are grouped 
#' using this parameter.
#'
#' @param code_destination (Type: character) country of destination variable (e.g. country ISO-3 codes). The variables are 
#' grouped using this parameter.
#'
#' @param robust (Type: logical) whether robust fitting should be used. By default this is set to \code{FALSE}.
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
#' fit <- fixed_effects(
#'   dependent_variable = "flow",
#'   distance = "distw",
#'   additional_regressors = c("rta", "comcur", "contig"),
#'   code_origin = "iso_o",
#'   code_destination = "iso_d",
#'   robust = FALSE,
#'   data = grav_small
#' )
#'
#' @return
#' The function returns the summary of the estimated gravity model as an
#' \code{\link[stats]{lm}}-object.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}},
#' \code{\link[sandwich]{vcovHC}}
#'
#' @export

fixed_effects <- function(dependent_variable,
                          distance,
                          additional_regressors = NULL,
                          code_origin,
                          code_destination,
                          robust = FALSE,
                          data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(robust))

  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)

  stopifnot(is.character(distance), distance %in% colnames(data), length(distance) == 1)

  if (!is.null(additional_regressors)) {
    stopifnot(is.character(additional_regressors), all(additional_regressors %in% colnames(data)))
  }

  valid_origin <- data %>% select(code_origin) %>% distinct() %>% as_vector()
  valid_destination <- data %>% select(code_destination) %>% distinct() %>% as_vector()
  
  stopifnot(is.character(code_origin), code_origin %in% colnames(data), length(code_origin) == 1)
  stopifnot(is.character(code_destination), code_destination %in% colnames(data), length(code_destination) == 1)

  # Discarding unusable observations ----------------------------------------
  d <- data %>%
    filter_at(vars(!!sym(distance)), any_vars(. > 0)) %>%
    filter_at(vars(!!sym(distance)), any_vars(is.finite(.))) %>%
    filter_at(vars(!!sym(dependent_variable)), any_vars(. > 0)) %>%
    filter_at(vars(!!sym(dependent_variable)), any_vars(is.finite(.)))

  # Transforming data, logging flows and distances --------------------------
  d <- data %>%
    mutate(
      dist_log = log(!!sym(distance)),
      y_log_fe = log(!!sym(dependent_variable))
    )

  # Model -------------------------------------------------------------------
  if (!is.null(additional_regressors)) {
    vars <- paste(c("dist_log", additional_regressors, code_origin, code_destination), collapse = " + ")
  } else {
    vars <- "dist_log"
  }
  
  form <- stats::as.formula(paste("y_log_fe", "~", vars))

  if (robust == TRUE) {
    model_fe <- MASS::rlm(form, data = d)
  } else {
    model_fe <- stats::lm(form, data = d)
  }
  
  model_fe$call <- form
  return(model_fe)
}
