#' @title Non-linear Least Squares (NLS)
#'
#' @description \code{nls} estimates gravity models in their
#' multiplicative form via Nonlinear Least Squares.
#'
#' @details \code{nls} is an estimation method for gravity models
#' belonging to generalized linear models. It is estimated via \code{\link[stats]{glm}} using the gaussian 
#' distribution and a log-link.
#'
#' As the method may not lead to convergence when poor
#' starting values are used, the linear predictions, fitted values,
#' and estimated coefficients resulting from a
#' \code{\link[gravity]{ppml}} estimation are used for the arguments
#' \code{etastart}, \code{mustart}, and \code{start}.
#'
#' For similar functions, utilizing the multiplicative form via the log-link,
#' but different distributions, see \code{\link[gravity]{ppml}}, \code{\link[gravity]{gpml}},
#' and \code{\link[gravity]{nbpml}}.
#'
#' \code{nls} estimation can be used for both, cross-sectional as well as
#' panel data, but its up to the user to ensure that the functions can be applied
#' to panel data.
#'
#' Depending on the panel dataset and the variables -
#' specifically the type of fixed effects -
#' included in the model, it may easily occur that the model is not computable.
#'
#' Also, note that by including bilateral fixed effects such as country-pair
#' effects, the coefficients of time-invariant observables such as distance
#' can no longer be estimated.
#'
#' Depending on the specific model, the code of the
#' respective function may has to be changed in order to exclude the distance
#' variable from the estimation.
#'
#' At the very least, the user should take special
#' care with respect to the meaning of the estimated coefficients and variances
#' as well as the decision about which effects to include in the estimation.
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
#' logged first. Interaction terms can be added.
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
#'     lgdp_o = log(gdp_o),
#'     lgdp_d = log(gdp_d)
#'   )
#' 
#' fit <- nls(
#'   dependent_variable = "flow",
#'   distance = "distw",
#'   additional_regressors = c("rta", "lgdp_o", "lgdp_d"),
#'   data = grav_small
#' )
#'
#' @return
#' The function returns the summary of the estimated gravity model similar to a
#' \code{\link[stats]{glm}}-object.
#'
#' @seealso \code{\link[stats]{glm}}, \code{\link[lmtest]{coeftest}},
#' \code{\link[sandwich]{vcovHC}}
#'
#' @export

nls <- function(dependent_variable,
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
    filter_at(vars(!!sym(distance)), any_vars(is.finite(.))) %>%
    filter_at(vars(!!sym(dependent_variable)), any_vars(. > 0)) %>%
    filter_at(vars(!!sym(dependent_variable)), any_vars(is.finite(.)))

  # Transforming data, logging distances ---------------------------------------
  d <- d %>%
    mutate(
      dist_log = log(!!sym(distance))
    ) %>%
    select(
      !!sym("dependent_variable"), !!sym("dist_log"), !!sym("additional_regressors")
    )

  # Model ----------------------------------------------------------------------
  if (!is.null(additional_regressors)) {
    vars <- paste(c("dist_log", additional_regressors), collapse = " + ")
  } else {
    vars <- "dist_log"
  }

  form <- stats::as.formula(
    sprintf(
      "%s ~ %s",
      dependent_variable,
      vars
    )
  )

  # For nls the starting values are retrieved from the resuts of PPML
  model_PPML <- stats::glm(form, data = d, family = stats::quasipoisson(link = "log"))
  model_PPML_eta <- model_PPML$linear.predictors
  model_PPML_mu <- model_PPML$fitted.values
  model_PPML_start <- model_PPML$coefficients

  model_nls <- stats::glm(form,
    data = d, family = stats::gaussian(link = "log"),
    control = list(maxit = 200, trace = FALSE),
    etastart = model_PPML_eta, # linear predictors
    mustart = model_PPML_mu, # fitted values
    start = model_PPML_start # estimated coefficients
  )

  model_nls$call <- form
  class(model_nls) <- c(class(model_nls), "gravity_nls")
  return(model_nls)
}
