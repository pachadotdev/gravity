#' @title Ordinary Least Squares (OLS)
#'
#' @description \code{ols} estimates gravity models in their traditional form
#' via Ordinary Least Squares (ols). It does not consider Multilateral
#' Resistance terms.
#'
#' @details \code{ols} estimates gravity models in their traditional, additive,
#' form via Ordinary Least Squares using the \code{lm} function.
#' Multilateral Resistance terms are not considered by this function.
#'
#' As the coefficients for the country's incomes were often found to be close to
#' unitary and unitary income elasticities are in line with some theoretical
#' foundations on international trade, it is sometimes assumed that the income
#' elasticities are equal to unity.
#'
#' In order to allow for the estimation with
#' and without the assumption of unitary income elasticities, the option
#' \code{uie} is built into \code{ols} with the default set to \code{FALSE}.
#'
#' \code{ols} estimation can be used for both, cross-sectional and
#' panel data. Nonetheless, the function is designed to be consistent with the
#' Stata code for cross-sectional data provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' The function \code{ols} was therefore tested for cross-sectional data. For the use with panel data
#' no tests were performed. Therefore, it is up to the user to ensure that the functions can be applied
#' to panel data.
#'
#' Depending on the panel dataset and the variables -
#' specifically the type of fixed effects -
#' included in the model, it may easily occur that the model is not computable.
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
#' @param dependent_variable (Type: character) name of the dependent variable.
#'
#' If \code{uie = TRUE} the dependent variable is divided by the product of
#' unilateral incomes (e.g. \code{income_origin} and \code{income_destination}) and logged afterwards.
#'
#' If \code{uie=FALSE} the dependent variable is logged directly. The transformed variable is then used as
#' the dependent variable and the logged income variables are used as independent variables in the
#' estimation.
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
#' @param income_origin (Type: character) origin income variable (e.g. GDP) in the dataset.
#'
#' @param income_destination (Type: character) destination income variable (e.g. GDP) in the dataset.
#'
#' @param code_origin (Type: character) country of origin variable (e.g. ISO-3 country codes). The variables are grouped 
#' using this parameter.
#'
#' @param code_destination (Type: character) country of destination variable (e.g. country ISO-3 codes). The variables are 
#' grouped using this parameter.
#'
#' @param uie (Type: logical) Dtermines whether the
#' parameters are to be estimated assuming unitary income elasticities. The default value is set
#' to \code{FALSE}.
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
#' fit <- ols(
#'   dependent_variable = "flow",
#'   distance = "distw",
#'   additional_regressors = c("rta", "contig", "comcur"),
#'   income_origin = "gdp_o",
#'   income_destination = "gdp_d",
#'   code_origin = "iso_o",
#'   code_destination = "iso_d",
#'   uie = FALSE,
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

ols <- function(dependent_variable,
                distance,
                additional_regressors = NULL,
                income_origin,
                income_destination,
                code_origin,
                code_destination,
                uie = FALSE,
                robust = FALSE,
                data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(uie))
  stopifnot(is.logical(robust))

  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)

  stopifnot(is.character(distance), distance %in% colnames(data), length(distance) == 1)

  if (!is.null(additional_regressors)) {
    stopifnot(is.character(additional_regressors), all(additional_regressors %in% colnames(data)))
  }

  stopifnot(is.character(income_origin), income_origin %in% colnames(data), length(income_origin) == 1)
  stopifnot(is.character(income_destination), income_destination %in% colnames(data), length(income_destination) == 1)
  
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

  # Transforming data, logging distances ---------------------------------------
  d <- d %>%
    mutate(
      dist_log = log(!!sym(distance))
    )

  # Income elasticities --------------------------------------------------------
  if (uie == TRUE) {
    # Transforming data, logging flows -----------------------------------------
    d <- d %>%
      mutate(
        y_log_ols = log(
          !!sym(dependent_variable) / (!!sym(income_origin) * !!sym(income_destination))
        )
      ) %>%
      select(
        !!sym("y_log_ols"), !!sym("dist_log"), !!sym("additional_regressors")
      )
  }

  if (uie == FALSE) {
    # Transforming data, logging flows -----------------------------------------
    d <- d %>%
      mutate(
        y_log_ols = log(!!sym(dependent_variable)),
        inc_o_log = log(!!sym(income_origin)),
        inc_d_log = log(!!sym(income_destination))
      ) %>%
      select(
        !!sym("y_log_ols"), !!sym("inc_o_log"), !!sym("inc_d_log"), !!sym("dist_log"), !!sym("additional_regressors")
      )
  }


  # Model -------------------------------------------------------------------
  if (!is.null(additional_regressors)) {
    vars <- paste(c("dist_log", "inc_o_log", "inc_d_log", additional_regressors), collapse = " + ")
  } else {
    vars <- paste(c("dist_log", "inc_o_log", "inc_d_log"), collapse = " + ")
  }
  
  form <- stats::as.formula(paste("y_log_ols", "~", vars))
  
  if (robust == TRUE) {
    model_ols <- MASS::rlm(form, data = d)
  } else {
    model_ols <- stats::lm(form, data = d)
  }

  model_ols$call <- form
  return(model_ols)
}
