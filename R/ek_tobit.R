#' @title \insertCite{Eaton2001;textual}{gravity} Tobit model (EK Tobit)
#'
#' @description \code{ek_tobit} estimates gravity models in their additive form
#' by conducting a censored regression.
#'
#' @details \code{ek_tobit} represents the \insertCite{Eaton2001;textual}{gravity} Tobit model where each country
#' is assigned specific censoring bounds.
#'
#' When taking the log of the gravity equation flows equal to zero
#' constitute a problem as their log is not defined. Therefore, in \code{ek_tobit} all values of
#' the dependent variable are redefined as intervals.
#'
#' The positive observations have both interval bounds equal to their original value.
#'
#' For zero flows the interval is left open. The right border of the interval is set to the
#' log of the minimum positive trade flow of the respective importing country.
#'
#' The defined data object of class \code{\link[survival]{Surv}} is then inserted
#' in \code{\link[survival]{survreg}} for the parameter estimation.
#'
#' \code{ek_tobit} is designed to be consistent with the Stata code provided at
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' For other Tobit functions, see \code{\link[gravity]{tobit}}
#' for a simple Tobit model where number \code{1} is added to all observations
#' and \code{\link[gravity]{et_tobit}} for the \insertCite{Eaton1995;textual}{gravity}
#' threshold Tobit model where instead of simply adding number \code{1}
#' to the data the threshold is estimated.
#'
#' The function is designed for cross-sectional data, but can be extended to panel data using the
#' \code{\link[survival]{survreg}} function.
#'
#' @param dependent_variable (Type: character) name of the dependent variable. This variable is logged and then used as the 
#' dependent variable in the estimation. As the log of zero is not defined, all flows equal to zero are replaced by a left 
#' open interval with the logged minimum trade flow of the respective importing country as right border.
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
#' grav_small <- grav_small %>%
#'   mutate(
#'     flow = ifelse(flow < 5, 0, flow), # cutoff for testing purposes
#'     lgdp_o = log(gdp_o),
#'    lgdp_d = log(gdp_d)
#'   )
#' 
#' fit <- ek_tobit(
#'   dependent_variable = "flow",
#'   distance = "distw",
#'   additional_regressors = c("distw", "rta", "lgdp_o", "lgdp_d"),
#'   code_destination = "iso_d",
#'   robust = FALSE,
#'   data = grav_small
#' )
#'
#' @return
#' The function returns the summary of the estimated gravity model as a
#' \code{\link[survival]{survreg}}-object.
#'
#' @seealso \code{\link[survival]{Surv}}, \code{\link[survival]{survreg}}, \code{\link[gravity]{tobit}}
#'
#' @export

ek_tobit <- function(dependent_variable,
                     distance,
                     additional_regressors = NULL,
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

  valid_destination <- data %>% select(code_destination) %>% distinct() %>% as_vector()
  
  stopifnot(is.character(code_destination), code_destination %in% colnames(data), length(code_destination) == 1)
  
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
    mutate(
      y_log_ek = ifelse(!!sym(dependent_variable) > 0, log(!!sym(dependent_variable)), NA)
    )

  # Minimum flows -----------------------------------------------------------
  d <- d %>%
    group_by(!!sym(code_destination)) %>%
    mutate(
      exportmin = min(!!sym("y_log_ek"), na.rm = TRUE)
    )

  # Transforming censored variables -----------------------------------------
  d <- d %>%
    mutate(
      flows_ek1 = ifelse(!!sym(dependent_variable) > 0, !!sym("y_log_ek"), -Inf)
    ) %>%
    mutate(
      flows_ek2 = ifelse(!!sym(dependent_variable) > 0, !!sym("flows_ek1"), !!sym("exportmin"))
    ) %>%
    ungroup()

  # Response variable -------------------------------------------------------
  f1 <- d %>% select(!!sym("flows_ek1")) %>% as_vector()
  f2 <- d %>% select(!!sym("flows_ek2")) %>% as_vector()

  y_cens_log_ek <- survival::Surv(f1, f2, type = "interval2") %>% as_vector()

  # Model -------------------------------------------------------------------
  if (!is.null(additional_regressors)) {
    vars <- paste(c("dist_log", additional_regressors), collapse = " + ")
  } else {
    vars <- "dist_log"
  }
  
  form <- stats::as.formula(paste("y_cens_log_ek", "~", vars))
  
  model_ek_tobit <- survival::survreg(
    form, 
    data = d, 
    dist = "gaussian", 
    robust = robust
  )
  
  model_ek_tobit$call <- form
  return(model_ek_tobit)
}
