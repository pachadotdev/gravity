#' @title Bonus vetus OLS (BVU)
#'
#' @description \code{bvu} estimates gravity models via Bonus vetus OLS with simple averages.
#'
#' @details Bonus vetus OLS is an estimation method for gravity models
#' developed by \insertCite{Baier2009,Baier2010;textual}{gravity} using simple averages to center a
#' Taylor-series.
#'
#' The \code{bvu} function considers Multilateral Resistance terms and allows to
#' conduct comparative statics. Country specific effects are subdued due
#' to demeaning. Hence, unilateral variables apart from incomes cannot be included 
#' in the estimation.
#'
#' \code{bvu} is designed to be consistent with the Stata code provided at
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' As, to our knowledge at the moment, there is no explicit literature covering
#' the estimation of a gravity equation by \code{bvu} using panel data,
#' we do not recommend to apply this method in this case.
#'
#' @param dependent_variable (Type: character) name of the dependent variable. This dependent variable is 
#' divided by the product of unilateral incomes such (i.e. \code{income_origin} and \code{income_destination}) 
#' and logged afterwards.
#'
#' @param distance (Type: character) name of the distance variable that should be taken as the key independent variable 
#' in the estimation. The distance is logged automatically when the function is executed.
#'
#' @param additional_regressors (Type: character) names of the additional regressors to include in the model (e.g. a dummy
#' variable to indicate contiguity). Unilateral metric variables such as GDP should be inserted via the arguments 
#' \code{income_origin} and \code{income_destination}. As country specific effects are subdued due to demeaning, no further unilateral variables apart from incomes can be added.
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
#' fit <- bvu(
#'   dependent_variable = "flow",
#'   distance = "distw",
#'   additional_regressors = c("rta", "contig", "comcur"),
#'   income_origin = "gdp_o",
#'   income_destination = "gdp_d",
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

bvu <- function(dependent_variable,
                distance,
                additional_regressors = NULL,
                income_origin,
                income_destination,
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

  stopifnot(is.character(income_origin), income_origin %in% colnames(data), length(income_origin) == 1)
  stopifnot(is.character(income_destination), income_destination %in% colnames(data), length(income_destination) == 1)
  
  stopifnot(is.character(code_origin), code_origin %in% names(data), length(code_origin) == 1)
  stopifnot(is.character(code_destination), code_destination %in% names(data), length(code_destination) == 1)

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

  # Transforming data, logging flows -------------------------------------------
  d <- d %>%
    mutate(
      y_log_bvu = log(
        !!sym(dependent_variable) / (!!sym(income_origin) * !!sym(income_destination))
      )
    )

  # Multilateral Resistance (MR) for distance ----------------------------------
  d <- d %>%
    group_by(!!sym(code_origin)) %>%
    mutate(mean_dist_log_1 = mean(!!sym("dist_log"), na.rm = TRUE)) %>%
    group_by(!!sym(code_destination), add = FALSE) %>%
    mutate(mean_dist_log_2 = mean(!!sym("dist_log"), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      mean_dist_log_3 = mean(!!sym("dist_log"), na.rm = TRUE),
      dist_log_mr = !!sym("dist_log") -
        (!!sym("mean_dist_log_1") + !!sym("mean_dist_log_2") - !!sym("mean_dist_log_3"))
    )

  # Multilateral Resistance (MR) for the other independent variables -----------
  if (!is.null(additional_regressors)) {
    d2 <- d %>%
      select(!!sym(code_origin), !!sym(code_destination), !!sym(distance), !!!syms(additional_regressors)) %>%
      gather(!!sym("key"), !!sym("value"), -!!sym(code_origin), -!!sym(code_destination)) %>%
      group_by(!!sym(code_origin), !!sym("key")) %>%
      mutate(mean_dist_log_1 = mean(!!sym("value"), na.rm = TRUE)) %>%
      group_by(!!sym(code_destination), !!sym("key")) %>%
      mutate(mean_dist_log_2 = mean(!!sym("value"), na.rm = TRUE)) %>%
      group_by(!!sym("key")) %>%
      mutate(
        mean_dist_log_3 = mean(!!sym("value"), na.rm = TRUE),
        dist_log_mr = !!sym("value") - (!!sym("mean_dist_log_1") + !!sym("mean_dist_log_2") - !!sym("mean_dist_log_3"))
      ) %>%
      ungroup() %>%
      mutate(key = paste0(!!sym("key"), "_mr")) %>%
      select(!!!syms(c(code_origin, code_destination, "key", "dist_log_mr"))) %>%
      spread(!!sym("key"), !!sym("dist_log_mr"))
  }

  # Model ----------------------------------------------------------------------
  if (!is.null(additional_regressors)) {
    d <- left_join(d, d2, by = c(code_origin, code_destination)) %>%
      select(!!sym("y_log_bvu"), ends_with("_mr"))
  } else {
    d <- select(d, !!sym("y_log_bvu"), ends_with("_mr"))
  }
  
  if (!is.null(additional_regressors)) {
    vars <- paste(c("dist_log_mr", paste0(additional_regressors, "_mr")), collapse = " + ")
  } else {
    vars <- "dist_log_mr"
  }
  
  form <- stats::as.formula(paste("y_log_bvu", "~", vars))

  if (robust == TRUE) {
    model_bvu <- MASS::rlm(form, data = d)
  } else {
    model_bvu <- stats::lm(form, data = d)
  }
  
  model_bvu$call <- form
  return(model_bvu)
}
