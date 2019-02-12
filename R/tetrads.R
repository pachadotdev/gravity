#' @title Tetrads
#'
#' @description \code{tetrads} estimates gravity models
#' by taking the ratio of the ratio of flows.
#'
#' @details \code{tetrads} is an estimation method for gravity models
#' developed by \insertCite{Head2010;textual}{gravity}.
#'
#' The function \code{tetrads} utilizes the multiplicative form of the
#' gravity equation. After choosing a reference exporter \code{A} and
#' importer \code{B} one can eliminate importer and exporter fixed effects
#' by taking the ratio of ratios.
#'
#' Only those exporters trading with the
#' reference importer and importers trading with the reference exporter will
#' remain for the estimation. Therefore, reference countries should
#' preferably be countries which trade with every other country in the dataset.
#'
#' After restricting the data in this way, \code{tetrads} estimates the gravity
#' equation in its additive form by OLS.
#'
#' By taking the ratio of ratios, all monadic effects diminish, hence no
#' unilateral variables such as GDP can be included as independent variables.
#'
#' \code{tetrads} estimation can be used for both, cross-sectional as well as
#' panel data. Nonetheless, the function is designed to be consistent with the
#' Stata code for cross-sectional data provided on the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' The function \code{tetrads} was therefore tested for cross-sectional data.
#'
#' If tetrads is used for panel data, the user may have to drop distance as an
#' independent variable as time-invariant effects drop.
#'
#' For applying \code{tetrads} to panel data see \insertCite{Head2010;textual}{gravity}.
#'
#' @param dependent_variable (Type: character) name of the dependent variable. This variable is logged and then used as 
#' the dependent variable in the estimation.
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
#' @param code_origin (Type: character) country of origin variable (e.g. ISO-3 country codes). The variables are grouped 
#' using this parameter.
#'
#' @param code_destination (Type: character) country of destination variable (e.g. country ISO-3 codes). The variables are 
#' grouped using this parameter.
#'
#' @param filter_origin (Type: character) Reference exporting country.
#'
#' @param filter_destination (Type: character) Reference importing country.
#'
#' @param multiway (Type: logical) in case \code{multiway = TRUE}, the
#' \code{\link[multiwayvcov]{cluster.vcov}} function is used for estimation following
#' \insertCite{Cameron2011;textual}{gravity} multi-way clustering of
#' variance-covariance matrices. The default value is set to \code{TRUE}.
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
#' fit <- tetrads(
#'   dependent_variable = "flow",
#'   distance = "distw",
#'   additional_regressors = "rta",
#'   code_origin = "iso_o",
#'   code_destination = "iso_d",
#'   filter_origin = countries_chosen[1],
#'   filter_destination = countries_chosen[2],
#'   data = grav_small
#' )
#'
#' @return
#' The function returns the summary of the estimated gravity model as an
#' \code{\link[stats]{lm}}-object.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}},
#'
#' @export

tetrads <- function(dependent_variable,
                    distance,
                    additional_regressors,
                    code_origin,
                    code_destination,
                    filter_origin = NULL,
                    filter_destination = NULL,
                    multiway = FALSE,
                    data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))

  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)

  stopifnot(is.character(distance), distance %in% colnames(data), length(distance) == 1)

  if (!is.null(additional_regressors)) {
    stopifnot(is.character(additional_regressors), all(additional_regressors %in% colnames(data)))
  }

  stopifnot(is.character(code_origin), code_origin %in% colnames(data), length(code_origin) == 1)
  stopifnot(is.character(code_destination), code_destination %in% colnames(data), length(code_destination) == 1)

  valid_origin <- data %>% select(code_origin) %>% distinct() %>% as_vector()
  valid_destination <- data %>% select(code_destination) %>% distinct() %>% as_vector()
  
  stopifnot(is.character(filter_origin), filter_origin %in% valid_origin, length(filter_origin) == 1)
  stopifnot(is.character(filter_destination), filter_destination %in% valid_destination, length(filter_destination) == 1)

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
      y_log_tetrads = log(!!sym(dependent_variable))
    )

  # Truncating dataset to reference importer and exporter partners -------------
  d2_filter_d <- d %>%
    filter_at(vars(!!sym(code_destination)), any_vars(. == filter_destination)) %>%
    select(!!sym(code_origin)) %>%
    distinct() %>%
    as_vector()

  d2 <- d %>%
    filter_at(vars(!!sym(code_origin)), any_vars(. %in% d2_filter_d))

  d2_filter_o <- d2 %>%
    filter_at(vars(!!sym(code_origin)), any_vars(. == filter_origin)) %>%
    select(!!sym(code_destination)) %>%
    distinct() %>%
    as_vector()

  d2 <- d2 %>%
    filter_at(vars(!!sym(code_destination)), any_vars(. %in% d2_filter_o))

  # Taking ratios, ratk --------------------------------------------------------
  d2 <- left_join(
    d2,
    d2 %>%
      filter_at(vars(!!sym(code_destination)), any_vars(. == filter_destination)) %>%
      select(!!sym(code_origin), y_log_tetrads_d = !!sym("y_log_tetrads"), dist_log_d = !!sym("dist_log")),
    by = code_origin
  ) %>%
    mutate(
      lXinratk = !!sym("y_log_tetrads") - !!sym("y_log_tetrads_d"),
      ldistratk = !!sym("dist_log") - !!sym("dist_log_d")
    ) %>%
    select(-!!sym("y_log_tetrads_d"), -!!sym("dist_log_d"))

  # Taking ratios, ratk, for the other independent variables -------------------
  d2 <- d2 %>%
    select(c(!!sym(code_origin), !!sym(code_destination), !!sym(distance), !!!syms(additional_regressors))) %>%
    gather(!!sym("key"), !!sym("value"), -!!sym(code_origin), -!!sym(code_destination)) %>%
    left_join(
      d2 %>%
        filter_at(vars(!!sym(code_destination)), any_vars(. == filter_destination)) %>%
        select(c(!!sym(code_origin), !!sym(distance), !!!syms(additional_regressors))) %>%
        gather(!!sym("key"), !!sym("value"), -!!sym(code_origin)),
      by = c(code_origin, "key")
    ) %>%
    mutate(value = !!sym("value.x") - !!sym("value.y")) %>%
    select(c(!!sym(code_origin), !!sym(code_destination), "key", "value")) %>%
    spread(!!sym("key"), !!sym("value")) %>%
    left_join(d2 %>% select(-one_of(distance, additional_regressors)), by = c(code_origin, code_destination))

  # Taking the ratio of ratios, rat --------------------------------------------
  d2 <- left_join(
    d2,
    d2 %>%
      filter_at(vars(!!sym(code_origin)), any_vars(. == filter_origin)) %>%
      select(!!sym(code_destination), lXinratk_o = !!sym("lXinratk"), ldistratk_o = !!sym("ldistratk")),
    by = code_destination
  ) %>%
    mutate(
      y_log_rat = !!sym("lXinratk") - !!sym("lXinratk_o"),
      dist_log_rat = !!sym("ldistratk") - !!sym("ldistratk_o")
    ) %>%
    select(-!!sym("lXinratk_o"), -!!sym("ldistratk_o"))

  # Taking the ratio of ratios, rat, for the other independent variables -------
  d2 <- d2 %>%
    select(c(!!sym(code_origin), !!sym(code_destination), !!sym(distance), !!!syms(additional_regressors))) %>%
    gather(!!sym("key"), !!sym("value"), -!!sym(code_origin), -!!sym(code_destination)) %>%
    left_join(
      d2 %>%
        filter_at(vars(!!sym(code_origin)), any_vars(. == filter_origin)) %>%
        select(c(!!sym(code_destination), additional_regressors)) %>%
        gather(!!sym("key"), !!sym("value"), -!!sym(code_destination)),
      by = c(code_destination, "key")
    ) %>%
    mutate(
      key = paste0(!!sym("key"), "_rat"),
      value = !!sym("value.x") - !!sym("value.y")
    ) %>%
    select(c(!!sym(code_origin), !!sym(code_destination), "key", "value")) %>%
    spread(!!sym("key"), !!sym("value")) %>%
    left_join(d2 %>% select(-one_of(distance, additional_regressors)), by = c(code_origin, code_destination)) %>%
    select(ends_with("_rat"), !!sym(code_origin), !!sym(code_destination))

  # Model ----------------------------------------------------------------------
  if (!is.null(additional_regressors)) {
    additional_regressors <- paste0(additional_regressors, "_rat")
    vars <- paste(c("dist_log_rat", additional_regressors), collapse = " + ")
  } else {
    vars <- "dist_log_rat"
  }
  
  form <- stats::as.formula(
    sprintf("y_log_rat ~ %s", vars)
  )
  
  form2 <- as.formula(sprintf("~ %s + %s", code_origin, code_destination))
  
  model_tetrads <- stats::lm(form, data = d2)

  # Multiway Variance-Covariance -----------------------------------------------
  if (multiway == TRUE) {
    model_tetrads_multiway <- lmtest::coeftest(
      model_tetrads,
      vcov = sandwich::vcovHC(model_tetrads, "HC1")
    )
  }
  
  # Return ---------------------------------------------------------------------
  if (multiway == FALSE) {
    model_tetrads$call <- form
    return(model_tetrads)
  } else {
    #model_tetrads_multiway$call <- form
    return(model_tetrads_multiway)
  }
}
