#' @title Structural Iterated Least Squares (SILS)
#'
#' @description \code{sils} estimates gravity models via
#' Structural Iterated Least Squares and an explicit inclusion
#' of the Multilateral Resistance terms.
#'
#' @details \code{sils} is an estimation method for gravity models
#' developed by \insertCite{Head2014;textual}{gravity}.
#'
#' The function \code{sils} utilizes the relationship between the Multilateral
#' Resistance terms and the transaction costs. The parameters are estimated by
#' an iterative procedure. The function executes loops until the parameters
#' stop changing significantly.
#'
#' \code{sils} is designed to be consistent with the Stata code provided at
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' As, to our knowledge at the moment, there is no explicit literature covering
#' the estimation of a gravity equation by \code{sils} using panel data,
#' and we do not recommend to apply this method in this case.
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
#' @param maxloop (Type: numeric) maximum number of outer loop iterations.
#' The default is set to 100. There will be a warning if the iterations
#' did not converge.
#'
#' @param decimal_places (Type: numeric) number of decimal places that should not change after a new
#' iteration for the estimation to stop. The default is set to 4.
#'
#' @param robust (Type: logical) whether robust fitting should be used. By default this is set to \code{FALSE}.
#'
#' @param verbose (Type: logical) determines whether the estimated coefficients
#' of each iteration should be printed in the console. The default is set
#' to \code{FALSE}.
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
#' fit <- sils(
#'   dependent_variable = "flow",
#'   distance = "distw",
#'   additional_regressors = "rta",
#'   income_origin = "gdp_o",
#'   income_destination = "gdp_d",
#'   code_origin = "iso_o",
#'   code_destination = "iso_d",
#'   maxloop = 50,
#'   dec_places = 3,
#'   robust = FALSE,
#'   verbose = FALSE,
#'   data = grav_small
#' )
#'
#' @return
#' The function returns the summary of the estimated gravity model as an
#' \code{\link[stats]{lm}}-object. It furthermore returns the resulting coefficients for each
#' iteration.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}},
#' \code{\link[sandwich]{vcovHC}}
#'
#' @export
#'
sils <- function(dependent_variable,
                 distance,
                 additional_regressors = NULL,
                 income_origin,
                 income_destination,
                 code_origin,
                 code_destination,
                 maxloop = 100,
                 decimal_places = 4,
                 robust = FALSE,
                 verbose = FALSE,
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

  valid_origin <- data %>% select(code_origin) %>% distinct() %>% as_vector()
  valid_destination <- data %>% select(code_destination) %>% distinct() %>% as_vector()
  
  stopifnot(is.character(code_origin), code_origin %in% colnames(data), length(code_origin) == 1)
  stopifnot(is.character(code_destination), code_destination %in% colnames(data), length(code_destination) == 1)
  
  stopifnot(maxloop > 0)

  stopifnot(decimal_places >= 1)

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

  # Setting starting values for the first iteration ----------------------------
  d <- d %>%
    mutate(
      P_i = 1,
      P_j = 1
    )

  loop <- 0
  dec_point <- 1 * 10^-decimal_places
  beta_distance <- 1
  beta_distance_old <- 0
  coef_dist <- 1

  beta <- vector(length = length(additional_regressors))
  names(beta) <- additional_regressors

  for (j in 1:length(additional_regressors)) {
    beta[j] <- 1
  }

  beta_old <- vector(length = length(additional_regressors))
  names(beta_old) <- additional_regressors

  for (j in 1:length(additional_regressors)) {
    beta_old[j] <- 0
  }

  coef_additional_regressors <- data.frame(matrix(nrow = 1, ncol = length(additional_regressors)))
  coef_additional_regressors[1, ] <- 1
  names(coef_additional_regressors) <- additional_regressors

  # Begin iterations -----------------------------------------------------------
  while (loop <= maxloop &
    abs(beta_distance - beta_distance_old) > dec_point &
    prod(abs(beta - beta_old) > dec_point) == 1) {

    # Updating betas -----------------------------------------------------------
    beta_distance_old <- beta_distance
    for (j in 1:length(additional_regressors)) {
      beta_old[j] <- beta[j]
    }

    # Updating transaction costs -----------------------------------------------
    costs <- data.frame(matrix(nrow = nrow(d), ncol = length(additional_regressors)))
    for (j in 1:length(additional_regressors)) {
      costs[, j] <- beta[j] * d[additional_regressors[j]][, 1]
    }
    costs <- apply(X = costs, MARGIN = 1, FUN = sum)

    d <- d %>%
      ungroup() %>%
      mutate(
        bd = beta_distance,
        co = costs,
        t_ij = exp(!!sym("bd") * !!sym("dist_log") + !!sym("co"))
      )

    # Contraction mapping ------------------------------------------------------
    d <- d %>%
      mutate(
        P_j_old = 0,
        P_i_old = 0
      )

    j <- 1

    while (j <= maxloop &
      sum(abs(d$P_j - d$P_j_old)) > dec_point &
      sum(abs(d$P_i - d$P_i_old)) > dec_point) {
      d <- d %>%
        mutate(
          P_j_old = !!sym("P_j"),
          P_i_old = !!sym("P_i")
        )

      # Inward MR --------------------------------------------------------------
      d <- d %>%
        group_by(!!sym(code_destination)) %>%
        mutate(
          P_j = sum(!!sym("t_ij") * !!sym(income_origin) / !!sym("P_i"))
        )

      # Outward MR -------------------------------------------------------------
      d <- d %>%
        group_by(!!sym(code_origin)) %>%
        mutate(
          P_i = sum(!!sym("t_ij") * !!sym(income_destination) / !!sym("P_j"))
        )

      j <- j + 1

      if (j == maxloop) {
        warning("The inner iteration did not converge before the inner loop reached maxloop=", maxloop, " iterations")
      }
    }

    # Model --------------------------------------------------------------------
    d <- d %>%
      mutate(
        y_log_sils = log(!!sym(dependent_variable)) -
          log((!!sym(income_origin) * !!sym(income_destination)) / (!!sym("P_i") * !!sym("P_j")))
      )

    if (!is.null(additional_regressors)) {
      vars <- paste(c("dist_log", additional_regressors), collapse = " + ")
    } else {
      vars <- "dist_log"
    }
    
    form <- stats::as.formula(paste("y_log_sils", "~", vars))

    if (robust == TRUE) {
      model_sils <- MASS::rlm(form, data = d)
    } else {
      model_sils <- stats::lm(form, data = d)
    }

    # Updating coefficients ----------------------------------------------------
    beta_distance <- stats::coef(model_sils)[2]
    for (j in 1:length(additional_regressors)) {
      beta[j] <- stats::coef(model_sils)[j + 2]
    }

    coef_dist <- c(coef_dist, beta_distance)
    coef_additional_regressors <- rbind(coef_additional_regressors, rep(0, times = length(additional_regressors)))

    for (j in 1:length(additional_regressors)) {
      coef_additional_regressors[additional_regressors[j]][loop + 2, ] <- beta[j]
    }

    # Coefficients -------------------------------------------------------------
    coef_sils <- cbind(
      loop = c(1:loop), dist = as.numeric(coef_dist)[2:(loop + 1)],
      coef_additional_regressors[2:(loop + 1), ]
    )

    if (verbose == TRUE) {
      cat("This is round", loop, "\n")
      cat("The coefficients are", beta_distance, beta, "\n")
    }

    loop <- loop + 1

    if (loop == maxloop) {
      warning("The outer iteration did not converge before the outer loop reached maxloop=", maxloop, " iterations")
    }
  }

  model_sils$call <- form
  return(model_sils)
}
