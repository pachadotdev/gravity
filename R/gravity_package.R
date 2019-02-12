#' gravity exported operators and S3 methods
#' The following functions are imported and then re-exported
#' from the gravity package to avoid listing Depends of gravity.
#' @importFrom dplyr filter filter_at select mutate group_by ungroup
#'   row_number left_join ends_with vars any_vars all_vars rowwise
#'   rename distinct one_of
#' @importFrom tidyr gather spread
#' @importFrom purrr as_vector
#' @importFrom rlang sym syms
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom survival Surv
#' @importFrom multiwayvcov cluster.vcov
#' @importFrom lmtest coeftest
#' @importFrom stats lm as.formula glm
#' @importFrom censReg censReg
#' @importFrom utils globalVariables
#' @importFrom Rdpack reprompt
#' @name gravity-exports
#' @keywords internal

globalVariables(c("."))
