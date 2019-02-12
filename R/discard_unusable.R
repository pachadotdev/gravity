#' @title Discard non-positive and/or non-finite observations in datasets
#'
#' @description \code{discard_unusable} drops observations that cannot be used
#' with models that convert columns to log scale, and therefore requiere non-negative
#' and finite observations.
#'
#' @description Consider that some of the functions within this package will drop
#' observations when required and it is not requiered to be run before fitting a model.
#'
#' @param data (Type: data.frame) the dataset to be used.
#'
#' @param columns The columns to be cleaned (e.g. \code{c("flow", "distw"))} in
#' the case of \code{ddm} when used with the example dataset \code{gravity_zeros})
#'
#' @examples
#' discard_unusable(gravity_zeros, "flow")
#' discard_unusable(gravity_zeros, c("flow", "distw"))
#' @return
#' The function returns the summary of the estimated gravity model as an
#' \code{\link[stats]{lm}}-object.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}},
#' \code{\link[sandwich]{vcovHC}}
#'
#' @export

discard_unusable <- function(data, columns) {
  data %>%
    filter_at(columns, all_vars(. > 0)) %>%
    filter_at(columns, all_vars(is.finite(.)))
}
