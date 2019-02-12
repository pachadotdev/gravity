#' @title Apply logarithm to distance column
#'
#' @description \code{log_distance} creates a new \code{dist_log}
#' column in log scale by taking the original distance column in the data
#'
#' @param data (Type: data.frame) the dataset to be used.
#'
#' @param distance The distance column
#'
#' @examples
#' log_distance(gravity_zeros, "distw")
#' @return
#' The function returns the summary of the estimated gravity model as an
#' \code{\link[stats]{lm}}-object.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}},
#' \code{\link[sandwich]{vcovHC}}
#'
#' @export

log_distance <- function(data, distance) {
  data %>%
    mutate(
      dist_log = log(!!sym(distance))
    )
}
