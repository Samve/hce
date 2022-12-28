#' A generic function for calculating win statistics
#'
#' @param x an object used to select a method.
#' @param ... further arguments passed to or from other methods.
#'
#' @return a data frame containing calculated values.
#' @export
#'
calcWINS <- function(x, ...) {
  UseMethod("calcWINS")
}
