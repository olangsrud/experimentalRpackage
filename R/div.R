#' Laplace noise
#'
#' @param n number of observations
#' @param b scale parameter(s)
#'
#' @return random deviates
#' @export
#'
#' @examples
#' Lap(4)
#' Lap(4, 20)
#' Lap(4, 1:4)
Lap <- function(n = 1, b = 1) {
  rexp(n, 1/b) * sample(c(-1, 1), n, replace = TRUE)
}