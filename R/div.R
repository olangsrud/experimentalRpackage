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



#' Model formula from vector of terms where intercept is specified as a string
#'
#' @param x Character vector
#' @param intercept String to specify intercept
#' @param env Parameter to \code{\link{as.formula}}
#'
#' @return formula
#' 
#' @importFrom stats as.formula update
#' @export
#'
#' @examples
#' Terms2formula(c("a", "b:c"))
#' Terms2formula(c("a", "b:c"), NULL)
#' Terms2formula(c("a", "b:c", "(Intercept)"))
#' Terms2formula(c("a", "b:c"), "1")
#' Terms2formula(c("a", "b:c", "1"), "1")
Terms2formula <- function(x, intercept = "(Intercept)", env = parent.frame()) {
  if (!is.null(intercept)) {
    x1 <- x == intercept
    x <- x[!x1]
    anyx1 <- any(x1)
  } else {
    anyx1 <- FALSE
  }
  if (!length(x))
    x <- "1"
  f <- as.formula(paste("~", paste(x, collapse = " + "), sep = ""), env = env)
  if (is.null(intercept) | anyx1) 
    return(f)
  update(f, ~. - 1)
}