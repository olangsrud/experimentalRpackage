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



#' Model formula from vector of terms where intercept is specified as \code{"1"}
#'
#' @param x Character vector
#' @param InterceptAs1 Whether to specify intercept as \code{"1"}
#' @param env Parameter to \code{\link{as.formula}}
#'
#' @return formula
#' @export
#'
#' @examples
#' Terms2formula(c("a", "b:c"))
#' Terms2formula(c("a", "b:c"), FALSE)
#' Terms2formula(c("a", "b:c", "1"))
Terms2formula <- function(x, InterceptAs1 = TRUE, env = parent.frame()) {
  x1 <- x == "1"
  if (InterceptAs1) 
    x <- x[!x1]
  f <- as.formula(paste("~", paste(x, collapse = " + "), sep = ""), env = env)
  if (!InterceptAs1 | sum(x1)) 
    return(f)
  update(f, ~. - 1)
}