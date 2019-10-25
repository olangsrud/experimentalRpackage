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
#' @importFrom stats as.formula 
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



#' Sequence within unique values 
#'
#' @param x vector 
#' @param sortdata matrix or vector to determine sequence order
#'
#' @return integer vector
#' @importFrom SSBtools  SortRows
#' @export
#'
#' @examples
#' # 1:4 within A and 1:2 within B
#' UniqueSeq(c("A", "A", "B", "B", "A", "A"))
#' 
#' # Ordered differently
#' UniqueSeq(c("A", "A", "B", "B", "A", "A"), c(4, 5, 20, 10, 3, 0))
UniqueSeq <- function(x, sortdata = matrix(1L, length(x), 0)) {
  ix <- SortRows(cbind(x, sortdata), index.return = TRUE)
  nr <- seq_len(length(x))
  sortx <- x[ix]
  (1L + nr - nr[match(sortx, sortx)])[order(ix)]
}
