#' A modified version of RegSDC::ReduceX
#' 
#' In this modified version of \code{\link{ReduceX}}, \code{digits=NULL} is default 
#' and use of \code{\link{Z2Yhat}} and \code{\link{RoundWhole}} is removed.
#'
#' @param z Z as a matrix
#' @param x X as a matrix
#' @param y Y as a matrix
#' @param digits When non-NULL and when NULL y input, output y estimates close to whole numbers will be rounded using 
#'        \code{digits} as input to \code{\link{RoundWhole}}.
#'
#' @return A list of four elements:
#'         \item{\code{x}}{Reduced \code{x}}
#'         \item{\code{z}}{Corresponding reduced \code{z} or NULL when no \code{z} in input}
#'         \item{\code{yKnown}}{Logical vector specifying elements of y that can be found directly as elements in z}
#'         \item{\code{y}}{As \code{y} in input (not reduced) or estimated \code{y} when NULL y in input}
#'         
#' @keywords internal
#' @importFrom  methods as
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' # Same data as in the paper
#' z <- RegSDCdata("sec7z")
#' x <- RegSDCdata("sec7x")
#' y <- RegSDCdata("sec7y")  # Now z is t(x) %*% y 
#' 
#' a <- ReduceXspes(x, z, y)
#' b <- ReduceXspes(x, z)
#' d <- ReduceXspes(x, z = NULL, y)  # No z in output
#' 
#' # Identical output for x and z
#' identical(a$x, b$x)
#' identical(a$x, d$x)
#' identical(a$z, b$z)
#' 
#' # Same y in output as input
#' identical(a$y, y)
#' identical(d$y, y)
#' 
#' # Estimate of y (yHat) when NULL y input
#' # NaN in this modified version of RegSDC::ReduceX
#' b$y
#' 
#' # These elements of y can be found directly in in z
#' y[a$yKnown, , drop = FALSE]
#' # They can be found by searching for unit colSums
#' colSums(x)[colSums(x) == 1]
#' 
#' # These trivial data rows can be omitted when processing data
#' x[!a$yKnown, ]
#' # Now several columns can be omitted since zero colSums
#' colSums0 <- colSums(x[!a$yKnown, ]) == 0
#' # The resulting matrix is output from the function
#' identical(x[!a$yKnown, !colSums0], a$x)
#' 
#' # Output z can be computed from this output x
#' identical(t(a$x) %*% y[!a$yKnown, , drop = FALSE], a$z)
ReduceXspes <- function(x, z = NULL, y = NULL, digits = NULL) {   #  ReduceX <- function(x, z = NULL, y = NULL, digits = 9) {
  Z2Yhat <- function(z, x, digits = NULL) {  ### sette til 0 istenen? skal jo ikke brukes 
    matrix(NaN, nrow(x), ncol(z))
  }
  RoundWhole <- function(x, digits = NULL, onlyZeros = FALSE) {
    x
  }
  yNULL <- is.null(y)
  zNULL <- is.null(z)
  if(yNULL & zNULL)
    stop("z or y must be supplied")
  colSums_1 <- which(colSums(x) == 1)
  x1 <- x[, colSums_1, drop = FALSE]
  x1dgT <- as(x1, "dgTMatrix")
  nonDub <- x1dgT@j[x1dgT@x != 0][!duplicated(x1dgT@i[x1dgT@x != 0])] + 1L
  x1 <- x1[, nonDub, drop = FALSE]
  if (!zNULL) 
    zA <- z[colSums_1[nonDub], , drop = FALSE]
  zA1 <- matrix(1, NCOL(x1), 1)
  yKnown1 <- round(x1 %*% zA1)
  yKnown1_0 <- which(yKnown1 == 0)
  if (yNULL) {
    yHat <- x1 %*% zA
  } else {
    yHat <- y
    yHat[yKnown1_0, ] <- 0
  }
  if (!zNULL) 
    z <- z - crossprod(x, yHat)
  x <- x[yKnown1_0, , drop = FALSE]
  colSums_ok <- which(colSums(x) != 0)
  if (!zNULL) 
    z <- z[colSums_ok, , drop = FALSE]
  x <- x[, colSums_ok, drop = FALSE]
  if (yNULL) {
    if (length(yKnown1_0))
      yHat[yKnown1_0, ] <- Z2Yhat(z, x, digits = NA)
    if (!is.null(digits))
      if (!is.na(digits)) 
        yHat <- RoundWhole(yHat, digits = digits)
  } else {
    yHat <- y
  }
  list(x = x, z = z, yKnown = yKnown1 != 0, y = yHat)
}


#' Reduce input to \code{\link{Mifp}} based on zeros in z and possible zeros in yStart
#' 
#' @param x a matrix 
#' @param z a single column matrix
#' @param yStart a starting estimate of \code{y}
#'
#' @return A list of three elements:
#'         \item{\code{x}}{Reduced \code{x}}
#'         \item{\code{z}}{Corresponding reduced \code{z}}
#'         \item{\code{yKnown}}{Logical vector specifying elements of y that can be found directly as 0}
#' @keywords internal
#' 
#' @importFrom Matrix rowSums
#' 
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' # Generate data with zeros
#' z2 <- EasyData("z2")
#' z2$ant[z2$hovedint == "trygd"] <- 0
#' z2$ant[1:3] <- 0
#' 
#' x <- FormulaSums(z2, ~~region + kostragr * hovedint - 1)
#' z <- t(x) %*% z2$ant
#' 
#' # Run ReduceBy0
#' a <- ReduceBy0(x, z)
#' 
#' # The rows known to be 0 and omitted from a$x
#' z2[a$yKnown, ]
#' 
#' # Dimension of x reduced
#' dim(x)
#' dim(a$x)
#' 
#' # Additional elements know to be 0
#' b <- ReduceBy0(x, z, c(0, 0, rep(1, 42)))
#' z2[b$yKnown, ]
ReduceBy0 <- function(x, z, yStart = NULL) {
  z0 <- as.vector(as.matrix(z)) == 0
  y0 <- rowSums(x[, z0, drop = FALSE]) > 0
  
  if (!is.null(yStart)) 
    y0 <- y0 | (yStart == 0)
  x <- x[!y0, !z0, drop = FALSE]
  return(list(x = x, z = z[!z0, , drop = FALSE], yKnown = y0))
}

