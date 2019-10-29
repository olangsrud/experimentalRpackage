
#' Special weighted and ridge penalized regression 
#' 
#' By using \code{\link{nnls}}, \code{\link{ginv}} or \code{\link{glmnet}}
#'
#' @param x Input matrix, each row is an ordinary observation 
#' @param y Ordinary response observation (vector or matrix)
#' @param xExact Input matrix, each row is a highly weighted observation 
#' @param yExact Highly weighted response observation (vector or matrix)
#' @param wExact Weight for highly weighted observations
#' @param lambda Ridge regression penalty parameter (sequence when glmnet)
#'
#' @return Output from \code{\link{nnls}}, \code{\link{glmnet}} or coefficient estimate calculated using \code{\link{ginv}} 
#' @importFrom nnls nnls 
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' x <- cbind(1:11, -11:-1, c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2))
#' y <- matrix(5 - sin(2 * (1:11)), dimnames = list(NULL, "y"))
#' x1 <- x[1:9, ]
#' x2 <- x[10:11, ]
#' y1 <- y[1:9, , drop = FALSE]
#' y2 <- y[10:11, , drop = FALSE]
#' 
#' 
#' # Generalized inverse estimation
#' ginvCoef <- WRRginv(x1, y1, x2, y2)
#' 
#' # Non-negative estimation by nnls 
#' # Coefficients in output element x 
#' nnlsCoef <- WRRnnls(x1, y1, x2, y2)$x
#' 
#' # Non-negative estimation by glmnet
#' # Take out best fit from matrix of coefficients 
#' gn <- WRRglmnet(x1, y1, x2, y2)
#' glmnetCoef <- coef(gn)[-1, which.max(gn$dev.ratio), drop = FALSE]
#' 
#' # Another estimation by glmnet (not non-negative)
#' # Take out best fit from matrix of coefficients
#' gnInf <- WRRglmnet(x1, y1, x2, y2, lower.limits = -Inf)
#' glmnetCoefInf <- coef(gnInf)[-1, which.max(gn$dev.ratio), drop = FALSE]
#' 
#' # All coefficients
#' coef4 <- as.matrix(cbind(ginvCoef, nnlsCoef, glmnetCoef, glmnetCoefInf))
#' colnames(coef4) <- c("ginv", "nnls", "glmnet", "glmnetInf")
#' print(coef4)
#' 
#' # Original y and fitted values. Close fit for last two observation.
#' cbind(y, x %*% coef4)
WRRnnls <- function(x, y, xExact = NULL, yExact = NULL, wExact = 1000, lambda = 0.0001^2) {
  if (NROW(xExact) == 0) {
    xExact <- NULL
    yExact <- NULL
  }
  nnls(rbind(x, wExact * xExact, sqrt(lambda) * diag(ncol(x))), rbind(y, wExact * yExact, matrix(rep(0, ncol(x)))))
}

#' @rdname WRRnnls
#' @importFrom MASS ginv
#' @export
WRRginv <- function(x, y, xExact = NULL, yExact = NULL, wExact = 1000, lambda = 0.0001^2) {
  if (NROW(xExact) == 0) {
    return(ginv(as.matrix(x)) %*% y)
  }
  ginv(as.matrix(rbind(x, wExact * xExact))) %*% rbind(y, wExact * yExact)
}


#' @rdname WRRnnls
#' @param intercept glmnet parameter
#' @param standardize glmnet parameter
#' @param thresh glmnet parameter
#' @param lower.limits glmnet parameter
#' @param ...  Further glmnet parameters
#' @importFrom glmnet glmnet
#' @export
WRRglmnet <- function(x, y, xExact = NULL, yExact = NULL, wExact = 1000, lambda = exp(((2 * 12):(-17 * 2))/2), 
                      intercept = FALSE, standardize = FALSE, thresh = 1e-10, lower.limits = 0, ...) {
  n <- nrow(x)
  if (is.null(xExact)) 
    nExact <- 0 else nExact <- nrow(xExact)
    gn <- glmnet(rbind(x, xExact), rbind(y, yExact), lambda = lambda, 
                 weights = c(rep(1, n), rep(wExact, nExact)), intercept = intercept, 
                 standardize = standardize, alpha = 0, family = "gaussian", 
                 lower.limits = lower.limits, thresh = thresh, ...)
}

