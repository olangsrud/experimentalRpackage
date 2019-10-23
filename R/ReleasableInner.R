#' Releasable inner cell frequencies by post-processing protected tabular data
#' 
#' Margins are first perturbed according to pTable or differential privacy (Laplace noise) by
#' \code{\link{NoisyCoverMargins}}.
#' Thereafter inner cells are estimated from the perturbed margins by
#' \code{\link{WRRnnls}}, \code{\link{WRRginv}} or \code{\link{WRRglmnet}}.
#'
#' @param data Input data as a data frame (inner cells)
#' @param freqVar Variable holding counts (name or number)
#' @param formula Model formula defining cells (margins) to be perturbed 
#' @param formulaExact Model formula defining cells (margins) to be exact
#' @param regMethod One of "nnls" (default),"ginv" or"glmnet". See \code{\link{WRRnnls}}.
#' @param returnDetails Whether to include all output from the underlying functions.
#' @param rerun Possibility for a variant that iterates. Ordinary margins are re-computed in each round while the exact margins are unchanged.
#' @param eps Differential privacy parameter
#' @param pMatrix Output from \code{\link{Pmatrix}}
#' @param keyVar Variable holding uniformly distributed keys (name or number).
#' @param wExact Weight for highly weighted margins 
#' @param intercept  Parameter to \code{\link{CoverMarginsExact}}
#' @param interceptExact Parameter to \code{\link{CoverMarginsExact}}
#' @param coverMargins Parameter to \code{\link{CoverMarginsExact}}
#' @param ... Parameters to or \code{\link{WRRnnls}}, \code{\link{WRRginv}} or \code{\link{WRRglmnet}}
#'
#' @return Estimated inner cells frequencies (yHat) or a list including all output from
#'         \code{\link{NoisyCoverMargins}} and \code{\link{WRRnnls}}, \code{\link{WRRginv}} or \code{\link{WRRglmnet}}.
#'         
#' @seealso  Another way to obtain releasable inner cell frequencies from protected tabular data is to 
#'           run \code{\link{SuppressDec}} with input from \code{\link{PTxyz}}.        
#' 
#' @export
#'
#' @examples
#' y <- c(3, 1, 5, 6, 4, 8, 2, 7, 27)
#' cols <- paste("col", col(matrix(1:9, 3, 3)), sep = "")
#' rows <- paste("row", row(matrix(1:9, 3, 3)), sep = "")
#' z <- data.frame(rows, cols, y)
#' pMatrix <- Pmatrix()
#' z$keys <- runif(9)
#' 
#' # Non-negative estimates without noise
#' ReleasableInner(z, "y", ~rows + cols, eps = Inf)
#' 
#' # Non-negative estimates with Laplace noise
#' ReleasableInner(z, "y", ~rows + cols, eps = 0.5)
#' 
#' # Using pTable and keys. Overall total perturbed.
#' ReleasableInner(z, "y", ~rows + cols, pMatrix = pMatrix, keyVar = "keys", intercept = TRUE)
#' 
#' # Generalized inverse (not non-negative)
#' ReleasableInner(z, "y", ~rows + cols, pMatrix = pMatrix, keyVar = "keys", intercept = TRUE, 
#'                 regMethod = "ginv")
#' 
#' # Preserving row totals
#' ReleasableInner(z, "y", ~rows + cols, ~rows)
#' 
#' # Other data and detailed output
#' z2 <- EasyData("z2")
#' a1 <- ReleasableInner(z2, "ant", ~region + kostragr * hovedint, ~kostragr + hovedint, 
#'                       returnDetails = TRUE)
#' 
#' # Other values of thresh and wExact
#' a2 <- ReleasableInner(z2, "ant", ~region + kostragr * hovedint, ~kostragr + hovedint, 
#'                       regMethod = "glmnet", thresh = 1E-8, wExact = 100)
#' 
#' # Not non-negative
#' a3 <- ReleasableInner(z2, "ant", ~region + kostragr * hovedint, ~kostragr + hovedint, 
#'                       regMethod = "glmnet", lower.limits = -Inf)
ReleasableInner <- function(data, freqVar, formula, formulaExact = NULL, regMethod = c("nnls", "ginv", "glmnet"), 
                            returnDetails = FALSE, rerun = 0, eps = 0.5, pMatrix = NULL, keyVar = NULL, 
                            wExact = 1000, intercept = FALSE, interceptExact = !is.null(formulaExact), 
                            coverMargins = TRUE, ...) {
  regMethod <- match.arg(regMethod)
  obj <- NoisyCoverMargins(data = data, freqVar = freqVar, formula = formula, formulaExact = formulaExact, 
                           eps = eps, pMatrix = pMatrix, keyVar = keyVar, 
                           intercept = intercept, interceptExact = interceptExact)
  NoisyCoverMarginsToInner(obj, regMethod = regMethod, returnDetails = returnDetails, rerun = rerun, ...)
}
    
   
#' @rdname ReleasableInner
#' @param obj Output from \code{\link{NoisyCoverMargins}}
#' @importFrom Matrix t
#' @importFrom stats coef
#' @export
NoisyCoverMarginsToInner <- function(obj, regMethod = c("nnls", "ginv", "glmnet"), returnDetails = FALSE, rerun = 0, ...) {
  regMethod <- match.arg(regMethod)
  
  if (regMethod == "nnls") 
    WRR <- WRRnnls
  if (regMethod == "ginv") 
    WRR <- WRRginv
  if (regMethod == "glmnet") 
    WRR <- WRRglmnet
  
  for (i in 0:rerun) {
    if (i == 0) 
      wrr <- WRR(t(obj$x), obj$zPerturbed, t(obj$xExact), obj$zExact, ...) 
    else 
      wrr <- WRR(t(obj$x), crossprod(obj$x, yHat), t(obj$xExact), obj$zExact, ...)
    if (regMethod == "nnls") 
      yHat <- matrix(wrr$x, ncol = 1)
    if (regMethod == "ginv") 
      yHat <- wrr
    if (regMethod == "glmnet") 
      yHat <- coef(wrr)[-1, which.max(wrr$dev.ratio), drop = FALSE]
  }
  
  if (!returnDetails) 
    return(yHat)
  
  list(yHat = yHat, noisyCoverMargins = obj, WRRoutput = wrr)
}
