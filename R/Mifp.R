#' Iterative proportional fitting from matrix input
#' 
#' The linear equation, \code{z = t(x) \%*\% y}, is (hopefully)  solved for \code{y} by
#' iterative proportional fitting
#' 
#' The algorithm will work similar to \code{\link{loglin}} when the input x-matrix is a overparameterized model matrix 
#' – as can be created by \code{\link{FormulaSums}}. See Examples.
#'
#' @param x a matrix 
#' @param z a single column matrix
#' @param iter maximum number of iterations
#' @param yStart a starting estimate of \code{y}
#' @param eps maximum allowed value of \code{max(abs(z - t(x) \%*\% yHat))} 
#' @param reduceBy0 When TRUE, \code{\link{ReduceBy0}}  used within the function 
#' @param reduceX When TRUE, \code{\link{ReduceBy0}} and  \code{\link{ReduceX}} used within the function (iteratively)
#'        
#'
#' @return \code{yHat}, the estimate of \code{y} 
#' 
#' @importFrom methods as
#' @importFrom utils flush.console
#' @importFrom Matrix drop0
#' 
#' @export
#'
#' @examples
#' # Generate input data for loglin
#' n <- 5:9
#' tab <- array(sample(1:prod(n)), n)
#' 
#' # Input parameters
#' iter <- 20
#' eps <- 1e-05
#' 
#' # Estimate yHat by loglin
#' out <- loglin(tab, list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3, 4), c(3, 4, 5)), 
#'               fit = TRUE, iter = iter, eps = eps)
#' yHatLoglin <- matrix(((out$fit)), ncol = 1)
#' 
#' # Transform the data for input to Pfifp
#' df <- as.data.frame.table(tab)
#' names(df)[1:5] <- c("A", "B", "C", "D", "E")
#' x <- FormulaSums(df, ~A:B + A:C + A:D + A:E + B:C:D + C:D:E - 1)
#' z <- t(x) %*% df$Freq
#' 
#' # Estimate yHat by Pfifp
#' yHatPfifp <- Mifp(x, z, iter = iter, eps = eps)
#' 
#' # Maximal absolute difference
#' max(abs(yHatPfifp - yHatLoglin))
#' 
#' 
#' #########################################
#' #  reduceBy0  and  reduceX examples 
#' #########################################
#' 
#' #' z3 <- EasyData("z3")
#' x <- FormulaSums(z3, ~region + kostragr * hovedint + region * mnd2 + fylke * mnd + 
#'                      mnd * hovedint + mnd2 * fylke * hovedint - 1)
#' 
#' 
#' t <- 360
#' y <- z3$ant
#' y[round((1:t) * 432/t)] <- 0
#' z <- t(x) %*% y
#' a1 <- Mifp(x, z, eps = 0.1)
#' a2 <- Mifp(x, z, reduceBy0 = TRUE, eps = 0.1)
#' a3 <- Mifp(x, z, reduceX = TRUE, eps = 0.1)
#' 
#' max(abs(a1 - a2))
#' max(abs(a1 - a3))
#' 
#' 
#' t <- 402
#' y <- z3$ant
#' y[round((1:t) * 432/t)] <- 0
#' z <- t(x) %*% y
#' a1 <- Mifp(x, z, eps = 1)
#' a2 <- Mifp(x, z, reduceBy0 = TRUE, eps = 1, iter = iter)
#' a3 <- Mifp(x, z, reduceX = TRUE, eps = 1, iter = iter)
#' max(abs(a1 - a2))
#' max(abs(a1 - a3))
#' 
#' 
#' t <- 411
#' y <- z3$ant
#' y[round((1:t) * 432/t)] <- 0
#' z <- t(x) %*% y
#' a1 <- Mifp(x, z)
#' a2 <- Mifp(x, z, reduceBy0 = TRUE)
#' a3 <- Mifp(x, z, reduceX = TRUE)
#' max(abs(a1 - a2))
#' max(abs(a1 - a3))
Mifp <- function(x, z, iter = 100, yStart = matrix(1, nrow(x), 1), eps = 0.01, 
                 reduceBy0 = FALSE, reduceX = FALSE) {
  
  if (reduceX) reduceBy0 <- TRUE
  
  if (reduceBy0) {
    cat("0")
    snx <- seq_len(nrow(x))
    flush.console()
    a <- ReduceBy0(x, z, yStart)
    aKnown <- as.vector(a$yKnown)
    yHat <- Matrix(0, length(aKnown), 1)
    yKnown <- aKnown
    rerun <- reduceX
    while (rerun) {
      rerun <- FALSE
      cat("x")
      flush.console()
      a <- ReduceXspes(a$x, a$z)
      aKnown <- as.vector(a$yKnown)
      if (any(aKnown)) {
        yHat[(snx[!yKnown])[seq_along(aKnown)[aKnown]], 1] <- a$y[seq_along(aKnown)[aKnown], 1]
        yKnown[!yKnown] <- aKnown
        cat("0")
        flush.console()
        a <- ReduceBy0(a$x, a$z)
        aKnown <- as.vector(a$yKnown)
        if (any(aKnown)) {
          yHat[(snx[!yKnown])[aKnown], 1] <- 0
          yKnown[!yKnown] <- aKnown
          rerun <- TRUE
        }
      }
    }
    
    cat("(",dim(x)[1],"*",dim(x)[2],"->", dim(a$x)[1],"*",dim(a$x)[2],")",sep="")
    
    if(any(!yKnown))
      yHat[seq_along(yKnown)[!yKnown], 1] <- Mifp(a$x, a$z, iter = iter, yStart = yStart[seq_along(yKnown)[!yKnown], 1], eps = eps)
    else
      cat("   0 iterations\n")
    
    deviation <- max(abs(crossprod(x, yHat) - z))
    
  
    if (!(deviation < eps)) 
      warning("Deviation limit exceeded")
    cat("Final deviation", deviation, "\n")
    return(yHat)
  }

  
  cat(":")
  flush.console()
  
  # Split input into components
  xT <- as(drop0(x), "dgTMatrix")
  cat("-")
  flush.console()
  us <- UniqueSeq(xT@i, xT@j)
  cat("-")
  n <- max(us)
  ma <- match(seq_len(n), us)
  cat("-")
  startCol <- c(xT@j[ma] + 1, ncol(x) + 1)
  xL <- vector("list", n)
  zL <- vector("list", n)
  cat("-")
  flush.console()
  for (i in seq_len(n)) {
    xL[[i]] <- xT[, startCol[i]:(startCol[i + 1] - 1), drop = FALSE]
    zL[[i]] <- z[startCol[i]:(startCol[i + 1] - 1), 1, drop = FALSE]
  }
  
  # Run iterative proportional fitting 
  t <- 0
  deviation <- max(abs(crossprod(x, yStart) - z))
  k1 <- -1  # Used for printing progress
  while (t < iter & deviation > eps) {
    t <- t + 1
    
    # Printing part
    k2 <- round(25 * sqrt(t/iter))
    if (k2 > k1) {
      cat(".")
      flush.console()
    }
    k1 <- k2
    
    # Computation part
    for (i in seq_len(n)) {
      faktorZ <- zL[[i]]/crossprod(xL[[i]], yStart)
      faktorZ[is.na(faktorZ)] <- 1
      faktor <- rep(1, nrow(x))
      faktor[xL[[i]]@i + 1] <- faktorZ[xL[[i]]@j + 1]
      yStart <- faktor * yStart
      # faktor2 = faktor*faktor2
    }
    deviation <- max(abs(crossprod(x, yStart) - z))
    
  }
  if (!(deviation < eps)) 
    warning("Iteration limit exceeded")
  cat("   ", t, "iterations: deviation", deviation, "\n")
  yStart
}

















