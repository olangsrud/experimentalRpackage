#'Iterative proportional fitting from matrix input
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
#' @param eps stopping criterion. Maximum allowed value of \code{max(abs(z - t(x) \%*\% yHat))}
#' @param tol Another stopping criterion. Maximum absolute difference between two iteration. 
#' @param reduceBy0 When TRUE, \code{\link{ReduceBy0}}  used within the function 
#' @param reduceX When TRUE, \code{\link{ReduceBy0}} and  \code{\link{ReduceXspes}} used within the function (iteratively)
#' @param returnDetails More output when TRUE. Quite similar to \code{\link{ReduceX}}
#' @param y It is possible to set \code{z} to NULL and supply original \code{y} instead  (\code{z = t(x) \%*\% y})
#' @param ordSplit Whether to order the pieces of the split x according to a rule
#' @param altSplit Whether to try alternative splitting when ordinary not working properly
#'        
#'
#' @return \code{yHat}, the estimate of \code{y} 
#' 
#' @importFrom methods as
#' @importFrom utils flush.console
#' @importFrom Matrix drop0
#' 
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' data2 <- EasyData("z2")
#' x <- FormulaSums(data2, ~fylke + kostragr * hovedint - 1)
#' z <- t(x) %*% data2$ant  # same as FormulaSums(data2, ant~fylke + kostragr * hovedint -1)
#' yHat <- Mipf(x, z)
#' 
#' #############################
#' # loglin comparison 
#' #############################
#' 
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
#' # Transform the data for input to Mipf
#' df <- as.data.frame.table(tab)
#' names(df)[1:5] <- c("A", "B", "C", "D", "E")
#' x <- FormulaSums(df, ~A:B + A:C + A:D + A:E + B:C:D + C:D:E - 1)
#' z <- t(x) %*% df$Freq
#' 
#' # Estimate yHat by Mipf
#' yHatPMipf <- Mipf(x, z, iter = iter, eps = eps)
#' 
#' # Maximal absolute difference
#' max(abs(yHatPMipf - yHatLoglin))
#' 
#' # Note: loglin reports one iteration extra 
#' 
#' # Another example. Only one iteration needed.
#' max(abs(Mipf(x = FormulaSums(df, ~A:B + C - 1), 
#'              z = FormulaSums(df, Freq ~ A:B + C -1)) 
#'              - matrix(loglin(tab, list(1:2, 3), fit = TRUE)$fit, ncol = 1)))
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
#' # Reduction by 0, but no iteration improvement. Identical results.
#' t <- 360
#' y <- z3$ant
#' y[round((1:t) * 432/t)] <- 0
#' z <- t(x) %*% y
#' a1 <- Mipf(x, z, eps = 0.1)
#' a2 <- Mipf(x, z, reduceBy0 = TRUE, eps = 0.1)
#' a3 <- Mipf(x, z, reduceX = TRUE, eps = 0.1)
#' max(abs(a1 - a2))
#' max(abs(a1 - a3))
#' 
#' 
#' # Improvement by reduceX. Changing eps and iter give more similar results.
#' t <- 402
#' y <- z3$ant
#' y[round((1:t) * 432/t)] <- 0
#' z <- t(x) %*% y
#' a1 <- Mipf(x, z, eps = 1)
#' a2 <- Mipf(x, z, reduceBy0 = TRUE, eps = 1)
#' a3 <- Mipf(x, z, reduceX = TRUE, eps = 1)
#' max(abs(a1 - a2))
#' max(abs(a1 - a3))
#' 
#' # Example with small eps and "Iteration stopped since tol reached"
#' t <- 384
#' y <- z3$ant
#' y[round((1:t) * 432/t)] <- 0
#' z <- t(x) %*% y
#' a1 <- Mipf(x, z, eps = 1e-14)
#' a2 <- Mipf(x, z, reduceBy0 = TRUE, eps = 1e-14)
#' a3 <- Mipf(x, z, reduceX = TRUE, eps = 1e-14)
#' max(abs(a1 - a2))
#' max(abs(a1 - a3))
#' 
#' # All y-data found by reduceX (0 iterations). 
#' t <- 411
#' y <- z3$ant
#' y[round((1:t) * 432/t)] <- 0
#' z <- t(x) %*% y
#' a1 <- Mipf(x, z)
#' a2 <- Mipf(x, z, reduceBy0 = TRUE)
#' a3 <- Mipf(x, z, reduceX = TRUE)
#' max(abs(a1 - y))
#' max(abs(a2 - y))
#' max(abs(a3 - y))
Mipf <- function(x, z, iter = 100, yStart = matrix(1, nrow(x), 1), eps = 0.01, tol = 1e-10, 
                  reduceBy0 = FALSE, reduceX = FALSE, returnDetails = FALSE, y = NULL,
                  ordSplit = FALSE, altSplit = FALSE) {
  
  if (is.null(z))
    z <- Matrix::crossprod(x, y)
  
  if (reduceX) reduceBy0 <- TRUE
  
  if (reduceBy0) {
    cat("-")
    snx <- seq_len(nrow(x))
    flush.console()
    a <- ReduceBy0(x, z, yStart)
    aKnown <- as.vector(a$yKnown)
    yHat <- Matrix(0, length(aKnown), 1)
    yKnown <- aKnown
    if (any(aKnown)) {
      cat("0")
      flush.console()
    }
    rerun <- reduceX
    while (rerun) {
      rerun <- FALSE
      a <- ReduceXspes(a$x, a$z)
      aKnown <- as.vector(a$yKnown)
      if (any(aKnown)) {
        cat("x")
        flush.console()
        yHat[(snx[!yKnown])[seq_along(aKnown)[aKnown]], 1] <- a$y[seq_along(aKnown)[aKnown], 1]
        yKnown[!yKnown] <- aKnown
        cat("-")
        flush.console()
        a <- ReduceBy0(a$x, a$z)
        aKnown <- as.vector(a$yKnown)
        rerun <- TRUE
        if (any(aKnown)) {
          cat("0")
          flush.console()
          yHat[(snx[!yKnown])[aKnown], 1] <- 0
          yKnown[!yKnown] <- aKnown
        }
      }
    }
    
    cat("(",dim(x)[1],"*",dim(x)[2],"->", dim(a$x)[1],"*",dim(a$x)[2],")",sep="")
    
    if(any(!yKnown))
      yHat[seq_along(yKnown)[!yKnown], 1] <- Mipf(a$x, a$z, iter = iter, yStart = yStart[seq_along(yKnown)[!yKnown], 1], 
                                                   eps = eps, tol = tol, ordSplit = ordSplit, altSplit = altSplit)
    else
      cat("   0 iterations\n")
    
    deviation <- max(abs(crossprod(x, yHat) - z))
    
    
    if (!(deviation < eps)) 
      warning("Deviation limit exceeded")
    cat("Final deviation", deviation, "\n")
    if(returnDetails){
      return(list(x = a$x, z = a$z, yKnown = yKnown, y = yHat))
    }
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
  
  
  
  if (max(sapply(xL, function(x) max(rowSums(x)))) > 1) {
    if (!altSplit) {
      warning("Split into component not working properly. Try altSplit=TRUE ?")
    } else {
      warning("Alternative algorithm - since split into component not working properly")
      
      xL <- vector("list", 0)
      zL <- vector("list", 0)
      xU <- xT
      
      us <- UniqueSeq(xT@i, xT@j)
      xU@x <- as.numeric(us)
      
      
      m0 <- 0
      ma <- 0
      while (!is.na(ma)) {
        ma <- match(2, xU@x)
        if (is.na(ma)) {
          mEnd <- ncol(xU)
        } else {
          mEnd <- min(xU@j[xU@x == 2])
        }
        xL1 <- xT[, m0 + (1:mEnd), drop = FALSE]
        rS1 <- rowSums(xL1)
        if (max(rS1) > 1) 
          stop("something wrong")
        if (!is.na(ma)) {
          xU <- xU[, -(1:mEnd), drop = FALSE]
          w1 <- which(rS1 == 1)
          i1 <- xU@i %in% (w1 - 1)
          xU@x[i1] <- xU@x[i1] - 1
        }
        xL <- c(xL, xL1)
        zL <- c(zL, z[m0 + (1:mEnd), 1, drop = FALSE])
        m0 <- m0 + mEnd
      }
      n <- length(xL)
    }
  }
  
  
  if (ordSplit) {
    ColSumsMean <- function(x) {
      sum(x)/dim(x)[2]
    }
    ord <- order(sapply(xL, ColSumsMean), decreasing = TRUE)
  } else {
    ord <- seq_len(n)
  }
  
  
  # Run iterative proportional fitting 
  t <- 0
  deviation <- max(abs(crossprod(x, yStart) - z))
  deviationLast  <- Inf 
  k1 <- -1  # Used for printing progress
  while (t < iter & deviation > eps & abs(deviation - deviationLast) > tol) {
    t <- t + 1
    
    # Printing part
    k2 <- round(25 * sqrt(t/iter))
    if (k2 > k1) {
      cat(".")
      flush.console()
    }
    k1 <- k2
    
    # Computation part
    for (i in ord) {
      faktorZ <- zL[[i]]/crossprod(xL[[i]], yStart)
      faktorZ[is.na(faktorZ)] <- 1
      faktor <- rep(1, nrow(x))
      faktor[xL[[i]]@i + 1] <- faktorZ[xL[[i]]@j + 1]
      yStart <- faktor * yStart
    }
    deviationLast <- deviation 
    deviation <- max(abs(crossprod(x, yStart) - z))
    
  }
  if (!(deviation < eps)){
    if(!(abs(deviation - deviationLast) > tol)){
      warning("Iteration stopped since tol reached")  
    } else {
      warning("Iteration limit exceeded")
    }
  }
  cat("   ", t, "iterations: deviation", deviation, "\n")
  if(returnDetails){
    return(list(x = x, z = z, yKnown = rep(FALSE, dim(yStart)[1]), y = yStart))
  }
  yStart
}
















