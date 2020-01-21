#' Join (merge) two data frames 
#'
#' @param x data frame
#' @param y data frame
#' @param by column names 
#' @param xLast Whether to reuse rows of x
#' @param yLast Whether to reuse rows of y
#' @param xAll Whether to match all rows of x by adding NAs 
#' @param yAll Whether to match all rows of y by adding NAs 
#' @param doSort  Whether to sort on the by columns
#'
#' @return data frame
#' @export
#' @examples 
#' x <- data.frame(ABC = c("A", "A", "A", "B", "B", "B", "B", "C", "C"), 
#'                 OneTwo = c(1, 1, 1, 1, 1, 2, 1, 1, 1), 
#'                 x = 10 * (1:9), stringsAsFactors = FALSE)
#' y <- data.frame(ABC = c("A", "A", "A", "B", "C", "C"), 
#'                 OneTwo = c(1, 1, 2, 1, 1, 2), 
#'                 y = 101:106, stringsAsFactors = FALSE)
#' x
#' y                 
#' 
#' # Inner join
#' SasJoin(x, y, xLast = FALSE, yLast = FALSE, xAll = FALSE, yAll = FALSE)
#' 
#' # Left join
#' SasJoin(x, y, xLast = FALSE, yLast = FALSE, xAll = TRUE, yAll = FALSE)
#' 
#' # Left join by reusing rows of y when possible
#' SasJoin(x, y, xLast = FALSE, yLast = TRUE, xAll = TRUE, yAll = FALSE)
#' 
#' # Reusing rows of y when possible but not a full left join
#' SasJoin(x, y, xLast = FALSE, yLast = TRUE, xAll = FALSE, yAll = FALSE)
#' 
#' # Reusing rows of both x and y
#' SasJoin(x, y, xLast = TRUE, yLast = TRUE, xAll = FALSE, yAll = FALSE)
#' 
#' # Outer join by reusing rows of both x and y
#' SasJoin(x, y, xLast = TRUE, yLast = TRUE, xAll = TRUE, yAll = TRUE)
#' 
#' # Outer join without reusing rows
#' SasJoin(x, y, xLast = FALSE, yLast = FALSE, xAll = TRUE, yAll = TRUE)
#'
SasJoin <- function(x, y, by = intersect(names(x), names(y)), xLast = TRUE, yLast = TRUE, xAll = TRUE, yAll = TRUE, doSort = TRUE) {
  if (class(x)[1] != "data.frame") 
    x <- as.data.frame(x)
  if (class(y)[1] != "data.frame") 
    y <- as.data.frame(y)
  
  if (length(by) == 0) {
    rg <- rep(1L, nrow(x) + nrow(y))
    doSort <- FALSE
  } else {
    rg <- SSBtools::RowGroups(rbind(x[, by, drop = FALSE], y[, by, drop = FALSE]))
  }
  
  xNr <- seq_len(nrow(x))
  yNr <- seq_len(nrow(y))
  
  xRg <- rg[xNr]
  yRg <- rg[-xNr]
  
  xRgSort <- sort(xRg, index.return = TRUE)
  yRgSort <- sort(yRg, index.return = TRUE)
  
  xSnr <- (xNr - xNr[match(xRgSort$x, xRgSort$x)])[order(xRgSort$ix)]
  ySnr <- (yNr - yNr[match(yRgSort$x, yRgSort$x)])[order(yRgSort$ix)]
  
  RevMatch <- function(x, y) (length(y):1)[match(x, rev(y))]
  
  if (xLast) 
    xL <- RevMatch(xRg[xSnr == 0], xRg) else xL <- integer(0)
  
  if (yLast) 
    yL <- RevMatch(yRg[ySnr == 0], yRg) else yL <- integer(0)
  
  ma <- SSBtools::Match(data.frame(a = xRg, b = xSnr), data.frame(a = yRg, b = ySnr))
  
  isnama <- is.na(ma)
  
  if (xAll | yLast) {
    ma1 <- match(xRg[isnama], yRg[yL])
    indX <- xNr
    indY <- ma
    indY[isnama] <- yL[ma1]
    if (!xAll) {
      xOK <- rep(TRUE, length(xNr))
      xOK[isnama][is.na(ma1)] <- FALSE
      indX <- indX[xOK]
      indY <- indY[xOK]
    }
  } else {
    indX <- xNr[!isnama]
    indY <- ma[!isnama]
  }
  
  z <- cbind(x[indX, , drop = FALSE], y[indY, setdiff(names(y), by), drop = FALSE])
  
  if (yAll | xLast) {
    ma2 <- match(yRg[-ma[!isnama]], xRg[xL])
    indX2 <- xL[ma2]
    indY2 <- yNr[-ma[!isnama]]
    if (!yAll) {
      indX2 <- indX2[!is.na(ma2)]
      indY2 <- indY2[!is.na(ma2)]
    }
    z <- rbind(z, cbind(x[indX2, setdiff(names(x), by), drop = FALSE], y[indY2, , drop = FALSE]))
  }
  
  if (doSort) {
    z <- z[SSBtools::SortRows(z[, by, drop = FALSE], index.return = TRUE), , drop = FALSE]
  }
  
  rownames(z) <- NULL
  z
}

