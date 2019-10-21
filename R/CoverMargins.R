#' Extract a set of model terms that covers a formula 
#' 
#' Other model terms depend on the extracted set. Hierarchical relationships seen in the data is considered. 
#'
#' @param data Input data as a data frame (inner cells)
#' @param formula Model formula defining cells to consider 
#' @param coverTable Whether to return a model term defining a cover table 
#'
#' @return Vector of extracted model terms or a two-component list when coverTable is TRUE
#' 
#' @note all terms retuned when coverTable=NA
#' 
#' @importFrom SSBtools HierarchicalGroups SeqInc
#' @importFrom stats rexp terms.formula
#' @export
#'
#' @examples
#' 
#' z3 <- EasyData("z3")
#' CoverMargins(z3, ~region + kostragr * hovedint)
#' CoverMargins(z3, ~region * mnd + hovedint * mnd + fylke * hovedint * mnd + 
#'                   kostragr * hovedint * mnd)
#' CoverMargins(z3, ~fylke * hovedint + kostragr * hovedint, TRUE)
#' CoverMargins(z3, ~hovedint * mnd2 + mnd + hovedint + kostragr * hovedint, TRUE)
CoverMargins <- function(data, formula, coverTable = FALSE) {
  x <- attr(terms.formula(formula), "factors")
  if (is.na(coverTable)) {
    return(colnames(x))
  }
  hg <- HierarchicalGroups(data[, all.vars(formula), drop = FALSE], eachName = TRUE)
  x[x > 0] <- 1L
  for (i in seq_along(hg)) {
    for (j in SeqInc(2, length(hg[[i]]))) {
      x[rownames(x) == hg[[i]][j], x[rownames(x) == hg[[i]][j - 1], ] == 1] <- 1
    }
  }
  aTerms <- Aterms(x)
  if (!coverTable) 
    return(aTerms)
  z <- character(0)
  for (i in seq_along(hg)) {
    hg[[i]] <- hg[[i]][hg[[i]] %in% rownames(x)]
  }
  hg <- hg[lapply(hg, length) > 0]
  hg <- lapply(hg, function(x) x[1])
  hg <- unlist(hg)
  dup <- names(hg) %in% names(hg)[duplicated(names(hg))]  # use main var when problem
  hg[dup] <- names(hg)[dup]
  hg <- unique(hg)
  list(coverMargins = aTerms, coverTable = paste(hg, collapse = ":"))
}


Aterms <- function(x) { # x = attr(terms.formula(~a*B*C + c*D),'factors')
  x <- x[, order(colSums(x), decreasing = TRUE), drop = FALSE]
  for (i in seq_len(ncol(x))) {
    z <- x[!as.logical(x[, i]), , drop = FALSE]
    cc <- colSums(z) > 0
    cc[seq_len(i)] <- TRUE
    x <- x[, cc, drop = FALSE]
    if (sum(cc) == i) 
      return(colnames(x))
  }
}



#' Extract two sets of model terms that cover two formulas 
#' 
#' One formula for cells to be perturbed and one formula for cells to be exact
#' 
#' The function make use of \code{\link{CoverMargins}}. Thus, hierarchical relationships seen in the data is considered. 
#' Each output element can be input to \code{\link{Terms2formula}}.
#'
#' @param data Input data as a data frame (inner cells)
#' @param formula Model formula defining cells to be perturbed 
#' @param formulaExact Model formula defining cells to be exact
#' @param intercept Whether to include intercept
#' @param interceptExact Whether to include intercept
#' @param coverMargins Whether to extract only the covering terms
#'
#' @return List of two sets of model terms
#' 
#' @importFrom CalibrateSSB MatchVarNames
#' @export
#'
#' @examples
#' z3 <- EasyData("z3")
#' CoverMarginsExact(z3, ~region + kostragr * hovedint)
#' CoverMarginsExact(z3, ~region + kostragr * hovedint * mnd, ~kostragr * hovedint)
#' CoverMarginsExact(z3, ~region + kostragr * hovedint * mnd, ~kostragr * hovedint + mnd)
#' CoverMarginsExact(z3, ~region, ~kostragr * hovedint, TRUE, FALSE)
#' CoverMarginsExact(z3, ~region * mnd, ~kostragr * hovedint, coverMargins = FALSE)
#' CoverMarginsExact(z3, ~region, interceptExact = TRUE)  
CoverMarginsExact <- function(data, formula, formulaExact = NULL, 
                              intercept = FALSE, interceptExact = !is.null(formulaExact), coverMargins = TRUE) {
  if (coverMargins) 
    coverTable <- FALSE else coverTable <- NA
    c1 <- CoverMargins(data, formula, coverTable)
    if (intercept) 
      c1 <- c("(Intercept)", c1)
    if (!is.null(formulaExact) | interceptExact) {
      if (!is.null(formulaExact)) 
        c0 <- CoverMargins(data, formulaExact, coverTable) else c0 <- NULL
        if (interceptExact) 
          c0 <- c("(Intercept)", c0)
        c1 <- c1[is.na(CalibrateSSB::MatchVarNames(c1, c0))]
    } else {
      c0 <- NULL
    }
    list(coverMargins = c1, coverMarginsExact = c0)
}




#' Laplace or pTable perturbed margins
#'
#' The function make use of \code{\link{CoverMarginsExact}}.
#' 
#' @param data Input data as a data frame (inner cells)
#' @param freqVar Variable holding counts (name or number)
#' @param formula Model formula defining cells to be perturbed 
#' @param formulaExact Model formula defining cells to be exact
#' @param eps Differential privacy parameter
#' @param pMatrix Output from \code{\link{Pmatrix}}
#' @param keyVar Variable holding uniformly distributed keys (name or number)
#' @param ... Further parameters to \code{\link{CoverMarginsExact}}
#'
#' @return A list 
#' 
#' @importFrom SSBtools FormulaSums
#' @importFrom Matrix Matrix crossprod
#' @importFrom stats runif
#' @export
#'
#' @examples
#' z2 <- EasyData("z2")
#' a1 <- NoisyCoverMargins(z2, "ant", ~region + kostragr * hovedint)
#' a2 <- NoisyCoverMargins(z2, "ant", ~region + kostragr * hovedint, ~kostragr + hovedint, eps = 0.25)
#' a3 <- NoisyCoverMargins(cbind(z2,keyVar = runif(44) * (z2$ant>0)), "ant", 
#'                         ~region + kostragr * hovedint, ~kostragr + hovedint, 
#'                         pMatrix = Pmatrix(), keyVar = "keyVar")
NoisyCoverMargins <- function(data, freqVar, formula, formulaExact = NULL, eps = 0.5, pMatrix = NULL, keyVar = NULL, ...) {
  cme <- CoverMarginsExact(data, formula, formulaExact, ...)
  f <- Terms2formula(cme$coverMargins)
  x <- FormulaSums(data, f)
  if (!is.null(cme$coverMarginsExact)) {
    fExact <- Terms2formula(cme$coverMarginsExact)
    xExact <- FormulaSums(data, fExact)
  } else {
    xExact <- x[, integer(0), drop = FALSE]
  }
  y <- Matrix(data[, freqVar], ncol = 1)
  
  z <- crossprod(x, y)
  zExact <- crossprod(xExact, y)
  
  nMargins <- length(cme$coverMargins)
  
  if (is.null(pMatrix)) {
    zPerturbed <- z + Lap(nrow(z), nMargins/eps)
  } else {
    if (!is.null(keyVar)) 
      rKey <- crossprod(x, data[, keyVar])[, 1, drop = TRUE]%%1 else {
        yNon0 <- y != 0
        y[yNon0] <- runif(sum(yNon0))
        rKey <- crossprod(x, y)[, 1, drop = TRUE]%%1
      }
    zPerturbed <- Pconvert(z, pMatrix, rKey)
  }
  list(coverMargins = cme$coverMargins, coverMarginsExact = cme$coverMarginsExact, nMargins = nMargins, 
       x = x, xExact = xExact, z = z, zExact = zExact, zPerturbed = zPerturbed)
}





