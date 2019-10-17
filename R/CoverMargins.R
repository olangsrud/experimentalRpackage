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

