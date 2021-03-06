% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SasJoin.R
\name{SasJoin}
\alias{SasJoin}
\title{Join (merge) two data frames}
\usage{
SasJoin(
  x,
  y,
  by = intersect(names(x), names(y)),
  xLast = TRUE,
  yLast = TRUE,
  xAll = TRUE,
  yAll = TRUE,
  doSort = TRUE
)
}
\arguments{
\item{x}{data frame}

\item{y}{data frame}

\item{by}{column names}

\item{xLast}{Whether to reuse rows of x}

\item{yLast}{Whether to reuse rows of y}

\item{xAll}{Whether to match all rows of x by adding NAs}

\item{yAll}{Whether to match all rows of y by adding NAs}

\item{doSort}{Whether to sort on the by columns}
}
\value{
data frame
}
\description{
Join (merge) two data frames
}
\examples{
x <- data.frame(ABC = c("A", "A", "A", "B", "B", "B", "B", "C", "C"), 
                OneTwo = c(1, 1, 1, 1, 1, 2, 1, 1, 1), 
                x = 10 * (1:9), stringsAsFactors = FALSE)
y <- data.frame(ABC = c("A", "A", "A", "B", "C", "C"), 
                OneTwo = c(1, 1, 2, 1, 1, 2), 
                y = 101:106, stringsAsFactors = FALSE)
x
y                 

# Inner join
SasJoin(x, y, xLast = FALSE, yLast = FALSE, xAll = FALSE, yAll = FALSE)

# Left join
SasJoin(x, y, xLast = FALSE, yLast = FALSE, xAll = TRUE, yAll = FALSE)

# Left join by reusing rows of y when possible
SasJoin(x, y, xLast = FALSE, yLast = TRUE, xAll = TRUE, yAll = FALSE)

# Reusing rows of y when possible but not a full left join
SasJoin(x, y, xLast = FALSE, yLast = TRUE, xAll = FALSE, yAll = FALSE)

# Reusing rows of both x and y
SasJoin(x, y, xLast = TRUE, yLast = TRUE, xAll = FALSE, yAll = FALSE)

# Outer join by reusing rows of both x and y
SasJoin(x, y, xLast = TRUE, yLast = TRUE, xAll = TRUE, yAll = TRUE)

# Outer join without reusing rows
SasJoin(x, y, xLast = FALSE, yLast = FALSE, xAll = TRUE, yAll = TRUE)

}
