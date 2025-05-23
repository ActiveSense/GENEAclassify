% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{impact}
\alias{impact}
\title{impact function}
\usage{
impact(x, threshold = 2.5)
}
\arguments{
\item{x}{vector of numeric values}

\item{threshold}{cutoff point above which acceleration is considered}
}
\value{
A single value data.
}
\description{
Called by \code{segmentation}.
}
\examples{
   Freq <- 100
   x <- c(1,1.5,2.6,2,3,3,1,1,1)
   impact(x)
}
\keyword{internal}
