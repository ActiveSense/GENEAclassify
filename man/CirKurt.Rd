% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{CirKurt}
\alias{CirKurt}
\title{CirKurt function}
\usage{
CirKurt(rotation)
}
\arguments{
\item{rotation}{vector of numeric values}

}
\value{
A single value data.
}
\description{
Called by \code{segmentation}.
}
\examples{
   tmp1 <- c(1,3,2,6,4,5,3,9,10)
   CirKurt(tmp1)
}
\keyword{internal}
