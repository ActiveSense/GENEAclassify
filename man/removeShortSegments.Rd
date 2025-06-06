% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/removeShortSegments.R
\name{removeShortSegments}
\alias{removeShortSegments}
\title{Remove Short Segments from Data}
\usage{
removeShortSegments(shortduration, changes, variance, time)
}
\arguments{
\item{shortduration}{integer vector, the elements that are determined to have short durations}

\item{changes}{integer vector, the estimated change points}

\item{variance}{numeric vector, the variance of the parameter estimate for the cpt object}

\item{time}{numeric vector, the time vector}
}
\value{
list with elements:\enumerate{
    \item changes, the new change points, 
    \item times, the new times at the start and end of each segment, and 
    \item duration, the new durations of each segment
}
}
\description{
Remove the short segments from the data based on the estimated 
variance in \code{\link{changeTimes}}. If variance is small, 
segments specified may be retained.
}
\examples{
library(changepoint)
set.seed(45265)
tm0 <- 1001:1060
d0 <- round(c(cumsum(runif(n = 20) * 2), 
    20:1 + rnorm(n = 20), 
    runif(n = 20) * 10))  
## identify changes in variance
c0 <- cpt.var(d0, penalty = "SIC", pen.value = 1e-3, method = "PELT")
## times of changepoints
cp0 <- cpts(c0)
t0 <- c(tm0[1], tm0[sort(cp0)], tm0[60])
sdur0 <- which(diff(t0) < 10)
v0 <- param.est(c0)$variance
## Note that variance for early changepoints is low
## so they are not removed.
GENEAclassify:::removeShortSegments(shortduration = sdur0, 
    changes = cp0, 
    variance = v0, 
    time = tm0)
}
\keyword{internal}
