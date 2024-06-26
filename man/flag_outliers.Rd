% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flag_outliers.R
\name{flag_outliers}
\alias{flag_outliers}
\title{Flag outliers based on quantile ranges}
\usage{
flag_outliers(x, qlow = 0, qhigh = 1, m, keep, exclude = NULL)
}
\arguments{
\item{x}{numeric; a vector of values}

\item{qlow}{numeric; lower quantile of a quantile range}

\item{qhigh}{numeric; higher quantile for a quantile range}

\item{m}{numeric; range multiplier}

\item{keep}{character; the value range you would like to keep. Can be
\code{"higher"}, \code{"lower"} or \code{"middle"}.}

\item{exclude}{numeric; values to exclude from the calculation}
}
\description{
This function flags outliers based on quantile ranges. It is a wrapper around
the \code{\link{quantile}} function. It can be used for flagging outliers
based on the interquartile range (IQR) or any other quantile range.
}
\examples{
set.seed(0)
data <- rnorm(100, 0, 1)
# Calculate interquartile range (IQR), keep middle range
qc(data, qlow = 0.25, qhigh = 0.75, m = 1.5, keep = "middle")
# Keep high values (flag values well below the 0.05 quantile)
qc(data, qlow = 0.05, qhigh = 1, m = 1, keep = "higher")
# Keep low value (flag values well above the 0.95 quantile)
qc(data, qlow = 0, qhigh = 0.95, m = 1, keep = "lower")
# Keep low value, but exclude values of 0
qc(data, qlow = 0, qhigh = 0.95, m = 1, keep = "lower", exclude = 0)
}
