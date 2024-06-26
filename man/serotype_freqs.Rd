% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serotype_freqs.R
\name{serotype_freqs}
\alias{serotype_freqs}
\title{Calculate serotype frequences}
\usage{
serotype_freqs(df, group_by = NULL)
}
\arguments{
\item{df}{data.frame, a table with one row for each sample. The table must
have a `serotype` variable and an optional grouping variable.}

\item{group_by}{character; the name of the grouping variable. If not `NULL`
the input data frame must contain a variable with this name.}
}
\value{
A data frame.
}
\description{
This function takes a data frame with a single grouping variable, e.g.
`region23` or `country` and prepares a table in which for each level of the
grouping variable serotypes are listed in decreasing frequence order.
}
\examples{
serotype_freqs(df, group_by = "region23")
}
