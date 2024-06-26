% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global_prevalent.R
\name{global_prevalent}
\alias{global_prevalent}
\title{Identify prevalent and global types}
\usage{
global_prevalent(df)
}
\arguments{
\item{df}{data.frame, a table with one row for each sample. The table must
have a `serotype` variable.}
}
\value{
A data frame.
}
\description{
This function takes a data frame which contains one row for each sample and
return a data frame of prevalent and global serotypes.
}
