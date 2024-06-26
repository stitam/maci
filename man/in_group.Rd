% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/in_group.R
\name{in_group}
\alias{in_group}
\title{Group samples using text data}
\usage{
in_group(x, include = character(), exclude = character(), cache = NULL)
}
\arguments{
\item{x}{character; a vector of strings used for grouping samples.}

\item{include}{character; a vector of inclusion terms which clearly identify
samples that belong to the group.}

\item{exclude}{character; a vector of exclusion terms which clearly identify
samples that do not belong to the group.}
}
\value{
A logical vector. Elements of the vector can be \code{TRUE} if an 
inclusion term was matched, \code{FALSE} if an exclusion term was matched,
or \code{NA} if neither inclusion nor exclusion terms were matched.
}
\description{
This function looks at inclusion and exclusion terms in a character vector to
determine whether samples belong to a group or not.
}
\examples{
\dontrun{
# Identify human associated samples using metadata terms from multiple
variables.

# Step 1. Look at variable called "host"
include <- "human"
exclude <- "dog"
meta$human_related <- in_group(
  meta$host,
  include = include,
  exclude = exclude
)

# Step 2. Look at variable "isolation_source"
# Only group samples that were not grouped in the previous step.
include <- "human blood"
exclude <- "wastewater"
meta$human_related <- in_group(
  meta$isolation_source,
  include = include,
  exclude = exclude,
  cache = meta$human_related)
)
}
}
