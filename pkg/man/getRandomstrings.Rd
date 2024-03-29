% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{getRandomstrings}
\alias{getRandomstrings}
\title{Generates random strings}
\usage{
getRandomstrings(center = NULL, maxdist = NULL, length = nchar(center),
  n = 1, method = "lv")
}
\arguments{
\item{center}{Reference string, i.e. the cluster center}

\item{maxdist}{The maximum allowed string distance}

\item{length}{The length of the string}

\item{n}{Number of strings to be generated}

\item{method}{The string distance method used to calculate the string, defaults to Levensthein distance}
}
\value{
A character string
}
\description{
Generates random strings
}
\examples{
getRandomstrings(center="hello", maxdist = 2, n = 5)
}

