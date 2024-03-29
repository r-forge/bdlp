% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.R
\name{saveSetup}
\alias{saveSetup}
\title{Saves a list of metadata objects to a new setup file}
\usage{
saveSetup(name, author, mail, inst, cit = "Unpublished", objects, table,
  seedinfo = list(100, paste(R.version$major, R.version$minor, sep = "."),
  RNGkind()), metaseedinfo = list(100, paste(R.version$major, R.version$minor,
  sep = "."), RNGkind()), custom_funcs = NULL, custom_name = NULL)
}
\arguments{
\item{name}{The name of the new setup (and thus also the filename)}

\item{author}{Full name of the author}

\item{mail}{Contact e-mail address of the author}

\item{inst}{Institution of the author}

\item{cit}{Reference to the publication where the setup was used, defaults to unpublished}

\item{objects}{List of metadata objects}

\item{table}{Info table for the setup}

\item{seedinfo}{Random number generator parameters for the data sets}

\item{metaseedinfo}{Random number generator parameters for the metadata}

\item{custom_funcs}{Custom functions that are needed to generate the meta(data)}

\item{custom_name}{Custom filename that deviates from the authorYear format}
}
\value{
A .R file that can be processed by \code{create.dataset}
}
\description{
Saves a list of metadata objects to a new setup file
}
\examples{
require(MASS)
a = new("metadata.metric", 
        clusters = list(c1 = list(n = 25, mu = c(4,5), Sigma=diag(1,2)),
                        c2 = list(n = 25, mu = c(-1,-2), Sigma=diag(1,2))),
                        genfunc = mvrnorm)
b = new("metadata.metric", 
        clusters = list(c1 = list(n = 44, mu = c(1,2), Sigma=diag(1,2)),
                        c2 = list(n = 66, mu = c(-5,-6), Sigma=diag(1,2))),
                        genfunc = mvrnorm)
saveSetup(name="test2002.R", author="Mister Twister", mail="mister.twister@edu.com",
           inst="Twister University", cit="Simple Data, pp. 23-24", objects=list(a, b),
           table=data.frame(n = c(50, 110), k = c(2,2), shape = c("spherical", "spherical")))  
}

