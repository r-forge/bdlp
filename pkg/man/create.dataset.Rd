\name{create.dataset}
\alias{create.dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a dataset by supplying a metadata file
}
\description{
This function creates one or more instances of a dataset 
defined by a metadata source file that either has been downloaded 
form BDLP or is supplied by the user.
}
\usage{
create.dataset(name, setnr = NULL, draws = 1, seedVersion = NULL, seedKind = NULL, seed, newsetup = NULL, file = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
the name of the setup that is supplied, is downloaded automatically from the BDLP
}
  \item{setnr}{
the number of the data set within the setup that is supplied. Can be 
retrieved by checking online at the BDLP or by running setupsummary(name)
}
  \item{draws}{
the number of instances that should be created
}
  \item{seedVersion}{
can be set if an older version of R should be used for the random number 
generator
}
  \item{seedKind}{
can be set to change the random number generator
}
  \item{seed}{
sets the random number seed
}
  \item{newsetup}{
needs to be supplied if 'name' is not supplied to indicate a setup that has been created 
locally or has been downloaded manually by the user
}
  \item{file}{
string that is to be used as a file name. The default file name is 'name'_set_'setnr'_seed_'seed'.db
}
}
\details{
Not really many details necessary at the moment.
}
\value{
The function returns an SQLite Database containing the requested number of datasets.
}
\references{
http://www.boku.ac.at/iasc/bdlp
}
\author{
Rainer Dangl
}
\note{
No notes
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (name, setnr = NULL, draws = 1, seedVersion = NULL, 
    seedKind = NULL, seed, newsetup = NULL, file = NULL) 
{
    if (!is.null(seedKind)) 
        RNGkind(seedKind)
    if (!is.null(seedVersion)) 
        RNGversion(seedVersion)
    if (!is.null(newsetup)) 
        source(newsetup)
    pb <- txtProgressBar(min = 0, max = draws, style = 3)
    if (is.null(file)) {
        dbname <- paste(name, "_set_", setnr, "_", "seed_", seed, 
            ".db", sep = "")
        file.create(dbname)
    }
    else {
        dbname <- file
        file.create(dbname)
    }
    for (i in 1:draws) {
        seed <- seed + 1
        metadata <- read.metadata(name = name, setnr = setnr, 
            seed = seed)
        output <- generate.data(metadata)
        writeToDatabase(output, dbname, i)
        Sys.sleep(0.1)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    cat("\n")
    cat(paste(draws, "version(s) of set no.", setnr, "of", name, 
        "experimental setup generated.\n"))
    cat(paste("Base seed ", seed - draws, " was used and is included in filename.\n", 
        sep = ""))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
