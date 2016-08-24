#' Generates a number of datasets from one metadata scenario
#'
#' @param name The path to the setup file
#' @param setnr The metadata scenario, as taken from the info table
#' @param draws The number of datasets that are drawn from the metadata scenario
#' @param seedinfo The random number generator seed parameters
#' @param metaseedinfo If necessary, a separate set of random number generator parameters for the metadata (e.g. cluster centers)
#' @param file A custom file name for the output database. Defaults to the pattern setupname_setnr_seed.db
#' @param increment The random number seed will increase by 1 for each draw from the base seed given in seedinfo
#' @return An SQLite database that contains the desired number of data sets drawn from a certain metadata scenario
# @examples
# create.dataset(name="dangl2014.R", setnr=1, draws=10)
#' @export
create.dataset <- function(name = NULL, setnr = NULL, draws = 1,  
                           seedinfo = list(100, 
                                           paste(R.version$major, R.version$minor, sep = "."),
                                           RNGkind()), 
                           metaseedinfo = list(100, 
                                           paste(R.version$major, R.version$minor, sep = "."),
                                           RNGkind()),
                           file = NULL, increment = T){

  if(!is.null(name)) source(name)
  spl <- strsplit(name, "/")
  spl <- unlist(spl)
  name <- spl[length(spl)]
  name <- unlist(strsplit(name, ".R"))
  
  pb <- txtProgressBar(min = 0, max = draws, style = 3)

  if(is.null(file)){ 
    dbname <- paste(name, "_set_", setnr, "_", "seed_", seedinfo[[1]] ,".db", sep="")
    file.create(dbname) 
  }
  else {
    dbname <- file
    file.create(dbname) 
  }
  
  for(i in 1:draws){
    if(increment == T) seedinfo[[1]] <- seedinfo[[1]] + 1
    metadata <- read.metadata(name = name, setnr = setnr, seedinfo = seedinfo, metaseedinfo = metaseedinfo)
    output <- generate.data(metadata)
    write.Database(output, dbname, i) 
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat("\n")
  cat(paste(draws, "version(s) of set no.", setnr, "of", name, "experimental setup generated.\n"))
  cat(paste("Base seed ", seedinfo[[1]]-draws, " was used and is included in filename.\n", sep=""))
}


#' Generate a dataset from a metadata object
#'
#' @param m A metadata object
#' @return A dataset as specified by the metadata object
#' @examples
#' m <- new("metadata.metric", 
#'          clusters = list(c1 = list(n = 25, mu = c(4,5), Sigma=diag(1,2)),
#'                          c2 = list(n = 25, mu = c(-1,-2), Sigma=diag(1,2))),
#'          dist = mvrnorm)
#' generate.data(m)
#' @export
setGeneric("generate.data", function(m) {standardGeneric("generate.data")})

#' Generate a dataset from a metadata object
#'
#' @param m A metadata object
#' @return A dataset as specified by the metadata object
#' @export
setMethod("generate.data", signature(m = "metadata.metric"),
function(m){
	
  set.seed(m@seedinfo[[1]])
  RNGversion(m@seedinfo[[2]])
  RNGkind(m@seedinfo[[3]][1], m@seedinfo[[3]][2])
  
  total_n <- sum(unlist(lapply(m@clusters, function(x) x$n)))
  k <- length(m@clusters)
  
  #zum Ermitteln der Variablenanzahl - noch zu Verbessern
  for(i in 1:length(m@clusters[[1]])){
    if(is.matrix(m@clusters[[1]][[i]]))
      vars <- ncol(m@clusters[[1]][[i]])
  }
  
  samp <- m@dist
  
  datamatrix <- matrix(0, nrow=total_n, ncol=vars)
  
  clus <- list()
  for(i in 1:k){
    clus[[i]] <- do.call(samp, m@clusters[[i]])
  }
  datamatrix <- do.call(rbind, clus)
  
  if(m@standardization == "NONE")
    as.data.frame(datamatrix)
  else {
    stdz <- get(m@standardization)
    as.data.frame(stdz(datamatrix))
  }
})

#' Generate a dataset from a metadata object
#'
#' @param m A metadata object
#' @return A dataset as specified by the metadata object
#' @export
setMethod("generate.data", signature(m = "metadata.functional"), 
function(m){
  
  set.seed(m@seedinfo[[1]])
  RNGversion(m@seedinfo[[2]])
  RNGkind(m@seedinfo[[3]][1], m@seedinfo[[3]][2])
  
  nf <- length(m@functions)
  gridMatrix <- m@gridMatrix
  xvalvector <- yvalvector <- cluster <-  vector()
  
  xvals <- seq(m@interval[1], m@interval[2], by=1/(m@resolution-1))
  pointsPerFunc <- rowSums(gridMatrix)
  curves <- rep(1:m@total_n, pointsPerFunc)
  sddist <- get(m@sd_distribution)
  
  for(i in 1:m@total_n) xvalvector <- c(xvalvector, xvals[gridMatrix[i,] == 1])
  
  for(i in 1:nf){
    upper <- floor(m@total_n*i/nf)
    for(j in (floor(m@total_n*(i-1)/nf) + 1):upper){
      yvals <- m@functions[[i]](xvals[gridMatrix[j,] == 1])
      yvalvector <- c(yvalvector, yvals + sddist(yvals, sd=m@sd))
      cluster[j] <- i
    }
  }
  
  data <- cbind(curves, xvalvector, yvalvector)
  as.data.frame(data)
})

#' Generate a dataset from a metadata object
#'
#' @param m A metadata object
#' @return A dataset as specified by the metadata object
#' @export
setMethod("generate.data", signature(m = "metadata.ordinal"),
function(m){
  
  set.seed(m@seedinfo[[1]])
  RNGversion(m@seedinfo[[2]])
  RNGkind(m@seedinfo[[3]][1], m@seedinfo[[3]][2])
  samp <- m@dist
  total_n <- sum(unlist(lapply(m@clusters, function(x) x$n)))
  k <- length(m@clusters)
  for(i in 1:length(m@clusters$c1)){
    if(is.matrix(m@clusters$c1[[i]]))
      vars <- ncol(m@clusters$c1[[i]])
  }
  
  datamatrix <- matrix(0, nrow=total_n, ncol=vars)
  
  clus <- list()
  for(i in 1:k){
    clus[[i]] <- do.call(samp, m@clusters[[i]])
  }
  datamatrix <- do.call(rbind, clus)
  
  as.data.frame(datamatrix)
})

#' Generate a dataset from a metadata object
#'
#' @param m A metadata object
#' @return A dataset as specified by the metadata object
#' @export
setMethod("generate.data", signature(m = "metadata.binary"),
function(m){
  
  set.seed(m@seedinfo[[1]])
  RNGversion(m@seedinfo[[2]])
  RNGkind(m@seedinfo[[3]][1], m@seedinfo[[3]][2])
  samp <- m@dist
  total_n <- sum(unlist(lapply(m@clusters, function(x) x$n)))
  k <- length(m@clusters)
  for(i in 1:length(m@clusters$c1)){
    if(is.matrix(m@clusters$c1[[i]]))
      vars <- ncol(m@clusters$c1[[i]])
  }
  
  datamatrix <- matrix(0, nrow=total_n, ncol=vars)
  
  clus <- list()
  for(i in 1:k){
    clus[[i]] <- do.call(samp, m@clusters[[i]])
  }
  datamatrix <- do.call(rbind, clus)
  
  as.data.frame(datamatrix)
})


#setMethod("generate.data", signature(m = "metadata.wordnet"),
#function(m){
#  set.seed(m@seedinfo[[1]])
#  RNGversion(m@seedinfo[[2]])
#  RNGkind(m@seedinfo[[3]][1], m@seedinfo[[3]][2])
#  
#  require("wordnet")
#  
#  clus <- list()
#  for(i in 1:length(m@clusters)) {
#     filter <- getTermFilter(m@filtertype, m@clusters[[i]]$word, m@case_ignore)
#     terms <- getIndexTerms(m@clusters[[i]]$wordtype, m@clusters[[i]]$n, filter)
#     clus[[i]] <- sapply(terms, getLemma)
#  }
#  
#  data.frame(wordlist = unlist(clus))
#})

#' Generate a dataset from a metadata object
#'
#' @param m A metadata object
#' @return A dataset as specified by the metadata object
#' @export
setMethod("generate.data", signature(m = "metadata.randomstring"),
function(m){
  set.seed(m@seedinfo[[1]])
  RNGversion(m@seedinfo[[2]])
  RNGkind(m@seedinfo[[3]][1], m@seedinfo[[3]][2])
  
  samp <- m@genfunc
  
  clus <- list()
  for(i in 1:length(m@clusters)) {
    clus[[i]] <- do.call(samp, m@clusters[[i]])
  }
  
  data.frame(stringlist = unlist(clus))
})

#' Generates random strings
#'
#' @param center Reference string, i.e. the cluster center
#' @param maxdist The maximum allowed string distance
#' @param length The length of the string
#' @param n Number of strings to be generated
#' @param method The string distance method used to calculate the string, defaults to Levensthein distance
#' @return A character string
#' @examples
#' get.randomstrings(center="hello", maxdist = 2, n = 5)
#' @export
get.randomstrings <- function(center = NULL, maxdist = NULL, length = nchar(center), n = 1, method = "lv"){
  l <- vector()
  l[1:n] <- ""
  for(i in 1:n){
    while(TRUE){
      str <- paste(sample(letters, length, T), collapse="")
      if(stringdist::stringdist(center, str, method) <= maxdist) {
        l[i] <- str
        break
      }
    }
  }
  return(l)
}

write.Database <- function(output, dbname, draw){
  driver <- dbDriver("SQLite")
  con <- dbConnect(driver, dbname)
  if(draw == 1)
    dbWriteTable(con, paste("draw", draw, sep="_"), output, row.names=FALSE)
  else
    dbWriteTable(con, paste("draw", draw, sep="_"), output, append=T, row.names=FALSE)
}


