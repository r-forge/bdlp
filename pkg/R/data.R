## main function to create a dataset

create.dataset <- function(name, setnr = NULL, draws = 1,  
                           seedinfo = list(100, 
                                           paste(R.version$major, R.version$minor, sep = "."),
                                           RNGkind()), 
                           metaseedinfo = list(100, 
                                           paste(R.version$major, R.version$minor, sep = "."),
                                           RNGkind()),
                           newsetup = NULL, file = NULL, increment = T){

  #if(!is.null(rng)) RNGkind(rng[1], rng[2])
  #if(!is.null(rv)) RNGversion(rv)
  if(!is.null(newsetup)) source(newsetup)
  
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
    writeToDatabase(output, dbname, i) 
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat("\n")
  cat(paste(draws, "version(s) of set no.", setnr, "of", name, "experimental setup generated.\n"))
  cat(paste("Base seed ", seedinfo[[1]]-draws, " was used and is included in filename.\n", sep=""))
}

## actual data generators

setGeneric("generate.data", function(m) standardGeneric("generate.data"))

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
  
  samp <- get(m@dist)
  
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

setMethod("generate.data", signature(m = "metadata.ordinal"),
function(m){

  require("GenOrd")
  
  set.seed(m@seedinfo[[1]])
  RNGversion(m@seedinfo[[2]])
  RNGkind(m@seedinfo[[3]][1], m@seedinfo[[3]][2])
  samp <- get(m@dist)
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


setMethod("generate.data", signature(m = "metadata.binary"),
function(m){

  require("MultiOrd")
  
  set.seed(m@seedinfo[[1]])
  RNGversion(m@seedinfo[[2]])
  RNGkind(m@seedinfo[[3]][1], m@seedinfo[[3]][2])
  samp <- get(m@dist)
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


setMethod("generate.data", signature(m = "metadata.wordnet"),
function(m){
  set.seed(m@seedinfo[[1]])
  RNGversion(m@seedinfo[[2]])
  RNGkind(m@seedinfo[[3]][1], m@seedinfo[[3]][2])
  
  require("wordnet")
  
  clus <- list()
  for(i in 1:length(m@clusters)) {
     filter <- getTermFilter(m@filtertype, m@clusters[[i]]$word, m@case_ignore)
     terms <- getIndexTerms(m@clusters[[i]]$wordtype, m@clusters[[i]]$n, filter)
     clus[[i]] <- sapply(terms, getLemma)
  }
  
  data.frame(wordlist = unlist(clus))
})


setMethod("generate.data", signature(m = "metadata.randomstring"),
function(m){
  set.seed(m@seedinfo[[1]])
  RNGversion(m@seedinfo[[2]])
  RNGkind(m@seedinfo[[3]][1], m@seedinfo[[3]][2])
  
  samp <- get(m@genfunc)
  
  require("stringdist")
  
  clus <- list()
  for(i in 1:length(m@clusters)) {
    clus[[i]] <- do.call(samp, m@clusters[[i]])
  }
  
  data.frame(stringlist = unlist(clus))
})



get_randomstrings <- function(center, maxdist, length = 10, n = 10, method = "lv"){
  l <- vector()
  l[1:n] <- ""
  for(i in 1:n){
    while(TRUE){
      str <- paste(sample(letters, length, T), collapse="")
      if(stringdist(center, str, method) <= maxdist) {
        l[i] <- str
        break
      }
    }
  }
  return(l)
}

