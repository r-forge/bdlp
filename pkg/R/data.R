## main function to create a dataset

create.dataset <- function(name, setnr = NULL, draws = 1, seedVersion = NULL, seedKind = NULL, seed, newsetup = NULL, file = NULL){

  if(!is.null(seedKind)) RNGkind(seedKind)
  if(!is.null(seedVersion)) RNGversion(seedVersion)
  if(!is.null(newsetup)) source(newsetup)
  pb <- txtProgressBar(min = 0, max = draws, style = 3)

  if(is.null(file)){ 
    dbname <- paste(name, "_set_", setnr, "_", "seed_", seed ,".db", sep="")
    file.create(dbname) 
  }
  else {
    dbname <- file
    file.create(dbname) 
  }
  
  for(i in 1:draws){
    seed <- seed + 1
    metadata <- read.metadata(name = name, setnr = setnr, seed = seed)
    output <- generate.data(metadata)
    writeToDatabase(output, dbname, i) 
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat("\n")
  cat(paste(draws, "version(s) of set no.", setnr, "of", name, "experimental setup generated.\n"))
  cat(paste("Base seed ", seed-draws, " was used and is included in filename.\n", sep=""))
}

## actual data generators

setGeneric("generate.data", function(m) standardGeneric("generate.data"))

setMethod("generate.data", signature(m = "metadata.metric"),
function(m){
  
  set.seed(m@seed)
  total_n <- sum(unlist(lapply(m@clusters, function(x) x$n)))
  k <- length(m@clusters)
  
  for(i in 1:length(clusters$cl1)){
    if(is.matrix(clusters$cl1[[i]]))
      vars <- ncol(clusters$cl1[[i]])
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
  
  set.seed(m@seed)
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
  samp <- get(m@dist)
  total_n <- sum(unlist(lapply(m@clusters, function(x) x$n)))
  k <- length(m@clusters)
  for(i in 1:length(clusters$cl1)){
    if(is.matrix(clusters$cl1[[i]]))
      vars <- ncol(clusters$cl1[[i]])
  }
  
  datamatrix <- matrix(0, nrow=total_n, ncol=vars)
  
  clus <- list()
  for(i in 1:k){
    clus[[i]] <- do.call(samp, m@clusters[[i]])
  }
  datamatrix <- do.call(rbind, clus)
  
  as.data.frame(datamatrix)
})




