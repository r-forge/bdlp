# misc functions 

list.setups <- function(){
  lib <- read.csv("data/library.csv")
  for(i in 1:nrow(lib)){
    source(paste("data/", lib[i,2], ".R", sep=""))
  }
  lib
}

searchfunc <- function(author, year, keyword){

  setupnames <- as.character(list.setups()[,2])
  results <- list()
  citationlist <- list()

  search <- c(author, year, keyword)
  args <- which(search != "")
  
  for(i in 1:length(setupnames)){
   func <- get(setupnames[i])
   citationlist[[i]] <- func(info=T)$reference 
  }
  
  if(length(args) == 1){
    for(i in 1:length(citationlist)){
      if(grepl(search[args], citationlist[[i]], ignore.case=T)){
        results[[i]] <- citationlist[[i]] 
      }
    }
  }
  if(length(args) == 2){
    for(i in 1:length(citationlist)){
      if(grepl(search[args[1]], citationlist[[i]], ignore.case=T) && grepl(search[args[2]], citationlist[[i]], ignore.case=T)){
        results[[i]] <- citationlist[[i]] 
      }
    }    
  }
  if(length(args) == 3){
    for(i in 1:length(citationlist)){
      if(grepl(search[args[1]], citationlist[[i]], ignore.case=T) && grepl(search[args[2]], citationlist[[i]], ignore.case=T) && grepl(search[args[3]], citationlist[[i]], ignore.case=T)){
        results[[i]] <- citationlist[[i]] 
      }
    } 
  }
  data.frame(Results=unlist(results))
}

check.setup <- function(file){

  cat("Sourcing input file ... \n")
  tryCatch({source(file); cat("Done.\n")}, 
           warning = function(w) {stop("File not found!\n")}, 
           error = function(e) {stop("Unsuccessful. Check stopped.\n")}
           )
  ##-----------------------
  
  cat("Checking consistency of function names ... \n")
  tryCatch({
    name <- strsplit(file, "/")
    name <- name[[1]][length(name[[1]])]
    name <- strsplit(name, ".R")[[1]]
    stopifnot(
      is.function(get(name)),
      regexpr(pattern = "[a-z]+[0-9]{4}", text=name) == T
    )
    cat("Done.\n")}, 
    warning = function(w) {stop("Inconsistencies found.\n")},
    error = function(e) {stop("Inconsistencies detected. Check stopped.\n")})
  
  ##-----------------------
  
  cat("Checking reference ... \n")
  tryCatch({
    unc <- get(name)
    s <- setupsummary(name)
    stopifnot(is.character(s$reference))
    cat("Done.\n")}, 
    warning = function(w) {stop("No reference.\n")},
    error = function(e) {stop("No reference found. Check stopped.\n")})
  
  ##-----------------------
  
  cat("Checking whether a summary is produced ... \n")
  tryCatch({
    func <- get(name)
    s <- setupsummary(name)
    stopifnot(is.data.frame(s$summary))
    cat("Done.\n")
  }, 
  warning = function(w) {stop("No summary.\n")},
  error = function(e) {stop("No summary output available. Check stopped.\n")})
  
  ##-----------------------
  
  cat("Checking whether metadata and datasets can be generated ... \n")
  tryCatch({
    func <- get(name)
    s <- setupsummary(name)
    runs <- nrow(s$summary)  
    
    
    for(i in 1:runs){
      meta <- func(setnr = i)
      data <- generate.data(meta)
    }
    cat("Done.\n")
  }, 
    warning = function(w) {stop("No data.\n")},
    error = function(e) {stop("Data generation failed. Check stopped.\n")})
  ##-----------------------
  
  cat("Check complete! You can upload your benchmarking setup!\n")
}


updateLibrary <- function(name){
  func <- get(name)
  m <- func(1)
  type <- strsplit(class(m)[1], "metadata.")[[1]][2]
  lib <- read.csv("data/library.csv")
  newlib <- rbind(lib, data.frame(type=type, name=name))
  write.csv(newlib, file = "data/library.csv",row.names=F)
}

#' Plot a metadata object
#' 
#' @param m A metadata object
#' @param option Only for metric data - 3d instead of 2d plot
#' @return A plot, created by generating an instance of the dataset from the metadata object
#' @examples
#' m <- dangl2014(1)
#' plot.metadata(m)
#' @export
setGeneric("plot.metadata", function(m, option) {standardGeneric("plot.metadata")})

setMethod("plot.metadata", signature(m = "metadata.metric", option = "character"),
function(m, option = "2d"){
  options(rgl.useNULL=TRUE)
  require(rgl)
  data <- generate.data(m)
  mems <- unlist(lapply(m@clusters, function(x) x$n))
  mems <- rep(1:length(m@clusters), mems)
  pr <- prcomp(data)
  prpred <- predict(pr, data)
  if(ncol(data) == 2) {
    plot.default(data, col= mems)
  } else if(ncol(data) == 3 && option == "3d"){
	plot3d(prpred[,1:3], col=mems)
  } else{
    plot.default(prpred[,1:2], col=mems) 
  }
})

setMethod("plot.metadata", signature(m = "metadata.functional"),
function(m){
  data <- generate.data(m)
  ppf <- rowSums(m@gridMatrix)
  plot.default(data[1:ppf[1],2:3], type="l", xlim=c(m@interval[1], m@interval[2]), ylim=c(min(data$yvalvector), max(data$yvalvector)))
  for(i in 1:m@total_n){
    lines(data$xvalvector[data$curves == i], data$yvalvector[data$curves == i], col = i)
  }
})

setMethod("plot.metadata", signature(m = "metadata.ordinal"),
function(m){
  data <- generate.data(m)
  k <- length(m@clusters)
  
  n <- vector()
  
  for(i in 1:length(m@clusters)){
    for(j in 1:length(m@clusters[[i]])){
      if(length(m@clusters[[i]][[j]]) == 1)
        n[i] <- m@clusters[[i]][[j]] 
    }
  }
  
  cu <- cumsum(n)
  
  if(k < 5) form <- c(2,2)
  if(k > 4 && k < 10) form <- c(3,3)
  if(k > 9 && k < 17) form <- c(4,4)
  
  
  par(mfrow=c(form))
  for(i in 1:k){
    if(i == 1) {
		mat <- sapply(data[1:cu[i],], function(x) table(x))
	} else {
		mat <- sapply(data[(cu[i-1]+1):(cu[i]),], function(x) table(x))
	}
    p <- prop.table(mat, margin=2)
    barplot(p, col = rainbow(nrow(mat)))
  }
})

setMethod("plot.metadata", signature(m = "metadata.binary"),
function(m){
  data <- generate.data(m)
  k <- length(m@clusters)
  
  n <- vector()
  
  for(i in 1:length(m@clusters)){
    for(j in 1:length(m@clusters[[i]])){
      if(length(m@clusters[[i]][[j]]) == 1)
        n[i] <- m@clusters[[i]][[j]] 
    }
  }
  
  cu <- cumsum(n)
  
  if(k < 5) form <- c(2,2)
  if(k > 4 && k < 10) form <- c(3,3)
  if(k > 9 && k < 17) form <- c(4,4)
  
  
  par(mfrow=c(form))
  for(i in 1:k){
    if(i == 1){
      mat <- sapply(data[1:cu[i],], function(x) table(x))
    } else {
      mat <- sapply(data[(cu[i-1]+1):(cu[i]),], function(x) table(x))
    }
    p <- prop.table(mat, margin=2)
    barplot(p, col = rainbow(2))
  }
})


create.fileskeleton <- function(newname, mail, inst, author, type, mat, cit, codefile = TRUE){
  
  if(codefile == TRUE) {
    filename <- paste(newname,".R", sep="")
    file.create(filename)
  }
  
  d <- vector("character", 10)
  
  d[1] <- "# Simulation setup created by"
  
  d[3] <- paste("# ", author, sep="")
  d[4] <- paste("# ", inst, sep="")
  d[5] <- paste("# ", mail, sep="")
  
  d[7] <- paste(newname, " <- function(setnr = NULL,
                              seedinfo = list(100, 
                                paste(R.version$major, R.version$minor, sep = '.'), RNGkind()), 
                              info = FALSE, 
                              metaseedinfo = list(100, 
                                paste(R.version$major, R.version$minor, sep = '.'),
                                RNGkind())){", sep="")
  
  l <- length(capture.output(dput(mat)))
  
  d[9] <- "  # setup info table"
  
  d[11:(11+l-1)] <- capture.output(dput(mat))
  
  d[11] <- paste("infotable <- ", d[11], sep="") 
  for(i in 1:l) d[11+i-1] <- paste("  ", d[11+i-1], sep="")
  
  
  d[11+l+1] <- paste("  reference <- ", "\"", as.character(cit), "\"", sep="")
  
  d[11+l+3] <- "  if(info==T) return(list(summary=infotable, reference=reference))"
  
  d[11+l+5] <- "  # calculate metadata for a given setnr and seed"
  d[11+l+6] <- "  # and create new metadata object"
  d[11+l+7] <- ""
  
  d[11+l+8] <-  "  set.seed(metaseedinfo[[1]])"
  d[11+l+9] <-  "  RNGversion(metaseedinfo[[2]])"
  d[11+l+10] <- "  RNGkind(metaseedinfo[[3]][1], metaseedinfo[[3]][2])"
  
  d[11+l+11] <- ""
  
  if(type=="metric")
    d[11+l+12] <- "  new(\"metadata.metric\", ...)"
  if(type=="functional")
    d[11+l+12] <- "  new(\"metadata.functional\", ...)"
  if(type=="ordinal")
    d[11+l+12] <- "  new(\"metadata.ordinal\", ...)"
  if(type=="binary")
    d[11+l+12] <- "  new(\"metadata.binary\", ...)"
  if(type=="randomstring")
    d[11+l+12] <- "  new(\"metadata.randomstring\", ...)"
  
  d[11+l+13] <- ""


  
  d[11+l+14] <- "}"
  d[11+l+15] <- ""
  d[11+l+16] <- "#------------------------------------------------------------"
  d[11+l+17] <- ""
  d[11+l+18] <- "# add additional function that are needed"
  d[11+l+19] <- "# during metadata generation here"
  d[11+l+20] <- "# appropriate function name: "
  
  d[11+l+20] <- paste(d[11+l+20], newname, "_myfunc()", sep="")
  
  d[is.na(d)] <- ""
  
  if(codefile == T) {
    con <- file(filename)
    writeLines(d,con)
    close(con)
  } else {
    return(d)
  }
}

read.metadata <- function(name, setnr, seedinfo = NULL, metaseedinfo = NULL){
    if(is.null(seedinfo) && is.null(metaseedinfo))
      do.call(name, list(setnr=setnr))
    else if(is.null(seedinfo) && !is.null(metaseedinfo))
      do.call(name, list(setnr=setnr, metaseedinfo=metaseedinfo))
    else if(!is.null(seedinfo) && is.null(metaseedinfo))
      do.call(name, list(setnr=setnr, seedinfo=seedinfo))
    else
      do.call(name, list(setnr=setnr, seedinfo=seedinfo, metaseedinfo=metaseedinfo))
}

#' Returns the setup summary
#'
#' @param name The name of the setup
#' @return The summary table of \code{name}
#' @examples
#' setupsummary(dangl2014)
#' @export
setupsummary <- function(name) {
  d <- do.call(name, list(info=T))
  rows <- nrow(d$summary)
  d$summary <- cbind(setnr = 1:rows, d$summary)
  d
}
