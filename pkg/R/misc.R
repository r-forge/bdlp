# misc functions 

writeToDatabase <- function(output, dbname, draw){
  require(RSQLite)
  driver <- dbDriver("SQLite")
  con <- dbConnect(driver, dbname)
  if(draw == 1)
    dbWriteTable(con, paste("draw", draw, sep="_"), output, row.names=FALSE)
  else
    dbWriteTable(con, paste("draw", draw, sep="_"), output, append=T, row.names=FALSE)
}

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

  cat("Sourcing input file ... ")
  tryCatch({
    source(file)
    cat("Done.<br>")
  }, error = function(e) {
    suppressMessages(stop())
    cat("Failed!<br>")
  })
  ##-----------------------
  
  cat("Checking consistency of function names ... ")
  tryCatch({
    name <- strsplit(file, "/")
    name <- name[[1]][length(name[[1]])]
    name <- strsplit(name, ".R")[[1]]
    stopifnot(
      is.function(get(name)),
      regexpr(pattern = "[a-z]+[0-9]{4}", text=name) == T
    )
    cat("Done.<br>")
  }, error = function(e) {
    suppressMessages(stop())
    cat("Failed.<br>")
  })
  
  ##-----------------------
  
  cat("Checking whether a summary is produced ... ")
  tryCatch({
    func <- get(name)
    s <- setupsummary(name)
    stopifnot(is.list(s))
    cat("Done.<br>")
  }, error = function(e) {
    cat("Failed.<br>")
  })
  
  ##-----------------------
  
  cat("Checking whether metadata and datasets can be generated ... ")
  tryCatch({
    func <- get(name)
    s <- setupsummary(name)
    runs <- nrow(s$summary)  
    
    
    for(i in 1:runs){
      meta <- func(setnr = i, seed = 1)
      data <- generate.data(meta)
    }
    cat("Done.<br>")
    #return(T)
  }, 
    error = function(e) {
    suppressMessages(stop())
    cat("Failed!<br>")
  })
  ##-----------------------
  cat("<p>Validity check successful ... upload completed. Setup available in library!")
  updateLibrary(name)
}


updateLibrary <- function(name){
  func <- get(name)
  m <- func(1,1)
  type <- strsplit(class(m)[1], "metadata.")[[1]][2]
  lib <- read.csv("data/library.csv")
  newlib <- rbind(lib, data.frame(type=type, name=name))
  write.csv(newlib, file = "data/library.csv",row.names=F)
}

setGeneric("plot", function(m, option) standardGeneric("plot"))
 
setMethod("plot", signature(m = "metadata.metric", option = "numeric"),
function(m, option){
  options(rgl.useNULL=TRUE)
  require(rgl)
  data <- generate.data(m)
  pr <- prcomp(data)
  prpred <- predict(pr, data)
  mems <- unlist(lapply(m@clusters, function(x) x$n))
  mems <- rep(1:m@k, mems)
  if(option == 1) plot.default(prpred[,1:2], col=mems) else plot3d(prpred[,1:3], col=mems)
})

setMethod("plot", signature(m = "metadata.functional"),
function(m){
  data <- generate.data(m)
  ppf <- rowSums(m@gridMatrix)
  plot.default(data[1:ppf[1],2:3], type="l", xlim=c(m@interval[1], m@interval[2]), ylim=c(min(data$yvalvector), max(data$yvalvector)))
  for(i in 1:m@total_n){
    lines(data$xvalvector[data$curves == i], data$yvalvector[data$curves == i], col = i)
  }
})

create.fileskeleton <- function(newname, mail, inst, author, type, mat, cit){

  #filename <- paste(newname,".R", sep="")

  #file.create(filename)
  
  d <- vector("character", 10)
  
  d[1] <- "# Simulation setup created by"
  
  d[3] <- paste("# ", author, sep="")
  d[4] <- paste("# ", inst, sep="")
  d[5] <- paste("# ", mail, sep="")
  
  d[7] <- paste(newname, " <- function(setnr=NULL, seed=NULL, info=F){", sep="")
  
  l <- length(capture.output(dput(mat)))
  
  d[9] <- "  # setup info table"
  
  d[11:(11+l-1)] <- capture.output(dput(mat))
  
  d[11] <- paste("infotable <- ", d[11], sep="") 
  for(i in 1:l) d[11+i-1] <- paste("  ", d[11+i-1], sep="")
  
  
  d[11+l+1] <- paste("  reference <- ", "\"", as.character(cit), "\"", sep="")
  
  d[11+l+3] <- "  if(info==T) return(list(summary=infotable, reference=reference))"
  
  d[11+l+5] <- "  # calculate metadata for a given setnr and seed"
  d[11+l+6] <- "  # and create new metadata object"
  
  
  if(type=="metric")
    d[11+l+8] <- "  new(\"metadata.metric\", ...)"
  if(type=="functional")
    d[11+l+8] <- "  new(\"metadata.functional\", ...)"
  if(type=="ordinal")
    d[11+l+8] <- "  new(\"metadata.ordinal\", ...)"

  d[(11+l+10):(11+l+16)] <- ""
  
  d[11+l+10] <- "}"
  d[11+l+11] <- ""
  d[11+l+12] <- "#------------------------------------------------------------"
  d[11+l+13] <- ""
  d[11+l+14] <- "# add additional function that are needed"
  d[11+l+15] <- "# during metadata generation here"
  d[11+l+16] <- "# appropriate function name: "
  
  d[11+l+16] <- paste(d[11+l+16], newname, "_myfunc()", sep="")
  
  d[is.na(d)] <- ""
  
  #con <- file(filename)
  #writeLines(d,con)
  #close(con)
  d
}

read.metadata <- function(name, setnr, seed){
    do.call(name, list(setnr=setnr, seed=seed))
}

setupsummary <- function(name) {
  d <- do.call(name, list(info=T))
  rows <- nrow(d$summary)
  d$summary <- cbind(setnr = 1:rows, d$summary)
  d
}



