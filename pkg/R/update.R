## modify metadata objects

update.metadata <- function(m, ...){
  args <- c(as.list(environment()), list(...))
  objname <- args$m
  a <- names(args)
  
  args <- args[!a == "m"]
  a <- a[!grepl("m", a)]
  
  if(all(a %in% slotNames(m))){
    for(i in 1:length(args)){
      slot(m, a[i]) <- args[[i]] 
    }
  } else {
    wrong <- a[!(a %in% slotNames(m))]
    print(paste("Argument", wrong, "not defined in the metadata object!", sep= " "))
  }
  m
}


