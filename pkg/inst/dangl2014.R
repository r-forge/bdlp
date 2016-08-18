##testsetup
require(MASS)
dangl2014 <- function(setnr = NULL, 
                      seedinfo = list(100, 
                                      paste(R.version$major, R.version$minor, sep = "."),
                                      RNGkind()), 
                      info = FALSE, 
                      metaseedinfo = list(100, 
                                          paste(R.version$major, R.version$minor, sep = "."),
                                          RNGkind())){

inf <- data.frame(n = c(50, 40), k = c(2,2), shape = c("spherical", "spherical"))
ref <- "Dangl R. (2014) A small simulation study. Journal of Simple Datasets 10(2), 1-10"
if(info == T) return(list(summary = inf, reference = ref))

if(is.null(metaseedinfo)) metaseedinfo <- seedinfo

set.seed(metaseedinfo[[1]])
RNGversion(metaseedinfo[[2]])
RNGkind(metaseedinfo[[3]][1], metaseedinfo[[3]][2])

if(setnr == 1) {
  return(new("metadata.metric", 
    clusters = list(c1 = list(n = 25, mu = c(4,5), Sigma=diag(1,2)),
                    c2 = list(n = 25, mu = c(-1,-2), Sigma=diag(1,2))),
    dist = "mvrnorm", seedinfo = seedinfo))
}
if(setnr == 2){
  return(new("metadata.metric", 
    clusters = list(c1 = list(n = 20, mu = c(0,2), Sigma=diag(1,2)),
                    c2 = list(n = 20, mu = c(-1,-2), Sigma=diag(1,2))),
    dist = "mvrnorm", seedinfo = seedinfo))
}

}
