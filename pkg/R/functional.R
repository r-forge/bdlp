# misc functions for functional data

timesamples <- function(ObsNr,  MaxTimePoints, MinTimePoints, int = c(Tstart,Tfin), NrGridPoints, reg.time = FALSE){
    T <- matrix(0,nrow = ObsNr,ncol = MaxTimePoints)
    t <- matrix(0,nrow = ObsNr,ncol = MaxTimePoints)
    N <- rep(0,ObsNr)
    seq <- seq(from = int[1],to=int[2], length.out = NrGridPoints)
    if(!reg.time) number <- sample(MinTimePoints:MaxTimePoints,ObsNr,replace=TRUE)
    else  number <- rep(MaxTimePoints, ObsNr)
    for(j in 1:ObsNr){
        if(!reg.time)  TimeIndex <- sort(sample(1:NrGridPoints,number[j]))
        else TimeIndex <-floor(seq(from= 1, to = NrGridPoints, by =(NrGridPoints-1)/(MaxTimePoints-1)))
        T[j,1:length(TimeIndex)] <- TimeIndex
        t[j,1:sum(T[j,]!=0)] <- seq[T[j,]]
    }
    N <- apply(T, 1, function(x) max(which(x>0)))
    return(list(TimeIndex=T,time=t, TimeNrIndivid=N))
}

samplegrid <- function(total_n,  minT, maxT, granularity, regular = FALSE){
    gridmat <- matrix(0, total_n, granularity)
    if(regular) points <- rep(maxT, total_n) else points <- sample(minT:maxT, total_n, replace=T)
    if(regular) positions <- sort(sample(1:granularity, maxT))
    for(i in 1:total_n){
      if(regular == F) positions <- sort(sample(1:granularity, points[i]))
      gridmat[i,positions] <- 1
    }    
    gridmat
}

formatFunc <- function(data){
    N <- table(data[,1])
    Obs <- matrix(0,length(N),max(N))
    for(i in 1:length(N)) {
        Obs[i,1:N[i]] <- data[,2][data[,1]==i]
     }
    T <-  matrix(0,length(N),max(N))
    for(i in 1:length(N)) T[i,1:N[i]] <- data[,3][data[,1]==i]
    return(structure(list(Obs, T, N), names = c("Obs", "T", "N")))
}

