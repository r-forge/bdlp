# metadata object definitions

setClass("metadata.common",
         representation(clusters = "list",
                        dist = "character",
                        seedinfo = "list"),
#                        vars = "numeric",
#                        total_n = "numeric",
#                        k = "numeric"),
         validity = function(object){
                      retval <- NULL
                      #checking whethere there are exactly as many entries in 'clusters' than the number in 'k'
#                      if(length(object@clusters) != object@k)
#                        retval <- "No. of list entries in slot 'clusters' and no. of clusters in slot 'k' do not match"
                      #checking whether k does not exceed the number of observations
#                      if(object@total_n < object@k)
#                        retval <- "Fewer observations than clusters"
                      #checking whether all lists stored in 'cluster' are really lists
#                      checkval <- sum(unlist(lapply(object@clusters, function(x) is.list(x))))
#                      if(checkval != object@k)
#                        retval <- "Number of list entries in slot 'clusters' and number of actual cluster do not match"
                      #checking whether all cluster metadata arguments are valid arguments for the data generating function
                      for(i in 1:length(object@clusters)) {
                        if(!all(names(object@clusters[[i]]) %in% names(formals(get(object@dist))))){
                          retval <- "Distribution function arguments and cluster metadata names do not match" 
                        }
                      }
                      #checking whether the number of metadata arguments is consistent over all clusters
                      argvals <- unlist(lapply(object@clusters, function(x) length(x)))
                      if(mean(argvals) != median(argvals))
                        retval <- "Inconsistencies in slot 'clusters' regarding number of arguments"
                      retval
                    },
         prototype = prototype(seedinfo = list(100, 
                                          paste(R.version$major, R.version$minor, sep = "."),
                                          RNGkind())
                              ))

setClass("metadata.metric",
         contains = "metadata.common",
         representation = representation(standardization = "character"),
         validity = function(object){
                      retval <- NULL
                    #  if(length(object@clusters) != object@k) 
                    #    retval <- "Cluster metadata and no. of clusters do not match"
                    #  for(i in 1:object@k) {
                    #    if(!all(names(object@clusters[[i]]) %in% names(formals(get(object@dist))))){
                    #      retval <- "Distribution function arguments and cluster metadata names do not match" 
                    #    }
                    #  }
                    #  if(!all(unlist(lapply(object@clusters, function(x) length(x$mu) == object@vars))))
                    #    retval <- "Dimensionality mismatch between slot 'vars' and at least one cluster center"
                    
                    #checking whether cluster center dimensionality of each cluster metadata list matches the overall number
#                      argnum <- sum(unlist(lapply(object@clusters, function(x) lapply(x, function(y) length(y) == object@vars))))
#                      if(argnum != object@k)
#                        retval <- "Dimensionality mismatch between slot 'vars' and at least one cluster center"
                    retval
                    },
         prototype = prototype(standardization = "NONE",
                   dist = "mvrnorm"))
           
setClass("metadata.functional",
         representation(functions = "list", 
                        minTimePoints = "numeric",
                        maxTimePoints = "numeric",
                        total_n = "numeric",
                        resolution = "numeric",
                        interval = "numeric",
                        sd = "numeric",
                        gridMatrix = "matrix",
                        seed = "numeric",
                        sd_distribution = "character",
                        regular = "logical"),
          validity = function(object){
                       retval = NULL
                       if(object@maxTimePoints <= object@minTimePoints)
                         retval <- "Number of maximum time points must be greater than the minimum number"
                       for(i in 1:length(object@functions)){
                         testvector <- seq(object@interval[1], object@interval[2], length.out=object@resolution)
                         FUN <- object@functions
                         vals <- FUN[[i]](testvector)
                         if(any(is.na(vals)) == T)
                           retval <- paste("Function ", i, " returned NA or NaN on a possible grid point")
                       }
                       retval
                     },
          prototype = prototype(sd_distribution = "rnorm",
                                sd = 0.2,
                                seed = 100,
                                minTimePoints = 2,
                                maxTimePoints = 10,
                                resolution = 100,
                                regular=T))


setClass("metadata.ordinal",
         contains = "metadata.common",
         validity = function(object){
                      retval = NULL
                      if(!all(unlist(lapply(clusters, function(x) names(x) %in% names(formals(get(object@dist)))))))
                        retval <- "Nonconforming arguments found in slot 'clusters'"
                      retval
         },
         prototype = prototype(dist = "ordsample"))
         
         
setClass("metadata.binary",
         contains = "metadata.common",
         validity = function(object){
			          retval = NULL
			          if(!all(unlist(lapply(clusters, function(x) names(x) %in% names(formals(get(object@dist)))))))
                      retval <- "Nonconforming arguments found in slot 'clusters'"
                      retval
	     },
		 prototype = prototype(dist = "rmvbin"))
		 

#setClass("metadata.mixed",
#         contains = "metadata.common",
#         validity = function(object){
#			 
#			 
#		 },
#		 prototype = prototype(

